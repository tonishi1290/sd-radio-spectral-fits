from __future__ import annotations

import copy
import dataclasses
import hashlib
import io
import json
from datetime import datetime, timezone
from functools import lru_cache
from importlib import metadata as importlib_metadata
from pathlib import Path
from typing import Any, Iterable, Sequence

import numpy as np
from astropy.io import fits
from astropy.table import Table

from .otf_bundle import OTFBundle

SCHEMA_VERSION = "otf-prov-v1"
PROV_EXTNAME = "PROVSTEP"
MAX_JSON_CHARS = 32760
MAX_TEXT_CHARS = 1024

_HEADER_KEYS = ("PRVSCHEM", "PRVNSTEP", "PRVMAIN", "PRVAUX", "PRVUTC", "BNDLID", "PRVCODE")


def _utc_now_iso() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")


@lru_cache(maxsize=1)
def _stable_code_version() -> str:
    for dist_name in ("sd-radio-spectral-fits", "sd_radio_spectral_fits"):
        try:
            ver = importlib_metadata.version(dist_name)
            ver = str(ver).strip()
            if ver:
                return f"sd-radio-spectral-fits={ver}"
        except importlib_metadata.PackageNotFoundError:
            pass
        except Exception:
            pass
    try:
        import sd_radio_spectral_fits as _pkg  # type: ignore

        ver = str(getattr(_pkg, "__version__", "") or "").strip()
        if ver:
            return f"sd-radio-spectral-fits={ver}"
    except Exception:
        pass
    return "unknown"


def _normalize_scalar(value: Any) -> Any:
    if isinstance(value, (np.bool_, bool)):
        return bool(value)
    if isinstance(value, (np.integer, int)):
        return int(value)
    if isinstance(value, (np.floating, float)):
        if np.isfinite(value):
            return float(value)
        if np.isnan(value):
            return "nan"
        if value > 0:
            return "inf"
        return "-inf"
    if isinstance(value, (np.str_, str)):
        return str(value)
    if value is None:
        return None
    return value


# Keep this deliberately conservative: provenance is for structured metadata,
# not for persisting bulky arrays.
def normalize_provenance_value(value: Any) -> Any:
    value = _normalize_scalar(value)
    if value is None or isinstance(value, (bool, int, float, str)):
        return value
    if dataclasses.is_dataclass(value) and not isinstance(value, type):
        return normalize_provenance_value(dataclasses.asdict(value))
    if isinstance(value, dict):
        return {str(k): normalize_provenance_value(v) for k, v in sorted(value.items(), key=lambda kv: str(kv[0]))}
    if isinstance(value, np.ndarray):
        arr = np.asarray(value)
        if arr.ndim == 0:
            return normalize_provenance_value(arr.item())
        if arr.size <= 64:
            return normalize_provenance_value(arr.tolist())
        return {
            "__ndarray__": True,
            "shape": [int(v) for v in arr.shape],
            "dtype": str(arr.dtype),
        }
    if isinstance(value, (list, tuple)):
        return [normalize_provenance_value(v) for v in value]
    if isinstance(value, set):
        return [normalize_provenance_value(v) for v in sorted(value, key=lambda x: repr(x))]
    if isinstance(value, fits.Header):
        return {str(k): normalize_provenance_value(value[k]) for k in value.keys()}
    if isinstance(value, Table):
        return {"__table__": True, "colnames": [str(c) for c in value.colnames], "nrow": int(len(value))}
    if hasattr(value, "tolist"):
        try:
            return normalize_provenance_value(value.tolist())
        except Exception:
            pass
    return repr(value)



def _summarize_unique_column(table: Any, name: str, *, max_values: int = 16) -> dict[str, Any] | None:
    try:
        if not hasattr(table, "columns") or name not in table.columns:
            return None
        series = table[name]
        values: list[Any] = []
        seen: set[str] = set()
        null_count = 0
        for raw in series:
            if raw is None:
                null_count += 1
                continue
            try:
                if np.isscalar(raw) and isinstance(raw, float) and not np.isfinite(raw):
                    null_count += 1
                    continue
            except Exception:
                pass
            norm = normalize_provenance_value(raw)
            key = json.dumps(norm, ensure_ascii=True, allow_nan=False, sort_keys=True, separators=(",", ":"))
            if key in seen:
                continue
            seen.add(key)
            values.append(norm)
            if len(values) >= max_values:
                break
        return {
            "present": True,
            "sample_values": values,
            "sample_truncated": len(seen) >= max_values,
            "null_count": int(null_count),
        }
    except Exception:
        return {"present": True, "sample_values": [], "sample_truncated": False, "error": "summary_failed"}



def summarize_scantable_input(scantables: Any, *, max_inputs: int = 16, max_unique_values: int = 16) -> dict[str, Any]:
    items = list(scantables) if isinstance(scantables, (list, tuple)) else [scantables]
    out: dict[str, Any] = {
        "n_input_scantables": int(len(items)),
        "per_input_truncated": False,
        "per_input": [],
        "total_nrow": 0,
        "nchan_values": [],
    }
    nchan_seen: list[int] = []
    key_meta = ("RESTFREQ", "RESTFRQ", "SPECSYS", "CTYPE1", "BUNIT", "TELESCOP", "OBJECT", "baseline_subtracted")
    col_keys = ("FDNUM", "IFNUM", "PLNUM", "OBSMODE", "SCAN", "IS_TURN")
    for idx, sc in enumerate(items):
        table = getattr(sc, "table", None)
        data = getattr(sc, "data", None)
        meta = getattr(sc, "meta", {}) or {}
        entry: dict[str, Any] = {"input_index": int(idx)}
        try:
            nrow = int(len(table)) if table is not None else None
        except Exception:
            nrow = None
        try:
            arr = np.asarray(data) if data is not None else None
            nchan = int(arr.shape[1]) if (arr is not None and arr.ndim == 2) else None
        except Exception:
            nchan = None
        entry["nrow"] = nrow
        entry["nchan"] = nchan
        if nrow is not None:
            out["total_nrow"] = int(out["total_nrow"]) + int(nrow)
        if nchan is not None and nchan not in nchan_seen:
            nchan_seen.append(int(nchan))
        entry["meta"] = {}
        for key in key_meta:
            if key in meta:
                entry["meta"][key] = normalize_provenance_value(meta.get(key))
        entry["columns"] = {}
        for key in col_keys:
            col_summary = _summarize_unique_column(table, key, max_values=max_unique_values)
            if col_summary is None:
                entry["columns"][key] = {"present": False}
            else:
                entry["columns"][key] = col_summary
        if idx < max_inputs:
            out["per_input"].append(entry)
        else:
            out["per_input_truncated"] = True
    out["nchan_values"] = sorted(int(v) for v in nchan_seen)
    return normalize_provenance_value(out)



def provenance_json_dumps(value: Any) -> str:
    txt = json.dumps(
        normalize_provenance_value(value),
        ensure_ascii=True,
        allow_nan=False,
        sort_keys=True,
        separators=(",", ":"),
    )
    if len(txt) > MAX_JSON_CHARS:
        raise ValueError(f"Provenance JSON payload is too long ({len(txt)} chars > {MAX_JSON_CHARS}).")
    return txt



def _clean_text(value: Any, *, field: str, max_chars: int = MAX_TEXT_CHARS) -> str:
    txt = "" if value is None else str(value)
    if len(txt) > max_chars:
        raise ValueError(f"Provenance text field {field!r} is too long ({len(txt)} chars > {max_chars}).")
    try:
        txt.encode("ascii")
    except UnicodeEncodeError as exc:
        raise ValueError(f"Provenance text field {field!r} must be ASCII-safe.") from exc
    return txt



def _input_closure(input_bundles: Sequence[OTFBundle]) -> tuple[list[dict[str, Any]], list[str]]:
    merged: list[dict[str, Any]] = []
    seen: set[str] = set()
    direct_input_step_ids: list[str] = []
    for bundle in input_bundles:
        prov = list(bundle.meta.get("provenance", []) or [])
        last_id = ""
        for rec in prov:
            sid = str(rec.get("step_id", "") or "")
            if sid:
                last_id = sid
            if sid and sid not in seen:
                merged.append(rec)
                seen.add(sid)
        direct_input_step_ids.append(last_id)
    return merged, direct_input_step_ids



def _step_hash_payload(*, op_id: str, kind: str, aux: str, params_input: Any, params_config: Any, params_resolved: Any, input_step_ids: Sequence[str]) -> str:
    payload = {
        "schema_version": SCHEMA_VERSION,
        "op_id": str(op_id),
        "kind": str(kind),
        "aux": str(aux),
        "params_input": normalize_provenance_value(params_input),
        "params_config": normalize_provenance_value(params_config),
        "params_resolved": normalize_provenance_value(params_resolved),
        "input_step_ids": [str(v) for v in input_step_ids],
    }
    return provenance_json_dumps(payload)



def _compute_step_id(*, op_id: str, kind: str, aux: str, params_input: Any, params_config: Any, params_resolved: Any, input_step_ids: Sequence[str]) -> str:
    payload = _step_hash_payload(
        op_id=op_id,
        kind=kind,
        aux=aux,
        params_input=params_input,
        params_config=params_config,
        params_resolved=params_resolved,
        input_step_ids=input_step_ids,
    )
    return hashlib.sha256(payload.encode("ascii")).hexdigest()



def derive_bundle_id(records: Sequence[dict[str, Any]]) -> str | None:
    if not records:
        return None
    sid = str(records[-1].get("step_id", "") or "")
    if not sid:
        return None
    txt = f"b:{sid[:32]}"
    return _clean_text(txt, field="bundle_id", max_chars=34)



def _safe_int(value: Any) -> int | None:
    try:
        if value is None or value == "":
            return None
        return int(value)
    except Exception:
        return None



def _header_summary_from_header(header: fits.Header | None) -> dict[str, Any]:
    if header is None:
        return {
            "schema": "",
            "nstep": None,
            "last_main_op": "",
            "last_aux_op": "",
            "last_utc": "",
            "bundle_id": "",
            "code_version": "",
        }
    return {
        "schema": str(header.get("PRVSCHEM", "") or ""),
        "nstep": _safe_int(header.get("PRVNSTEP", None)),
        "last_main_op": str(header.get("PRVMAIN", "") or ""),
        "last_aux_op": str(header.get("PRVAUX", "") or ""),
        "last_utc": str(header.get("PRVUTC", "") or ""),
        "bundle_id": str(header.get("BNDLID", "") or ""),
        "code_version": str(header.get("PRVCODE", "") or ""),
    }



def update_header_provenance_summary(header: fits.Header, records: Sequence[dict[str, Any]], *, bundle_id: str | None = None) -> None:
    for key in _HEADER_KEYS:
        if key in header:
            del header[key]
    if not records:
        return
    last = records[-1]
    last_main = None
    last_aux = None
    for rec in records:
        if str(rec.get("kind", "main")) == "aux":
            last_aux = rec
        else:
            last_main = rec
    header["PRVSCHEM"] = (SCHEMA_VERSION, "OTF provenance schema")
    header["PRVNSTEP"] = (int(len(records)), "Number of provenance steps")
    if last_main is not None:
        header["PRVMAIN"] = (_clean_text(last_main.get("op_id", ""), field="PRVMAIN", max_chars=68), "Last main provenance op")
    if last_aux is not None:
        header["PRVAUX"] = (_clean_text(last_aux.get("op_id", ""), field="PRVAUX", max_chars=68), "Last aux provenance op")
    header["PRVUTC"] = (_clean_text(last.get("utc", ""), field="PRVUTC", max_chars=68), "UTC of last provenance step")
    codever = _clean_text(last.get("code_version", ""), field="PRVCODE", max_chars=68)
    if codever:
        header["PRVCODE"] = (codever, "Code version for last step")
    if bundle_id:
        header["BNDLID"] = (_clean_text(bundle_id, field="BNDLID", max_chars=68), "Bundle provenance identifier")



def append_bundle_provenance_step(
    bundle: OTFBundle,
    *,
    input_bundles: Sequence[OTFBundle] | None,
    op_id: str,
    module: str,
    function: str,
    kind: str = "main",
    aux: str = "",
    params_input: Any = None,
    params_config: Any = None,
    params_resolved: Any = None,
    results_summary: Any = None,
    duration_sec: float | None = None,
    utc: str | None = None,
    code_version: str | None = None,
) -> dict[str, Any]:
    if input_bundles is None:
        base_records = list(bundle.meta.get("provenance", []) or [])
        direct_input_step_ids = [str(base_records[-1].get("step_id", ""))] if base_records else []
    else:
        base_records, direct_input_step_ids = _input_closure(list(input_bundles))
    record = {
        "step_id": _compute_step_id(
            op_id=op_id,
            kind=kind,
            aux=aux,
            params_input=params_input,
            params_config=params_config,
            params_resolved=params_resolved,
            input_step_ids=direct_input_step_ids,
        ),
        "op_id": _clean_text(op_id, field="op_id"),
        "kind": _clean_text(kind, field="kind", max_chars=32),
        "aux": _clean_text(aux, field="aux", max_chars=64),
        "utc": _clean_text(utc or _utc_now_iso(), field="utc", max_chars=32),
        "module": _clean_text(module, field="module"),
        "function": _clean_text(function, field="function"),
        "input_step_ids": [str(v) for v in direct_input_step_ids],
        "params_input": normalize_provenance_value(params_input),
        "params_config": normalize_provenance_value(params_config),
        "params_resolved": normalize_provenance_value(params_resolved),
        "results_summary": normalize_provenance_value(results_summary),
        "duration_sec": (None if duration_sec is None else float(duration_sec)),
        "code_version": _clean_text(code_version or _stable_code_version(), field="code_version", max_chars=128),
    }
    records = list(base_records)
    if not records or str(records[-1].get("step_id", "")) != record["step_id"]:
        records.append(record)
    bundle.meta["provenance_schema"] = SCHEMA_VERSION
    bundle.meta["provenance"] = records
    bundle_id = derive_bundle_id(records)
    if bundle_id is not None:
        bundle.meta["bundle_id"] = bundle_id
    update_header_provenance_summary(bundle.header, records, bundle_id=bundle.meta.get("bundle_id"))
    return record



def _stringify_cell(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, bytes):
        return value.decode("ascii")
    return str(value)



def provenance_records_to_table(records: Sequence[dict[str, Any]]) -> Table:
    rows: list[dict[str, Any]] = []
    for rec in records:
        rows.append({
            "STEPID": _clean_text(rec.get("step_id", ""), field="STEPID", max_chars=128),
            "OPID": _clean_text(rec.get("op_id", ""), field="OPID", max_chars=128),
            "KIND": _clean_text(rec.get("kind", ""), field="KIND", max_chars=32),
            "AUX": _clean_text(rec.get("aux", ""), field="AUX", max_chars=64),
            "UTC": _clean_text(rec.get("utc", ""), field="UTC", max_chars=32),
            "MODULE": _clean_text(rec.get("module", ""), field="MODULE", max_chars=128),
            "FUNC": _clean_text(rec.get("function", ""), field="FUNC", max_chars=128),
            "INSTEP": provenance_json_dumps(list(rec.get("input_step_ids", []) or [])),
            "PINPUT": provenance_json_dumps(rec.get("params_input", {})),
            "PCONFIG": provenance_json_dumps(rec.get("params_config", {})),
            "PRESOLV": provenance_json_dumps(rec.get("params_resolved", {})),
            "RESULT": provenance_json_dumps(rec.get("results_summary", {})),
            "CODEVER": _clean_text(rec.get("code_version", ""), field="CODEVER", max_chars=128),
            "DSEC": (np.nan if rec.get("duration_sec", None) is None else float(rec.get("duration_sec"))),
        })
    return Table(rows=rows)



def provenance_table_to_records(table: Table) -> list[dict[str, Any]]:
    if len(table) == 0:
        return []
    records: list[dict[str, Any]] = []
    names = {str(n).upper(): n for n in table.colnames}
    for row in table:
        def col(name: str, default: Any = "") -> Any:
            key = names.get(name)
            return default if key is None else row[key]

        rec = {
            "step_id": _stringify_cell(col("STEPID", "")),
            "op_id": _stringify_cell(col("OPID", "")),
            "kind": _stringify_cell(col("KIND", "main")) or "main",
            "aux": _stringify_cell(col("AUX", "")),
            "utc": _stringify_cell(col("UTC", "")),
            "module": _stringify_cell(col("MODULE", "")),
            "function": _stringify_cell(col("FUNC", "")),
            "input_step_ids": json.loads(_stringify_cell(col("INSTEP", "[]")) or "[]"),
            "params_input": json.loads(_stringify_cell(col("PINPUT", "{}")) or "{}"),
            "params_config": json.loads(_stringify_cell(col("PCONFIG", "{}")) or "{}"),
            "params_resolved": json.loads(_stringify_cell(col("PRESOLV", "{}")) or "{}"),
            "results_summary": json.loads(_stringify_cell(col("RESULT", "{}")) or "{}"),
            "code_version": _stringify_cell(col("CODEVER", "")),
            "duration_sec": None,
        }
        dsec = col("DSEC", np.nan)
        try:
            dsec_f = float(dsec)
            if np.isfinite(dsec_f):
                rec["duration_sec"] = dsec_f
        except Exception:
            pass
        records.append(rec)
    return records



def get_bundle_provenance(bundle: OTFBundle, *, copy_records: bool = True) -> list[dict[str, Any]]:
    records = list(bundle.meta.get("provenance", []) or [])
    return copy.deepcopy(records) if copy_records else records



def get_bundle_provenance_summary(bundle: OTFBundle) -> dict[str, Any]:
    summary = _header_summary_from_header(bundle.header)
    if "provenance_schema" in bundle.meta:
        summary["schema"] = str(bundle.meta.get("provenance_schema", "") or summary.get("schema", ""))
    if "bundle_id" in bundle.meta:
        summary["bundle_id"] = str(bundle.meta.get("bundle_id", "") or summary.get("bundle_id", ""))
    records = list(bundle.meta.get("provenance", []) or [])
    if records:
        summary["nstep"] = int(len(records))
        last = records[-1]
        summary["last_utc"] = str(last.get("utc", "") or summary.get("last_utc", ""))
        main_ops = [str(r.get("op_id", "")) for r in records if str(r.get("kind", "main")) != "aux"]
        aux_ops = [str(r.get("op_id", "")) for r in records if str(r.get("kind", "main")) == "aux"]
        if main_ops:
            summary["last_main_op"] = main_ops[-1]
        if aux_ops:
            summary["last_aux_op"] = aux_ops[-1]
        codever = str(last.get("code_version", "") or "")
        if codever:
            summary["code_version"] = codever
    return summary



def _open_fits_source(source: str | Path | fits.HDUList):
    if isinstance(source, fits.HDUList):
        return source, False
    hdul = fits.open(source, memmap=True)
    return hdul, True



def read_fits_provenance(source: str | Path | fits.HDUList, *, copy_records: bool = True) -> list[dict[str, Any]]:
    hdul, should_close = _open_fits_source(source)
    try:
        records: list[dict[str, Any]] = []
        for hdu in hdul[1:]:
            if isinstance(hdu, fits.BinTableHDU) and str(hdu.name).upper() == PROV_EXTNAME:
                records = provenance_table_to_records(Table(hdu.data))
                break
        return copy.deepcopy(records) if copy_records else records
    finally:
        if should_close:
            hdul.close()



def read_fits_provenance_summary(source: str | Path | fits.HDUList | fits.Header) -> dict[str, Any]:
    if isinstance(source, fits.Header):
        return _header_summary_from_header(source)
    hdul, should_close = _open_fits_source(source)
    try:
        return _header_summary_from_header(hdul[0].header)
    finally:
        if should_close:
            hdul.close()



def get_provenance(source: OTFBundle | str | Path | fits.HDUList, *, copy_records: bool = True) -> list[dict[str, Any]]:
    if isinstance(source, OTFBundle):
        return get_bundle_provenance(source, copy_records=copy_records)
    return read_fits_provenance(source, copy_records=copy_records)



def get_provenance_summary(source: OTFBundle | str | Path | fits.HDUList | fits.Header) -> dict[str, Any]:
    if isinstance(source, OTFBundle):
        return get_bundle_provenance_summary(source)
    return read_fits_provenance_summary(source)



def extract_provenance_table(source: OTFBundle | str | Path | fits.HDUList) -> Table:
    return provenance_records_to_table(get_provenance(source, copy_records=True))



def _truncate_text(text: str, max_chars: int) -> str:
    if max_chars <= 0 or len(text) <= max_chars:
        return text
    if max_chars <= 3:
        return text[:max_chars]
    return text[: max_chars - 3] + "..."



def _format_json_block(value: Any, *, indent: int, max_chars: int) -> list[str]:
    txt = json.dumps(normalize_provenance_value(value), ensure_ascii=True, sort_keys=True, indent=2)
    lines = txt.splitlines() or [txt]
    prefix = " " * indent
    return [prefix + _truncate_text(line, max_chars) for line in lines]



def format_provenance(
    source: OTFBundle | str | Path | fits.HDUList | fits.Header | Sequence[dict[str, Any]],
    *,
    include_params: bool = False,
    include_results: bool = False,
    max_chars: int = 160,
) -> str:
    if isinstance(source, fits.Header):
        summary = get_provenance_summary(source)
        records: list[dict[str, Any]] = []
    elif isinstance(source, Sequence) and not isinstance(source, (str, bytes, bytearray, OTFBundle, fits.HDUList, Path)):
        records = copy.deepcopy(list(source))
        summary = {
            "schema": SCHEMA_VERSION,
            "nstep": len(records),
            "last_main_op": next((str(r.get("op_id", "")) for r in reversed(records) if str(r.get("kind", "main")) != "aux"), ""),
            "last_aux_op": next((str(r.get("op_id", "")) for r in reversed(records) if str(r.get("kind", "main")) == "aux"), ""),
            "last_utc": str(records[-1].get("utc", "")) if records else "",
            "bundle_id": "",
            "code_version": str(records[-1].get("code_version", "")) if records else "",
        }
    else:
        summary = get_provenance_summary(source)  # type: ignore[arg-type]
        records = get_provenance(source, copy_records=True) if not isinstance(source, fits.Header) else []

    lines = [
        f"schema={summary.get('schema', '') or ''}",
        f"bundle_id={summary.get('bundle_id', '') or ''}",
        f"nstep={summary.get('nstep', None)}",
        f"last_main_op={summary.get('last_main_op', '') or ''}",
        f"last_aux_op={summary.get('last_aux_op', '') or ''}",
        f"last_utc={summary.get('last_utc', '') or ''}",
        f"code_version={summary.get('code_version', '') or ''}",
    ]
    if not records:
        return "\n".join(lines)

    lines.append("steps:")
    for idx, rec in enumerate(records):
        op = str(rec.get("op_id", "") or "")
        kind = str(rec.get("kind", "main") or "main")
        aux = str(rec.get("aux", "") or "")
        utc = str(rec.get("utc", "") or "")
        module = str(rec.get("module", "") or "")
        func = str(rec.get("function", "") or "")
        step_id = str(rec.get("step_id", "") or "")
        input_ids = list(rec.get("input_step_ids", []) or [])
        dsec = rec.get("duration_sec", None)
        head = f"[{idx}] op_id={op} kind={kind}"
        if aux:
            head += f" aux={aux}"
        if utc:
            head += f" utc={utc}"
        lines.append(head)
        lines.append(f"  function={module}::{func}")
        lines.append(f"  step_id={step_id}")
        lines.append(f"  input_step_ids={input_ids}")
        if dsec is not None:
            lines.append(f"  duration_sec={dsec}")
        if include_params:
            lines.append("  params_input:")
            lines.extend(_format_json_block(rec.get("params_input", {}), indent=4, max_chars=max_chars))
            lines.append("  params_config:")
            lines.extend(_format_json_block(rec.get("params_config", {}), indent=4, max_chars=max_chars))
            lines.append("  params_resolved:")
            lines.extend(_format_json_block(rec.get("params_resolved", {}), indent=4, max_chars=max_chars))
        if include_results:
            lines.append("  results_summary:")
            lines.extend(_format_json_block(rec.get("results_summary", {}), indent=4, max_chars=max_chars))
    return "\n".join(lines)



def show_provenance(
    source: OTFBundle | str | Path | fits.HDUList | fits.Header | Sequence[dict[str, Any]],
    *,
    include_params: bool = False,
    include_results: bool = False,
    max_chars: int = 160,
    file: io.TextIOBase | None = None,
) -> str:
    txt = format_provenance(
        source,
        include_params=include_params,
        include_results=include_results,
        max_chars=max_chars,
    )
    if file is None:
        print(txt)
    else:
        file.write(txt)
        if not txt.endswith("\n"):
            file.write("\n")
    return txt
