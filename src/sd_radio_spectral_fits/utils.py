from __future__ import annotations
from dataclasses import dataclass
import builtins
from typing import Iterator, Sequence
import numpy as np
import pandas as pd

@dataclass(frozen=True)
class FailPolicy:
    mode: str = "exit"  # exit|skip|mask

def iter_blocks(n: int, block_size: int) -> Iterator[tuple[int,int]]:
    if block_size <= 0:
        yield (0, n)
        return
    i = 0
    while i < n:
        j = builtins.min(n, i + block_size)
        yield (i, j)
        i = j

def parse_windows(specs: Sequence[str]) -> list[tuple[float,float]]:
    out: list[tuple[float,float]] = []
    for s in specs:
        a, b = s.split(":")
        a = float(a); b = float(b)
        if b < a:
            a, b = b, a
        out.append((a, b))
    return out

def in_any_windows(v: np.ndarray, windows: Sequence[tuple[float,float]]) -> np.ndarray:
    m = np.zeros(v.shape, dtype=bool)
    for a, b in windows:
        m |= (v >= a) & (v <= b)
    return m

def subtract_windows(base: Sequence[tuple[float,float]], sub: Sequence[tuple[float,float]]) -> list[tuple[float,float]]:
    if not base:
        return []
    def merge(ws):
        ws = sorted(ws, key=lambda x: x[0])
        merged=[]
        for a,b in ws:
            if not merged or a > merged[-1][1]:
                merged.append([a,b])
            else:
                merged[-1][1] = builtins.max(merged[-1][1], b)
        return [(a,b) for a,b in merged]
    base_m = merge(base)
    sub_m = merge(sub)
    if not sub_m:
        return base_m
    out=[]
    for a,b in base_m:
        cur = [(a,b)]
        for c,d in sub_m:
            new=[]
            for x,y in cur:
                if d <= x or c >= y:
                    new.append((x,y))
                else:
                    if c > x:
                        new.append((x, builtins.min(c,y)))
                    if d < y:
                        new.append((builtins.max(d,x), y))
            cur = new
            if not cur:
                break
        out.extend([(x,y) for x,y in cur if y > x])
    return out

_ALLOWED = {"icrs", "j2000", "fk5", "galactic"}


def normalize_frame(name: str) -> str:
    n = (name or "").strip().lower()
    if n in ("j2000", "fk5"):
        return "icrs"
    return n

def validate_mapping_frame(meta: dict, mapping: pd.DataFrame) -> str:
    frame = normalize_frame(str(meta.get("coord_frame", "icrs")))
    if frame not in _ALLOWED:
        raise ValueError(f"Unsupported coord_frame={meta.get('coord_frame')!r}. Allowed: icrs/j2000/fk5/galactic.")
    cols = {str(c).upper() for c in mapping.columns}

    if frame == "galactic":
        if not ({"GLON", "GLAT"} <= cols):
             raise ValueError("coord_frame=galactic but mapping lacks GLON/GLAT (or lowercase aliases).")
    else:
        has_ra = any(k in cols for k in ("RA", "RA_DEG", "OBSRA"))
        has_dec = any(k in cols for k in ("DEC", "DEC_DEG", "OBSDEC"))
        if not (has_ra and has_dec):
            raise ValueError("coord_frame=icrs but mapping lacks RA/DEC (or standard aliases).")
            
    return frame
