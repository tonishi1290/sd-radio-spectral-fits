from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any
import copy

import numpy as np
from astropy.io import fits
from astropy.table import Table


@dataclass
class OTFBundle:
    """
    Standard in-memory cube bundle for OTF gridding / coadd / FFT-PLAIT.

    Parameters
    ----------
    data
        Science cube, shape = (nchan, ny, nx).
    header
        FITS header carrying WCS and science metadata for ``data``.
    variance
        Variance cube with shape broadcastable to ``data``.
    valid_mask
        Bool mask, shape = (nchan, ny, nx) or (ny, nx).
    support_mask
        Bool spatial support mask, shape = (ny, nx).
    unit
        Data unit. If omitted, inferred from header ``BUNIT``.
    family_label
        User-level family label such as ``'X'`` or ``'Y'``.
    image_ext
        Additional 2D/3D image-like extensions, keyed by EXTNAME.
    table_ext
        Additional table-like extensions, keyed by EXTNAME.
    meta
        Arbitrary runtime metadata not suitable for FITS header cards.
    """

    data: np.ndarray
    header: fits.Header
    variance: np.ndarray | None = None
    valid_mask: np.ndarray | None = None
    support_mask: np.ndarray | None = None
    unit: str | None = None
    family_label: str | None = None
    image_ext: dict[str, np.ndarray] = field(default_factory=dict)
    table_ext: dict[str, Table] = field(default_factory=dict)
    meta: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        self.data = np.asarray(self.data)
        if self.data.ndim != 3:
            raise ValueError(f"OTFBundle.data must be 3D (nchan, ny, nx); got shape={self.data.shape}")
        if not isinstance(self.header, fits.Header):
            self.header = fits.Header(self.header)
        if self.variance is not None:
            self.variance = np.asarray(self.variance)
        if self.valid_mask is not None:
            self.valid_mask = np.asarray(self.valid_mask, dtype=bool)
        if self.support_mask is not None:
            self.support_mask = np.asarray(self.support_mask, dtype=bool)
        if self.unit is None:
            self.unit = str(self.header.get("BUNIT", "")) or None
        if self.family_label is None:
            fam = self.header.get("FAMILY", self.meta.get("family_label"))
            self.family_label = None if fam is None else str(fam)
        if self.family_label is not None:
            self.meta.setdefault("family_label", str(self.family_label))
        for key, arr in list(self.image_ext.items()):
            self.image_ext[key] = np.asarray(arr)
        for key, tab in list(self.table_ext.items()):
            if isinstance(tab, Table):
                continue
            if isinstance(tab, dict):
                self.table_ext[key] = Table(tab)
            else:
                self.table_ext[key] = Table(tab)

    @property
    def shape(self) -> tuple[int, int, int]:
        return self.data.shape

    @property
    def nchan(self) -> int:
        return int(self.data.shape[0])

    @property
    def ny(self) -> int:
        return int(self.data.shape[1])

    @property
    def nx(self) -> int:
        return int(self.data.shape[2])

    def copy(self, deep: bool = True) -> "OTFBundle":
        if deep:
            return OTFBundle(
                data=np.array(self.data, copy=True),
                header=self.header.copy(),
                variance=None if self.variance is None else np.array(self.variance, copy=True),
                valid_mask=None if self.valid_mask is None else np.array(self.valid_mask, copy=True),
                support_mask=None if self.support_mask is None else np.array(self.support_mask, copy=True),
                unit=None if self.unit is None else str(self.unit),
                family_label=None if self.family_label is None else str(self.family_label),
                image_ext={k: np.array(v, copy=True) for k, v in self.image_ext.items()},
                table_ext={k: v.copy(copy_data=True) for k, v in self.table_ext.items()},
                meta=copy.deepcopy(self.meta),
            )
        return OTFBundle(
            data=self.data,
            header=self.header,
            variance=self.variance,
            valid_mask=self.valid_mask,
            support_mask=self.support_mask,
            unit=self.unit,
            family_label=self.family_label,
            image_ext=dict(self.image_ext),
            table_ext=dict(self.table_ext),
            meta=dict(self.meta),
        )
