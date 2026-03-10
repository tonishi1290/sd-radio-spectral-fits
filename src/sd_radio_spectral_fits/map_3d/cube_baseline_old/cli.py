# src/sd_radio_spectral_fits/cube_baseline/cli.py
import argparse
import logging
import numpy as np
from astropy.io import fits

from .session import BaselineSession
from .orchestrator import run_one_iteration
from ..map.cube_analysis import calculate_moment0, append_analysis_hdus_to_fits

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")

def build_v_axis_kms(header: fits.Header, nchan: int) -> np.ndarray:
    """FITSヘッダーから速度軸(v_axis)を再構築し、必ず km/s 単位に変換する"""
    crval = header.get('CRVAL3', 0.0)
    cdelt = header.get('CDELT3', 1.0)
    crpix = header.get('CRPIX3', 1.0)
    v_axis = crval + (np.arange(nchan) + 1 - crpix) * cdelt
    
    ctype3 = header.get("CTYPE3", "").upper()
    cunit3 = header.get("CUNIT3", "km/s").strip().lower().replace(" ", "")
    
    if cunit3 in ("m/s", "ms-1", "m/sec"):
        v_axis = v_axis / 1000.0
    elif cunit3 == "hz" or "FREQ" in ctype3:
        restfreq = header.get("RESTFREQ", header.get("RESTFRQ", 0.0))
        if restfreq <= 0:
            raise ValueError("CUNIT3=Hz or CTYPE3=FREQ, but RESTFREQ is missing. Cannot convert v_axis to km/s.")
        # Radio convention velocity: v = (rest - obs) / rest * c
        v_axis = (restfreq - v_axis) / restfreq * 299792.458
        
    return v_axis

def get_dv_kms(header: fits.Header) -> float:
    """ヘッダーの単位情報を解釈して km/s 単位の分解能を計算する"""
    ctype3 = header.get("CTYPE3", "").upper()
    cunit3 = header.get("CUNIT3", "km/s").strip().lower().replace(" ", "")
    dv = abs(header.get("CDELT3", 1.0))
    
    if cunit3 in ("m/s", "ms-1", "m/sec"):
        return dv / 1000.0
    elif cunit3 == "hz" or "FREQ" in ctype3:
        restfreq = header.get("RESTFREQ", header.get("RESTFRQ", 0.0))
        if restfreq <= 0:
            raise ValueError("CUNIT3=Hz or CTYPE3=FREQ, but RESTFREQ is missing. Cannot convert to km/s.")
        return (dv / restfreq) * 299792.458
    else:
        return dv

def run_cli_pipeline(
    input_fits: str, 
    output_fits: str, 
    iterations: int = 2, 
    poly_order: int = 1,
    method: str = 'derivative'
):
    logging.info(f"Loading FITS: {input_fits}")
    with fits.open(input_fits, memmap=True) as hdul:
        # np.array() で明示的にメモリに載せることで、withを抜けた後のmemmap起因のクラッシュを防ぐ
        cube_fits = np.array(hdul[0].data)
        header = hdul[0].header.copy()
        
    if cube_fits.ndim != 3:
        raise ValueError("Primary HDU must contain a 3D cube.")
        
    cube_original = np.transpose(cube_fits, (1, 2, 0))
    ny, nx, nchan = cube_original.shape
    
    logging.info(f"Cube shape: {cube_original.shape}")
    
    # v_axis を最初から km/s で生成することで、手動窓の単位と完全に整合させる
    v_axis = build_v_axis_kms(header, nchan)
    dv_kms = get_dv_kms(header)
    
    session = BaselineSession(cube_original, v_axis)
    
    for i in range(iterations):
        logging.info(f"--- Iteration {i+1} / {iterations} ---")
        run_one_iteration(
            session=session,
            auto_mask_method=method,
            auto_mask_kwargs={'deriv_snr': 3.0} if method == 'derivative' else {},
            manual_v_windows=None,
            poly_order=poly_order,
            target_mask_2d=None
        )
        
        current_rms = np.nanmedian(session.fit_stats['rms_map'])
        logging.info(f"  -> Median RMS: {current_rms:.4f}")

    logging.info("Extracting final results...")
    cube_work = session.get_full_cube_work()
    mask_3d = session.mask_3d_current
    
    mom0_map = calculate_moment0(cube_work, mask_3d, dv_kms=dv_kms)
    
    cube_out_fits = np.transpose(cube_work, (2, 0, 1))
    hdul_out = fits.HDUList([fits.PrimaryHDU(data=cube_out_fits.astype(np.float32), header=header)])
    
    append_analysis_hdus_to_fits(hdul_out, mask_3d, mom0_map, header)
    
    # RMSマップの安全な書き出し (3次元目のWCS行列を徹底パージ)
    rms_hdr = header.copy()
    rms_hdr['BTYPE'] = 'BaselineRMS'
    rms_hdr['WCSAXES'] = 2
    
    keys_to_remove = ["CTYPE3", "CUNIT3", "CRPIX3", "CRVAL3", "CDELT3", "NAXIS3", "CROTA3"]
    for k in list(rms_hdr.keys()):
        if (k.startswith("PC") or k.startswith("CD")) and "3" in k:
            keys_to_remove.append(k)
            
    for k in set(keys_to_remove):
        rms_hdr.pop(k, None)
        
    hdul_out.append(fits.ImageHDU(data=session.fit_stats['rms_map'].astype(np.float32), header=rms_hdr, name='RMS'))

    logging.info(f"Saving to {output_fits}")
    hdul_out.writeto(output_fits, overwrite=True)
    logging.info("Done!")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Cube Baseline Iterative CLI")
    parser.add_argument("input", help="Input 3D FITS file")
    parser.add_argument("output", help="Output 3D FITS file")
    parser.add_argument("--iter", type=int, default=2, help="Number of iterations")
    parser.add_argument("--order", type=int, default=1, help="Polynomial order")
    parser.add_argument("--method", type=str, default="derivative", help="Masking method")
    
    args = parser.parse_args()
    run_cli_pipeline(args.input, args.output, iterations=args.iter, poly_order=args.order, method=args.method)
