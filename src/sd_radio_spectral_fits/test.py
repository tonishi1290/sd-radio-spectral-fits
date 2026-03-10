import os
import sd_radio_spectral_fits as sdrsf

# ==========================================
# CLIコマンドに基づいたパラメータ設定
# ==========================================
raw_file = "orikl_raw_dump_hot_off_on_1.fits"
cal_file = "orikl_tastar_1.fits"
coadd_file = "orikl_final_vlsrk_1.fits"
fit_file = "orikl_final_vlsrk_1_fitted.fits"

# 共有パラメータ
tamb_k = 300.0
vwin_list = ["-20:0", "17:40"]
poly_order = 3

def main():
    print("=== Full Pipeline Verification (CLI match) ===")

    if not os.path.exists(raw_file):
        print(f"Error: {raw_file} が見つかりません。")
        return

    # 1. sdrsf-make-tastar 相当
    print("\n[Step 1] Calibration (Ta*)...")
    sdrsf.run_tastar_calibration(
        input_data=raw_file,
        output_path=cal_file,
        tamb_k=tamb_k,
        overwrite=True
    )

    # 2. sdrsf-coadd-fits-vlsrk 相当
    # 積算時にベースライン補正とRMSウェイトを適用
    print("\n[Step 2] Velocity Coadd (LSRK + RMS Weighting)...")
    sdrsf.run_velocity_coadd(
        inputs=[cal_file],
        output_path=coadd_file,
        mode="rms_weight",
        vmin=-20,
        vmax=40,
        pos_tol_arcsec=1,
        baseline_vwin=vwin_list,
        baseline_poly=poly_order,
        rms_vwin=["-20:0"], # CLIの指定に合わせる
        rms_poly=3,
        overwrite=True
    )

    # 3. sdrsf-baseline-fit 相当
    # 積算済みのデータに対して再度フィッティング（BSL_COEF保存の確認用）
    print("\n[Step 3] Final Baseline Fitting (for BSL_COEF)...")
    sdrsf.run_baseline_fit(
        input_data=coadd_file,
        output_path=fit_file,
        vwin=vwin_list,
        poly_order=poly_order,
        overwrite=True
    )

    # 4. 可視化による最終チェック
    print("\n[Step 4] Final Verification with Viewer...")
    st = sdrsf.read_scantable(fit_file)
    
    if "BSL_COEF" in st.table.columns:
        print(f"✅ 'BSL_COEF' column verified. Sample: {st.table.iloc[0]['BSL_COEF']}")
        
        print("\nLaunching Interactive Viewer...")
        # 前回修正した plotting.py により、v_corr_kms がなくても 
        # CTYPE1='VRAD' や SPECSYS='LSRK' を見て正常に起動するはずです
        sdrsf.view_spectra(st, axis_type="vel")
    else:
        print("❌ Error: 'BSL_COEF' がテーブルに見当たりません。")

if __name__ == "__main__":
    main()
