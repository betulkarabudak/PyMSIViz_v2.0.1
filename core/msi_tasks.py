import numpy as np
from core.utils import getionimage, bisect_spectrum
from core.volcano_analytics import (
    build_global_spectrum_for_features,
    select_volcano_features)
from pyimzml.ImzMLParser import ImzMLParser
try:
    from scipy.stats import (
        gaussian_kde, norm,
        mannwhitneyu,
        kruskal, f_oneway, levene, shapiro
    )
    _HAS_SCIPY_STATS = True
except Exception:
    _HAS_SCIPY_STATS = False
def task_compute_ionimage( parser, mz, tol, z, reduce_func=sum, *, progress_cb=None, cancel_cb=None):
    tol = abs(float(tol))
    mz = float(mz)
    z = int(z)

    tol = abs(float(tol))
    mz = float(mz)
    z = int(z)

    max_y = parser.imzmldict["max count of pixels y"]
    max_x = parser.imzmldict["max count of pixels x"]
    im = np.zeros((max_y, max_x), dtype=float)

    coords = parser.coordinates
    total = len(coords)

    if total == 0:
        return {
            "image": im,
            "mz": mz,
            "tol": tol,
            "z": z,
            "cancelled": False
        }

    step = max(1, total // 100)

    for i, (x, y, z_) in enumerate(coords):
        if cancel_cb and cancel_cb():
            return {
                "image": im,
                "mz": mz,
                "tol": tol,
                "z": z,
                "cancelled": True
            }

        if z_ == z:
            mzs, ints = map(np.asarray, parser.getspectrum(i))
            min_i, max_i = bisect_spectrum(mzs, mz, tol)
            if max_i >= min_i:
                im[y - 1, x - 1] = reduce_func(ints[min_i:max_i + 1])

        if progress_cb and (i % step == 0):
            progress_cb(int((i / total) * 100))

    if progress_cb:
        progress_cb(100)

    return {
        "image": im,
        "mz": mz,
        "tol": tol,
        "z": z,
        "cancelled": False
    }

def task_ai_assist_top_peak(
    parser,
    *,
    max_pixels=1200,
    n_bins=6000
):
    """
    Global mean spectrum üzerinden en yüksek piki bulur.
    SADECE hesap yapar, UI YOK.
    """

    mz, inten = build_global_spectrum_for_features(
        parser,
        max_pixels=max_pixels,
        n_bins=n_bins
    )

    if mz is None or mz.size == 0:
        return None

    return float(mz[int(np.argmax(inten))])

def task_roi_compare_violin(
    parser,
    roi_mgr,
    selected,
    mz1, tol1,
    mz2, tol2,
    mz3, tol3,
    use_multi,
    *,
    max_pixels=None,
    progress_cb=None,
    cancel_cb=None
):
    data = []
    labels = []

    infoA = roi_mgr.get_roi(selected[0])
    infoB = roi_mgr.get_roi(selected[1])

    if infoA is None or infoB is None:
        return {"data": [], "labels": []}

    maskA, zA = infoA["mask"], int(infoA["z"])
    maskB, zB = infoB["mask"], int(infoB["z"])

    if zA != zB:
        return {"data": [], "labels": []}

    z = zA

    mz_list = [(mz1, tol1)]
    if use_multi:
        if mz2 is not None and tol2 is not None:
            mz_list.append((mz2, tol2))
        if mz3 is not None and tol3 is not None:
            mz_list.append((mz3, tol3))

    for i, (name, mask) in enumerate(
        [(selected[0], maskA), (selected[1], maskB)]
    ):
        if cancel_cb and cancel_cb():
            return {"data": [], "labels": []}

        vals_all = []

        for mz, tol in mz_list:
            if cancel_cb and cancel_cb():
                return {"data": [], "labels": []}

            img = getionimage(parser, float(mz), float(tol), z)

            if img is None or img.shape != mask.shape:
                continue

            vals = img[mask].astype(float)
            vals = vals[np.isfinite(vals)]

            if vals.size > 0:
                vals_all.append(vals)

        if vals_all:
            merged = np.concatenate(vals_all)
            data.append(merged)
            labels.append(name)

        if progress_cb:
            progress_cb(int(((i + 1) / 2) * 100))

    if len(data) != 2:
        return {"data": [], "labels": []}

    return {
        "data": data,
        "labels": labels
    }


def _volcano_roi_intensity_matrix(parser, features, tol, z, flat_indices, *, cancel_cb=None, progress_cb=None):
    """Integrate all selected m/z windows per spectrum in a single ROI-only scan.

    Returns a matrix (selected spatial pixels x features), retaining zero intensity
    for pixels absent from the imzML spatial coordinate list. Window endpoints are
    inclusive, matching getionimage. No full imaging data cube is allocated.
    """
    features = np.asarray(features, dtype=float)
    selected = np.asarray(flat_indices, dtype=np.int64)
    width = int(parser.imzmldict["max count of pixels x"])
    matrix = np.zeros((selected.size, features.size), dtype=float)
    location = {int(v): j for j, v in enumerate(selected)}
    tolerance = abs(float(tol))
    lo_mz = features - tolerance
    hi_mz = features + tolerance
    coords = parser.coordinates
    total = len(coords)
    for i, (x, y, plane) in enumerate(coords):
        if cancel_cb is not None and cancel_cb():
            return None
        if int(plane) != int(z):
            continue
        idx = (int(y) - 1) * width + int(x) - 1
        row = location.get(idx)
        if row is None:
            continue
        mz, intensity = parser.getspectrum(i)
        mz = np.asarray(mz)
        intensity = np.asarray(intensity, dtype=float)
        starts = np.searchsorted(mz, lo_mz, side="left")
        stops = np.searchsorted(mz, hi_mz, side="right")
        valid = stops > starts
        # Python sum matches the legacy ion-image reducer's operation ordering.
        matrix[row, valid] = [float(sum(intensity[a:b])) for a, b in
                              zip(starts[valid], stops[valid])]
        if progress_cb and (i % max(1, total // 100) == 0):
            progress_cb(int((i + 1) * 100 / total))
    if progress_cb:
        progress_cb(100)
    return matrix


def task_roi_volcano(
    parser,
    saved_rois,
    nameA,
    nameB,
    tol,
    max_pixels=1500,
    alpha=0.05,
    *,
    random_state=0,
    progress_cb=None,
    cancel_cb=None
):
    import numpy as np
    from scipy.stats import mannwhitneyu
    import logging

    logger = logging.getLogger(__name__)

    empty_payload = {
        "mz": np.array([]),
        "logFC": np.array([]),
        "p_raw": np.array([]),
        "q": np.array([]),
        "reason": ""
    }

    try:
        # ---------- ROI CHECK ----------
        infoA = saved_rois.get(nameA)
        infoB = saved_rois.get(nameB)
        if infoA is None or infoB is None:
            empty_payload["reason"] = "roi_not_found"
            return empty_payload

        maskA, zA = infoA["mask"], int(infoA["z"])
        maskB, zB = infoB["mask"], int(infoB["z"])

        # ROIs must refer to the same acquisition plane and full image grid.
        if zA != zB:
            empty_payload["reason"] = "roi_z_mismatch"
            return empty_payload
        expected_shape = (parser.imzmldict["max count of pixels y"],
                          parser.imzmldict["max count of pixels x"])
        if maskA.shape != expected_shape or maskB.shape != expected_shape:
            empty_payload["reason"] = "roi_mask_shape_mismatch"
            return empty_payload
        if max_pixels < 3:
            empty_payload["reason"] = "insufficient_max_pixels"
            return empty_payload

        # ---------- PIXEL INDEX (DOWNSAMPLING) ----------
        idxA = np.where(maskA.ravel())[0]
        idxB = np.where(maskB.ravel())[0]

        rng = np.random.default_rng(random_state)   # reproducible pixel subsampling
        if idxA.size > max_pixels:
            idxA = rng.choice(idxA, max_pixels, replace=False)
        if idxB.size > max_pixels:
            idxB = rng.choice(idxB, max_pixels, replace=False)

        # ---------- FEATURE SELECTION ----------
        mz_grid, inten = build_global_spectrum_for_features(
            parser,
            max_pixels=max_pixels,
            random_state=random_state
        )
        # ---------- FEATURE LIMIT ----------
        MAX_FEATURES = 120
        if hasattr(parser, "large_dataset_mode") and parser.large_dataset_mode:
            MAX_FEATURES = 10
            logger.info(
                "Large dataset mode ENABLED → feature limit set to 10 (Volcano)"
            )

        # the MAX_FEATURES most intense peaks (returned sorted by m/z)
        features = select_volcano_features(mz_grid, inten, max_features=MAX_FEATURES)

        if features is None or len(features) < 5:
            empty_payload["reason"] = "not_enough_features"
            return empty_payload

        logfc = []
        pvals = []
        valid_mz = []

        # ---------- SINGLE-PASS ROI SPECTRAL EXTRACTION ----------
        combined_indices = np.union1d(idxA, idxB)
        data = _volcano_roi_intensity_matrix(
            parser, features, tol, zA, combined_indices,
            cancel_cb=cancel_cb, progress_cb=progress_cb
        )
        if data is None:
            empty_payload["reason"] = "cancelled"
            return empty_payload
        rows_A = np.searchsorted(combined_indices, idxA)
        rows_B = np.searchsorted(combined_indices, idxB)
        logfc, pvals, valid_mz = [], [], []
        for j, mzv in enumerate(features):
            vA = data[rows_A, j]
            vB = data[rows_B, j]
            vA = vA[np.isfinite(vA)]
            vB = vB[np.isfinite(vB)]
            if vA.size < 3 or vB.size < 3:
                continue
            medA, medB = np.median(vA), np.median(vB)
            logfc.append(np.log2((medB + 1e-9) / (medA + 1e-9)))
            _, p = mannwhitneyu(vA, vB, alternative="two-sided")
            pvals.append(p)
            valid_mz.append(mzv)

        # ---------- FINAL CHECK ----------
        if len(pvals) < 3:
            empty_payload["reason"] = "not_enough_valid_features"
            return empty_payload

        logfc = np.asarray(logfc)
        pvals = np.asarray(pvals)

        # ---------- BH-FDR ----------
        order = np.argsort(pvals)
        ranked = pvals[order]
        m = len(ranked)

        qvals = np.empty_like(ranked)
        prev = 1.0
        for i in reversed(range(m)):
            q = ranked[i] * m / (i + 1)
            q = min(q, prev)
            qvals[i] = q
            prev = q

        qvals_corrected = np.empty_like(qvals)
        qvals_corrected[order] = qvals


        return {
            "mz": np.asarray(valid_mz),
            "logFC": logfc,
            "p_raw": pvals,
            "q": qvals_corrected,
            "random_state": random_state,
            "max_pixels": max_pixels,
            "feature_selection": "global_peak_height_descending",
            "extraction_mode": "single_pass_selected_roi_pixels"
        }

    except Exception as e:
        empty_payload["reason"] = f"exception: {str(e)}"
        return empty_payload





def task_load_imzml_safe(filename):
    """
    SAFE loader for large imzML datasets.
    - No full spectrum traversal
    - No heavy RAM usage
    - Estimates m/z range from representative spectra
    """
    parser = ImzMLParser(filename)

    mz_min = None
    mz_max = None

    coords = parser.coordinates
    if coords:
        for i in range(min(10, len(coords))):
            mzs, _ = parser.getspectrum(i)
            if mzs is not None and len(mzs) > 0:
                mz_min = float(mzs[0])
                mz_max = float(mzs[-1])
                break

    if mz_min is None or mz_max is None:
        mz_min, mz_max = 0.0, 0.0

    info = {
        "mz_min": mz_min,
        "mz_max": mz_max
    }

    return parser, info
