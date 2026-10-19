import numpy as np
from scipy.signal import savgol_filter, find_peaks


# ============================================================
# Global spectrum for volcano feature selection
# ============================================================

def build_global_spectrum_for_features(
    parser,
    max_pixels=1500,
    n_bins=6000,
    random_state=0
):
    """
    Build a global mean spectrum used ONLY for feature selection.
    """
    n = len(parser.coordinates)
    if n == 0:
        raise ValueError("No pixels found.")

    idxs = (
        np.arange(n)
        if n <= max_pixels
        else np.random.default_rng(random_state).choice(n, max_pixels, replace=False)
    )

    mzs_list, int_list = [], []

    for idx in idxs:
        mzs, ints = parser.getspectrum(int(idx))
        if mzs is not None and len(mzs) > 0:
            mzs_list.append(np.asarray(mzs))
            int_list.append(np.asarray(ints))

    if not mzs_list:
        raise ValueError("No valid spectra found in dataset.")

    lengths = [len(m) for m in mzs_list]

    # fast path: identical grids
    if len(set(lengths)) == 1:
        mz = np.mean(np.vstack(mzs_list), axis=0)
        inten = np.mean(np.vstack(int_list), axis=0)
        return mz.astype(float), inten.astype(float)

    # fallback: bin to common grid
    all_mz = np.concatenate(mzs_list)
    all_int = np.concatenate(int_list)

    bins = np.linspace(all_mz.min(), all_mz.max(), n_bins + 1)
    dig = np.digitize(all_mz, bins)

    mz, inten = [], []
    for i in range(1, len(bins)):
        m = (dig == i)
        if m.any():
            mz.append(bins[i])
            inten.append(all_int[m].mean())

    mz = np.asarray(mz, dtype=float)
    inten = np.asarray(inten, dtype=float)

    if mz.size == 0:
        raise ValueError("Empty global spectrum.")

    return mz, inten


# ============================================================
# Volcano feature selection
# ============================================================

def select_volcano_features(mz, inten, max_features=800, min_dist_mz=0.05):
    """Return actual local maxima ranked by smoothed peak height, highest first.

    No artificial features are generated for flat or peak-free spectra.
    """
    mz = np.asarray(mz, dtype=float)
    inten = np.asarray(inten, dtype=float)
    valid = np.isfinite(mz) & np.isfinite(inten)
    mz, inten = mz[valid], inten[valid]
    if mz.size < 3 or max_features <= 0:
        return np.array([], dtype=float)
    order = np.argsort(mz)
    mz, inten = mz[order], inten[order]
    if mz.size >= 7:
        window = min(31, mz.size if mz.size % 2 else mz.size - 1)
        inten = savgol_filter(inten, window, 3, mode="interp")
    inten = np.maximum(inten, 0.0)
    if np.ptp(inten) <= 1e-12 or np.max(inten) <= 0:
        return np.array([], dtype=float)
    step = np.median(np.diff(mz))
    distance = max(1, int(round(min_dist_mz / max(step, 1e-12))))
    peaks, props = find_peaks(inten, distance=distance,
                              height=np.percentile(inten, 75),
                              prominence=max(1e-12, 0.01 * np.ptp(inten)))
    if not peaks.size:
        return np.array([], dtype=float)
    heights = props["peak_heights"]
    ranked = np.argsort(-heights, kind="stable")[:max_features]
    return mz[peaks[ranked]]
