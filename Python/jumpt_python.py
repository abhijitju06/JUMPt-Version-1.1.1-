#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
JUMPt Python — Setting-2 Python Edition 



Run:
  python jumpt_python.py --params JUMPt_python.params
"""

import os, re, math, argparse, warnings
from dataclasses import dataclass
from typing import List, Tuple

import numpy as np
import pandas as pd
from scipy.optimize import least_squares, minimize
from scipy.linalg import expm
from scipy.stats import t as student_t

warnings.filterwarnings("ignore", category=DeprecationWarning)

LN2 = math.log(2.0)

# -----------------------------
# Param parsing
# -----------------------------
def parse_params(path: str) -> dict:
    P = {}
    with open(path, "r", encoding="utf-8") as f:
        for raw in f:
            s = raw.strip()
            if not s or "=" not in s or s.startswith(("#","%",";")):
                continue
            k, v = s.split("=", 1)
            k = k.strip()
            v = re.split(r"[#%]", v, 1)[0].strip().rstrip(";").strip().strip("'").strip('"')
            # numeric if possible
            try:
                if re.fullmatch(r"[-+]?\d+", v):
                    v = int(v)
                elif re.fullmatch(r"[-+]?(?:\d*\.\d+|\d+\.\d*)(?:[eE][-+]?\d+)?", v):
                    v = float(v)
            except:
                pass
            P[k] = v

    # sensible defaults 
    P.setdefault("number_of_timepoints", 3)
    P.setdefault("bin_size", 30)
    P.setdefault("purity_of_SILAC_food", 99.0)  # heavy % in food -> impurity = 1 - purity/100
    P.setdefault("alphaUB", 0.1)
    P.setdefault("hl_tmin_days", 0.03)
    P.setdefault("hl_tmax_days", 100.0)
    P.setdefault("hl_A_min_days", 1.0)
    P.setdefault("hl_A_max_days", 100.0)
    P.setdefault("disable_physics_filter", 0)
    P.setdefault("physics_filter_mode", "clip")  # or 'drop'
    P.setdefault("min_observations_per_protein", 3)
    P.setdefault("free_lys_weight_mode", "auto")  # 'auto' => wF = M; 'matlab' => wF = M/5; or numeric
    P.setdefault("initial_condition_mode", "ones")  # 'ones' or 'data'
    P.setdefault("apparent_T50_calculation", 1)
    P.setdefault("random_seed", 12345)
    P.setdefault("n_starts", 20)
    P.setdefault("max_nfev_global", 800)
    P.setdefault("max_nfev_refine", 400)
    P.setdefault("sort_by_apparent_t50", 1)  # 1=sort (default), 0=keep input order
    return P


# -----------------------------
# Data container
# -----------------------------
@dataclass
class DataStruct:
    t: np.ndarray
    LysRatio: np.ndarray           # (T,)
    ProtInfo: pd.DataFrame         # first 1–2 cols: Protein, Gene (if present)
    SILAC_data: np.ndarray         # (T, N)
    timepoints: List[str]
    protBins: List[int]
    Ansatz_: np.ndarray            # (n_starts, N+2) — kept for determinism/compat
    physics_filter_mode: str
    free_lys_weight_mode: str
    y0_mode: str


# -----------------------------------------
# Robust Excel loader (matches sheet style)
# -----------------------------------------
def _find_time_row(df: pd.DataFrame, tcols: List[str]) -> int:
    for r in range(min(5, len(df))):
        vals = pd.to_numeric(df.loc[r, tcols], errors="coerce").to_numpy(float)
        if np.all(np.isfinite(vals)):
            return r
    return 0

def _find_free_lys_row(df: pd.DataFrame, tcols: List[str], after_row: int) -> int:
    prefix_cols = df.columns[:df.columns.get_loc(tcols[0])]
    for r in range(after_row+1, min(after_row+10, len(df))):
        s = " ".join(str(df.loc[r, c]) for c in prefix_cols if pd.notna(df.loc[r, c]))
        s_low = s.lower()
        if ("free" in s_low and "lys" in s_low) or ("lysine" in s_low):
            return r
    return min(after_row+1, len(df)-1)

def _read_input_like_matlab(xlsx_path: str):
    All = pd.read_excel(xlsx_path, header=0, dtype=str)
    tcols = [c for c in All.columns if str(c).lower().startswith("time_point")]
    if not tcols:
        raise ValueError("No 'time_point*' columns found in the input sheet.")

    time_row = _find_time_row(All, tcols)
    t = pd.to_numeric(All.loc[time_row, tcols], errors="coerce").to_numpy(float)

    free_row = _find_free_lys_row(All, tcols, time_row)
    LysRatio = pd.to_numeric(All.loc[free_row, tcols], errors="coerce").to_numpy(float)

    first_tp_idx = All.columns.get_loc(tcols[0])
    ProtInfo = All.iloc[free_row+1:, :first_tp_idx].copy().reset_index(drop=True)

    SILAC_cells = All.iloc[free_row+1:, All.columns.get_indexer(tcols)]
    SILAC_data = SILAC_cells.apply(pd.to_numeric, errors="coerce").to_numpy(float).T

    # sanitize to [0,1] or NaN
    SILAC_data[(SILAC_data < 0) | (SILAC_data > 1)] = np.nan
    LysRatio[(LysRatio < 0) | (LysRatio > 1)] = np.nan

    return t, LysRatio, ProtInfo, SILAC_data, tcols


# -----------------------------------------
# Apparent T50 with positivity
# -----------------------------------------
def apparent_t50_single(t: np.ndarray, y: np.ndarray) -> float:
    m = np.isfinite(y) & np.isfinite(t)
    tt = t[m].astype(float)
    yy = y[m].astype(float)
    if tt.size < 2 or np.nanmin(yy) < 0 or np.nanmax(yy) > 1:
        return float("nan")
    # b = exp(theta) > 0
    def obj(theta):
        b = math.exp(theta[0])
        return float(np.linalg.norm(yy - np.exp(-b * tt), ord=2))
    res = minimize(obj, x0=np.array([0.0]), method="Nelder-Mead",
                   options={"xatol":1e-12,"fatol":1e-12,"maxiter":5000,"maxfev":5000})
    b = math.exp(float(res.x[0]))
    if not np.isfinite(b) or b <= 0:
        return float("nan")
    return float(LN2 / b)


# -----------------------------------------
# Physics filter
# -----------------------------------------
def apply_physics_filter(Lys: np.ndarray, Prot: np.ndarray, mode: str, min_obs: int):
    """
    Lys: (T,), Prot: (T, N)
    Returns Prot_filtered and keep_mask (N,)
    """
    T, N = Prot.shape
    keep = np.ones(N, dtype=bool)
    eps = 1e-6

    if mode.lower() == "clip":
        for j in range(N):
            y = Prot[:, j].copy()
            bad = np.isfinite(y) & np.isfinite(Lys) & (y <= Lys + 1e-12)
            y[bad] = np.minimum(np.maximum(Lys[bad] + eps, 0), 1.0)
            if np.isfinite(y).sum() < min_obs:
                keep[j] = False
            Prot[:, j] = y

    elif mode.lower() == "drop":
        for j in range(N):
            y = Prot[:, j]
            if np.any(np.isfinite(y) & np.isfinite(Lys) & (y <= Lys + 1e-12)) or (np.isfinite(y).sum() < min_obs):
                keep[j] = False
        Prot = Prot[:, keep]

    else:
        # no physics, just enforce min_obs
        for j in range(N):
            if np.isfinite(Prot[:, j]).sum() < min_obs:
                keep[j] = False

    return Prot[:, keep], keep


# -----------------------------------------
# Binning + multistart seeds
# -----------------------------------------
@dataclass
class BinBuild:
    data: 'DataStruct'
    free_lys_weight_fn: callable

def build_binning(P: dict) -> BinBuild:
    t, LysRatio, ProtInfo, SILAC, tp_cols = _read_input_like_matlab(P["input_file"])

    # physics filter
    phys_enabled = (int(P.get("disable_physics_filter", 0)) == 0)
    min_obs = int(P.get("min_observations_per_protein", 3))
    if phys_enabled:
        SILAC, keep_mask = apply_physics_filter(LysRatio, SILAC.copy(), str(P.get("physics_filter_mode","clip")), min_obs)
        ProtInfo = ProtInfo.loc[keep_mask].reset_index(drop=True)
    else:
        SILAC, keep_mask = apply_physics_filter(LysRatio, SILAC.copy(), "none", min_obs)
        ProtInfo = ProtInfo.loc[keep_mask].reset_index(drop=True)

    # optional sorting
    sort_by_t50 = int(P.get("sort_by_apparent_t50", 1)) == 1
    if sort_by_t50:
        appr = np.array([apparent_t50_single(t, SILAC[:, j]) for j in range(SILAC.shape[1])], float)
        order = np.argsort(appr, kind="mergesort")
        SILAC = SILAC[:, order]
        ProtInfo = ProtInfo.iloc[order, :].reset_index(drop=True)

    # deterministic bins
    N = SILAC.shape[1]
    bin_size = int(P.get("bin_size", 30))
    if N < bin_size:
        bin_size = N
    protBins = [0]
    while protBins[-1] + bin_size < N:
        protBins.append(protBins[-1] + bin_size)
    if protBins[-1] != N:
        protBins.append(N)

    # multistart seed matrix (deterministic)
    rng = np.random.RandomState(int(P.get("random_seed", 12345)))
    Ansatz_ = rng.rand(int(P.get("n_starts", 20)), N + 2)

    y0_mode = str(P.get("initial_condition_mode","ones")).lower()
    wmode   = str(P.get("free_lys_weight_mode","auto")).lower()

    def free_lys_weight_for_bin(m: int) -> float:
        # m = proteins in current bin
        if wmode == "auto":          # strong anchor to free-Lys (previous default)
            return float(m)
        if wmode == "matlab" or wmode == "auto_div5":
            return max(1.0, m/5.0)   # wF ≈ M/5
        # numeric override
        try:
            return float(wmode)
        except:
            return float(m)

    data = DataStruct(
        t=t.astype(float),
        LysRatio=LysRatio.astype(float),
        ProtInfo=ProtInfo,
        SILAC_data=SILAC.astype(float),
        timepoints=tp_cols,
        protBins=protBins,
        Ansatz_=Ansatz_,
        physics_filter_mode=str(P.get("physics_filter_mode","clip")),
        free_lys_weight_mode=wmode,
        y0_mode=y0_mode
    )
    return BinBuild(data=data, free_lys_weight_fn=free_lys_weight_for_bin)


# -----------------------------------------
# Exact LTI propagation via matrix exponential
# -----------------------------------------
def simulate_lti_exact(t: np.ndarray, A: np.ndarray, b: np.ndarray, y0: np.ndarray) -> np.ndarray:
    n = A.shape[0]
    Y = np.empty((t.size, n), float)
    Y[0] = y0
    M = np.zeros((n+1, n+1), float)
    M[:n, :n] = A
    M[:n, -1] = b
    for k in range(1, t.size):
        dt = float(t[k] - t[k-1])
        E = expm(M * dt)
        F = E[:n, :n]
        g = E[:n, -1]
        Y[k] = F @ Y[k-1] + g
    return Y


# -----------------------------------------
# Model matrices (γ_A, γ_i, α)
# -----------------------------------------
def build_A_b_from_gamma(p_gamma: np.ndarray, theta_food_impurity: float) -> Tuple[np.ndarray, np.ndarray]:
    """
    p_gamma = [γ_A, γ_1..γ_m, alpha]
    """
    gamma_A = float(p_gamma[0])
    m = p_gamma.size - 2
    gammas = np.asarray(p_gamma[1:1+m], float)
    alpha = float(p_gamma[-1])

    A = np.zeros((1+m, 1+m), float)
    b = np.zeros(1+m, float)

    A[0, 0] = -(gamma_A + alpha*np.sum(gammas))
    if m:
        A[0, 1:] = alpha * gammas
        A[1:, 0] = gammas
        A[1:, 1:] = -np.diag(gammas)
    b[0] = gamma_A * theta_food_impurity
    return A, b


def _initial_y0(mode: str, Y_obs_first: np.ndarray) -> np.ndarray:
    if mode == "data":
        y0 = np.nan_to_num(Y_obs_first)
        y0[~np.isfinite(y0)] = 1.0
        return y0
    return np.ones_like(Y_obs_first, dtype=float)  # LIGHT = 1 at t0


def weighted_residuals_gamma(p: np.ndarray, t: np.ndarray, Y_obs: np.ndarray,
                             theta_food_imp: float, ub: np.ndarray, wF: float, y0_mode: str) -> np.ndarray:
    p = np.clip(p, 0, ub)
    A, b = build_A_b_from_gamma(p, theta_food_imp)
    y0 = _initial_y0(y0_mode, Y_obs[0, :])
    Y_sim = simulate_lti_exact(t, A, b, y0)
    R = (Y_sim - Y_obs)
    R[:, 0] *= wF   # weight free-Lys residuals
    return R[np.isfinite(R)]


# -----------------------------------------
# Fit corrected half-lives (Setting-2)
# -----------------------------------------
def fit_corrected(binbuild: BinBuild, P: dict) -> pd.DataFrame:
    data = binbuild.data
    wF_of_bin = binbuild.free_lys_weight_fn

    theta_food = float(P["purity_of_SILAC_food"]) / 100.0
    impurity   = 1.0 - theta_food

    N = data.SILAC_data.shape[1]

    # bounds (convert HL bounds to γ bounds)
    HL_p_min = float(P["hl_tmin_days"]); HL_p_max = float(P["hl_tmax_days"])
    HL_A_min = float(P["hl_A_min_days"]); HL_A_max = float(P["hl_A_max_days"])
    alphaUB  = float(P["alphaUB"])

    gA_min = LN2 / HL_A_max
    gA_max = LN2 / HL_A_min
    gP_min = LN2 / HL_p_max
    gP_max = LN2 / HL_p_min

    ub_all = np.r_[np.full(1, gA_max), np.full(N, gP_max), np.full(1, alphaUB)]

    rows = []
    for bi in range(len(data.protBins)-1):
        lo = data.protBins[bi]; hi = data.protBins[bi+1]
        m = hi - lo
        if m <= 0:
            continue
        print(f"  Optimizing proteins {lo+1} to {hi} (of {data.protBins[-1]})")

        Y_obs = np.c_[data.LysRatio, data.SILAC_data[:, lo:hi]]  # (T, 1+m)
        ub = ub_all[:1+m+1]
        wF = wF_of_bin(m)

        # deterministic multistart in γ-space within bounds
        rng = np.random.RandomState(int(P.get("random_seed", 12345)))
        gA0 = rng.uniform(gA_min, gA_max, size=(int(P["n_starts"]), 1))
        gP0 = rng.uniform(gP_min, gP_max, size=(int(P["n_starts"]), m))
        a0  = rng.uniform(0.0, alphaUB, size=(int(P["n_starts"]), 1))
        starts = np.hstack([gA0, gP0, a0])

        # global: TRF in γ, pick by free-Lys SSE (like MATLAB tie-breaker)
        best = None; best_free_sse = np.inf
        for s in starts:
            res = least_squares(
                weighted_residuals_gamma, s,
                args=(data.t, Y_obs, impurity, ub, wF, data.y0_mode),
                method="trf", bounds=(np.zeros_like(ub), ub),
                xtol=1e-10, ftol=1e-10, gtol=1e-10,
                max_nfev=int(P.get("max_nfev_global", 800)), verbose=0
            )
            pg = np.clip(res.x, 0, ub)
            A, b = build_A_b_from_gamma(pg, impurity)
            y0 = _initial_y0(data.y0_mode, Y_obs[0, :])
            Ysim = simulate_lti_exact(data.t, A, b, y0)
            free_sse = float(np.nansum((Ysim[:, 0] - data.LysRatio)**2))
            if free_sse < best_free_sse:
                best_free_sse = free_sse
                best = pg

        # refine: LM in log-γ (alpha fixed)
        alpha = float(best[-1])
        q0 = np.log(np.maximum(best[:-1], 1e-12))

        def resid_q(q):
            pg = np.exp(q)            # positive γ
            p_full = np.r_[pg, alpha] # + α (fixed in refine)
            return weighted_residuals_gamma(p_full, data.t, Y_obs, impurity, ub, wF, data.y0_mode)

        try:
            res2 = least_squares(resid_q, q0, method="lm",
                                 xtol=1e-12, ftol=1e-12, gtol=1e-12,
                                 max_nfev=int(P.get("max_nfev_refine", 400)))
        except Exception:
            res2 = least_squares(resid_q, q0, method="trf",
                                 xtol=1e-12, ftol=1e-12, gtol=1e-12,
                                 max_nfev=int(P.get("max_nfev_refine", 400)))

        qhat = res2.x
        pg_hat = np.exp(qhat)              # γ_A, γ_i...
        p_full = np.r_[pg_hat, alpha]      # + α

        # CI in log-γ space -> delta to HL
        try:
            J = res2.jac
            mres = res2.fun.size; ppar = res2.x.size
            dof = max(mres - ppar, 1)
            s2 = float(res2.fun @ res2.fun) / dof
            cov_q = s2 * np.linalg.pinv(J.T @ J)
            se_q = np.sqrt(np.maximum(np.diag(cov_q), 0.0))
            tval = student_t.ppf(0.975, dof)
            ci_q_half = tval * se_q
        except Exception:
            ci_q_half = np.full_like(qhat, np.nan)

        # simulate once to get per-protein residuals
        A, b = build_A_b_from_gamma(p_full, impurity)
        y0 = _initial_y0(data.y0_mode, Y_obs[0, :])
        Ysim = simulate_lti_exact(data.t, A, b, y0)
        Prot_err = np.nansum((Y_obs[:, 1:] - Ysim[:, 1:])**2, axis=0)

        # to half-lives
        gA = float(pg_hat[0]); gP = np.asarray(pg_hat[1:], float)
        HL_A = float(LN2 / max(gA, 1e-12))
        HLs  = LN2 / np.maximum(gP, 1e-12)

        # delta: dHL/dq = -HL (since HL = ln2 * exp(-q))
        ci_q = ci_q_half if np.ndim(ci_q_half)==1 and ci_q_half.size==qhat.size else np.full_like(qhat, np.nan)
        ci_halfwidth = np.abs(np.r_[HL_A, HLs]) * ci_q  # half-widths (95% ~ t*SE)

        # pretty display & extremes flag like MATLAB sheet
        for k in range(m):
            hl = float(HLs[k])
            if hl < 0.5:
                disp = "<0.5"
                extreme = 1
            elif hl > 100:
                disp = ">100"
                extreme = 1
            else:
                disp = f"{hl:.3f}"
                extreme = 0

            rows.append([
                *list(data.SILAC_data[:, lo+k]),
                hl,
                disp,
                int(extreme),
                float(ci_halfwidth[1+k] if np.isfinite(ci_halfwidth[1+k]) else np.nan),
                float(Prot_err[k])
            ])

    # compose result DataFrame
    SILAC_T_names = list(data.timepoints)
    df_SILAC_T = pd.DataFrame(
        rows,
        columns=SILAC_T_names + [
            "HalfLife_in_days",
            "HalfLife_display",
            "Is_Extreme",
            "CI_halfwidth",
            "residual_error"
        ]
    )
    out = pd.concat([data.ProtInfo.reset_index(drop=True), df_SILAC_T], axis=1)
    return out


# -----------------------------------------
# Apparent export
# -----------------------------------------
def export_apparent(binbuild: BinBuild, input_file: str):
    data = binbuild.data
    N = data.SILAC_data.shape[1]
    t50 = np.array([apparent_t50_single(data.t, data.SILAC_data[:, j]) for j in range(N)], float)
    order = np.argsort(t50, kind="mergesort")

    cols = list(data.ProtInfo.columns)
    base = {}
    if len(cols) >= 1:
        base["Protein"] = data.ProtInfo.iloc[order, 0].astype(str).to_numpy()
    if len(cols) >= 2:
        base["Gene"]    = data.ProtInfo.iloc[order, 1].astype(str).to_numpy()

    out_app = pd.DataFrame(base)
    out_app["Apparent T50"] = t50[order]
    out_path = f"results_Apparent_T50_{os.path.basename(input_file)}"
    out_app.to_excel(out_path, index=False)
    print(f"  Completed exporting apparent half-lives → {out_path}")


# -----------------------------------------
# Main
# -----------------------------------------
def main(params_path: str):
    P = parse_params(params_path)
    input_file = str(P.get("input_file","")).strip()
    if not input_file:
        raise ValueError("input_file missing in params.")
    base = os.path.dirname(os.path.abspath(params_path))
    if not os.path.isabs(input_file):
        input_file = os.path.join(base, input_file)
    if not os.path.exists(input_file):
        raise FileNotFoundError(input_file)

    print("\n\n  Welcome to JUMPt Python (matched mode)...\n")
    binbuild = build_binning(P)
    print("  Completed reading input file.\n  Now fitting the protein data and calculating corrected half-lives...\n")

    corrected = fit_corrected(binbuild, P)

    out_corr = f"results_Corrected_T50_{os.path.basename(input_file)}"
    with pd.ExcelWriter(out_corr) as w:
        # parameter sheet for provenance
        pd.DataFrame([{
            "input_file": os.path.basename(input_file),
            "bin_size": int(P["bin_size"]),
            "purity_of_SILAC_food": float(P["purity_of_SILAC_food"]),
            "number_of_timepoints": int(P["number_of_timepoints"]),
            "apparent_T50_calculation": int(P["apparent_T50_calculation"]),
            "n_starts": int(P["n_starts"]),
            "random_seed": int(P["random_seed"]),
            "max_nfev_global": int(P["max_nfev_global"]),
            "max_nfev_refine": int(P["max_nfev_refine"]),
            "alphaUB": float(P["alphaUB"]),
            "hl_tmin_days": float(P["hl_tmin_days"]),
            "hl_tmax_days": float(P["hl_tmax_days"]),
            "hl_A_min_days": float(P["hl_A_min_days"]),
            "hl_A_max_days": float(P["hl_A_max_days"]),
            "disable_physics_filter": int(P["disable_physics_filter"]),
            "physics_filter_mode": str(P["physics_filter_mode"]),
            "min_observations_per_protein": int(P["min_observations_per_protein"]),
            "free_lys_weight_mode": str(P["free_lys_weight_mode"]),
            "initial_condition_mode": str(P["initial_condition_mode"]),
            "sort_by_apparent_t50": int(P["sort_by_apparent_t50"]),
        }]).to_excel(w, sheet_name="parameter_file", index=False)

        corrected.to_excel(w, sheet_name="results", index=False)

    print(f"\n  Completed exporting corrected half-lives → {out_corr}\n")
    if int(P.get("apparent_T50_calculation", 1)) == 1:
        export_apparent(binbuild, input_file)
    print("\n  *******  JUMPt program is complete *******\n")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--params", default="JUMPt_python.params")
    args = ap.parse_args()
    main(args.params)
