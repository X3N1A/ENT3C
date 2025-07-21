import numpy as np
import pandas as pd


def vN_entropy(M, SUB_M_SIZE_FIX, CHRSPLIT, PHI_MAX, phi, BIN_TABLE):
    S = []
    BIN_TABLE_NEW = []
    N = M.shape[0]

    if np.isnan(SUB_M_SIZE_FIX):
        SUB_M_SIZE = np.round(N / CHRSPLIT)
    else:
        SUB_M_SIZE = SUB_M_SIZE_FIX

    SUB_M_SIZE = int(SUB_M_SIZE)
    WN = 1 + np.floor((N - SUB_M_SIZE) / phi)
    while WN > PHI_MAX:
        phi = phi + 1
        WN = 1 + np.floor((N - SUB_M_SIZE) / phi)
    WN = int(WN)

    R1 = np.array(range(0, N, phi))
    R1 = R1[0:WN]
    R2 = R1 + SUB_M_SIZE - 1
    R2 = R2[0:WN]
    R = np.stack((R1, R2), axis=1).astype(int)

    M = M.astype(float)  # cannot fill with nan otherwise
    M[M == 0] = np.nan
    M = np.log(M)
    # with open("outputPY.txt", "w") as FN:
    for rr in range(0, WN):
        # Python slicing is EXCLUSIVE of the stop index

        m = M[R[rr, 0] : R[rr, 1] + 1, R[rr, 0] : R[rr, 1] + 1].copy()

        # A view is a new array object that looks like a separate array but
        # shares the same data in memory as the original array. If you modify m,
        # you’re also modifying M in-place.
        # np.savetxt("matrix_python.isolated.txt", m, fmt="%.15g", delimiter="\t")

        # see GIT: NormM=0 and ChrY
        # print(np.any(np.sum((np.isnan(m) | (m == 0)), axis=1) >= (m.shape[0] - 1)))

        if np.all(np.isnan(m) | (m == 0)):
            ENT = np.nan

        else:
            # if any row/column is nan/0, ignore this column (0 comes from log(1)!)
            mask = np.sum((np.isnan(m) | (m == 0)), axis=1) < (SUB_M_SIZE)
            m = m[mask, :]
            m = m[:, mask]

            m[np.isnan(m)] = np.nanmin(m.flatten())
            P = np.corrcoef(m, rowvar=False)
            # problems with corrcoef and const. columns.
            # m = np.array([[12, 3, 2], [2, 3, 11], [1, 3, 7]])
            # m = np.array([[12, 0.3, 2], [2, 0.3, 11], [1, 0.3, 7]])
            # print(np.std(m, axis=0, ddof=1))
            # NumPy computes population std (N in denom), matlab the sample std ((N-1) in denom)!
            # I need np.std(m, axis=0, ddof=1)
            SDs = np.std(m, axis=0, ddof=1) < np.finfo(float).eps
            P[SDs, :] = 0
            P[:, SDs] = 0
            P[np.isnan(P)] = 0  #
            np.fill_diagonal(P, 1.0)
            rho = P / P.shape[0]
            # Compute the eigenvalues of a complex Hermitian or real symmetric matrix.
            vals = np.linalg.eigvalsh(rho)
            vals = vals[vals > np.finfo(float).eps]
            vals = np.real(vals)
            ENT = -np.sum(vals * np.log(vals))

        S.append(ENT)
        # print(f"Step {rr}: S = {ENT:.7f}")  # , file=FN)

        BIN_TABLE_NEW.append(
            [
                BIN_TABLE.iloc[R[rr, 0], BIN_TABLE.columns.get_loc("binNr")],
                BIN_TABLE.iloc[R[rr, 1], BIN_TABLE.columns.get_loc("binNr")],
                BIN_TABLE.iloc[R[rr, 0], BIN_TABLE.columns.get_loc("start")],
                BIN_TABLE.iloc[R[rr, 1], BIN_TABLE.columns.get_loc("end")],
            ]
        )

    S = np.array(S, dtype=np.float64)

    BIN_TABLE_NEW = pd.DataFrame(
        BIN_TABLE_NEW, columns=["binNr1", "binNr2", "start", "end"]
    )
    return S, SUB_M_SIZE, WN, phi, BIN_TABLE_NEW
