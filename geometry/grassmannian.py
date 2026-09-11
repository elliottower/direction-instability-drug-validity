"""Grassmannian geometry for subspace-level validity corrections.

Provides geodesic distance, subspace overlap, Frechet mean on Gr(k,d),
and the essential-response correction (projecting out a confounding
subspace before recomputing distances).
"""
import numpy as np
from scipy import linalg


def principal_angles(U1: np.ndarray, U2: np.ndarray) -> np.ndarray:
    """Principal angles between two subspaces given as orthonormal bases."""
    M = U1.T @ U2
    svals = linalg.svdvals(M)
    svals = np.clip(svals, -1.0, 1.0)
    return np.arccos(svals)


def geodesic_distance(U1: np.ndarray, U2: np.ndarray) -> float:
    """Grassmannian geodesic distance = L2 norm of principal angles."""
    return float(np.linalg.norm(principal_angles(U1, U2)))


def subspace_overlap(U: np.ndarray, V: np.ndarray) -> float:
    """Mean cosine of principal angles -- 1.0 = aligned, 0.0 = orthogonal."""
    angles = principal_angles(U, V)
    return float(np.mean(np.cos(angles)))


def frechet_mean_subspace(subspaces: list[np.ndarray], n_iter: int = 50) -> np.ndarray:
    """Frechet mean on Gr(k, d) via iterative projection.

    Starting from the first subspace, repeatedly projects all subspaces
    onto the tangent space at the current estimate, averages, and maps back.

    Args:
        subspaces: list of (d, k) orthonormal basis matrices.
        n_iter: maximum iterations.

    Returns:
        (d, k) orthonormal basis of the Frechet mean subspace.
    """
    mean = subspaces[0].copy()
    for _ in range(n_iter):
        P = mean @ mean.T
        tangent_sum = np.zeros_like(mean)
        for U in subspaces:
            proj = P @ U
            Q, R = np.linalg.qr(proj)
            tangent_sum += Q
        new_mean, _ = np.linalg.qr(tangent_sum)
        if np.linalg.norm(new_mean @ new_mean.T - P, 'fro') < 1e-8:
            break
        mean = new_mean
    return mean


def project_out_subspace(
    U: np.ndarray,
    confound: np.ndarray,
) -> np.ndarray:
    """Project confounding subspace out of U on Gr(k, d).

    Removes the component of U that lies in the confound subspace,
    then re-orthogonalizes.

    Args:
        U: (d, k) orthonormal basis of the target subspace.
        confound: (d, k_c) orthonormal basis of the confounding subspace.

    Returns:
        (d, k) orthonormal basis of the residual subspace. If U lies
        entirely within confound, returns a random orthogonal complement.
    """
    d, k = U.shape
    P_confound = confound @ confound.T
    residual = U - P_confound @ U
    Q, R = np.linalg.qr(residual)
    rank = np.sum(np.abs(np.diag(R)) > 1e-8)
    if rank < k:
        complement = np.eye(d) - P_confound - Q[:, :rank] @ Q[:, :rank].T
        eigvals, eigvecs = np.linalg.eigh(complement)
        extra = eigvecs[:, -k + rank:]
        Q = np.column_stack([Q[:, :rank], extra])
        Q, _ = np.linalg.qr(Q)
    return Q[:, :k]
