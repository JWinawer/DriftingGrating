"""Ring-toy front end: the analytic energy of Models/ring_toy/ring_toy.py on a ring at 6 deg.

Locations are L points of polar angle on a ring of radius ECC = 6 deg, so the distance between
two locations is the arc length 6 * (difference in polar angle, radians). A spatial pool of
s.d. s_deg (visual degrees) is a circular Gaussian in that arc length; the toy's s_sp in degrees
of polar angle is s_deg = s_sp * pi / 180 * 6 (its best value, 15 deg of polar angle, is 1.57 deg).

Energy per channel is the squared von Mises tuning of ring_toy.tuning(): E = c^2 exp(2 kappa_t
(cos 2(phi - theta) - 1)), so the toy's numerator En = (c * tuning)^n with n = 2 is E^1 here and
its semi-saturation sigma^n is sigma here. With kappa_t = 2 the energy tuning has a half-width at
half-height of 17 deg, close to the 8-orientation pyramid's 18 deg (cos^14).
"""
import numpy as np

from backend import Geometry, WEDGES, circdiff, cosd

ECC = 6.0


class RingGeometry(Geometry):
    def __init__(self, L=72, K=36):
        self.alpha = np.arange(L) * 360 / L
        self.theta = np.arange(K) * 180 / K
        self.masks = np.stack([np.abs(circdiff(self.alpha, w, 360)) <= 22.5 for w in WEDGES])
        self._ws = {}

    def blur(self, X, s_deg):
        if s_deg == 0:
            return X
        if np.isinf(s_deg):
            return np.broadcast_to(self.spatial_mean(X), X.shape)
        if s_deg not in self._ws:
            arc = np.deg2rad(circdiff(self.alpha[:, None], self.alpha[None, :], 360)) * ECC
            w = np.exp(-0.5 * (arc / s_deg) ** 2)
            self._ws[s_deg] = w / w.sum(1, keepdims=True)
        return X @ self._ws[s_deg].T


def stimuli(geom, contrast=0.5, kappa_t=2.0):
    """Energies (1, K, L) and labels for the 8 stimuli, in ring_toy.stimuli() order."""
    Es, stims = [], []
    for exp in ('dg', 'da'):
        for s in (0, 45, 90, 135):
            phi = np.full(len(geom.alpha), float(s)) if exp == 'dg' else (geom.alpha + s) % 180
            E = contrast ** 2 * np.exp(2 * kappa_t * (cosd(2 * (phi[None, :] - geom.theta[:, None])) - 1))
            Es.append(E[None])
            bar = s * np.ones(8) if exp == 'dg' else (WEDGES + s) % 180
            stims.append(dict(exp=exp, name=f'{exp}{s}', bar=bar))
    return Es, stims
