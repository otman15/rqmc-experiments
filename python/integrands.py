#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 23 11:29:56 2026

@author: otman
"""

import numpy as np
from scipy.special import erf, ndtr, ndtri


def sumueu(x):
    x = np.asarray(x)
    return np.sum(x * np.exp(x), axis=-1) - x.shape[-1]


def smoothperb4(x, expw):
    x = np.asarray(x)
    s = x.shape[-1]
    omega = np.arange(1, s + 1, dtype=np.float64) ** (-expw)
    terms = 1.0 + (30.0 * x**2 * (1.0 - x)**2 - 1.0) * omega
    return np.prod(terms, axis=-1) - 1.0


def mc2(x):
    x = np.asarray(x)
    s = x.shape[-1]
    return np.prod((s - x) / (s - 0.5), axis=-1) - 1.0


def polynomial(x):
    x = np.asarray(x)
    s = x.shape[-1]
    a = np.arange(1, s + 1, dtype=np.float64) / s
    terms = 1.0 + (x - 0.5) * a
    return np.prod(terms, axis=-1) - 1.0


def oscillatory(x):
    """SSJ Oscillatory model, centered to match the Genz oscillatory form."""
    x = np.asarray(x)
    s = x.shape[-1]
    a = np.arange(1, s + 1, dtype=np.float64) / s
    mean = np.cos((s + 1.0) / 4.0) * np.prod(2.0 * np.sin(a / 2.0) / a)
    return np.cos(np.sum(a * x, axis=-1)) - mean


def gaussian(x):
    """SSJ Gaussian model with Genz peak sign convention and exact centering."""
    x = np.asarray(x)
    s = x.shape[-1]
    mean = ((np.sqrt(np.pi) / 2.0) * erf(1.0)) ** s
    return np.exp(-np.sum(x**2, axis=-1)) - mean


def _normal_scores(x):
    """Map uniforms to normal scores, matching SSJ NormalDist.inverseF01."""
    x = np.asarray(x, dtype=np.float64)
    low = np.nextafter(0.0, 1.0)
    high = np.nextafter(1.0, 0.0)
    return ndtri(np.clip(x, low, high))


def smoothgauss(x):
    """SSJ SmoothGauss benchmark: smooth centered function of normal scores."""
    z = _normal_scores(x)
    s = z.shape[-1]
    return ndtr(1.0 + np.sum(z, axis=-1) / np.sqrt(s)) - ndtr(
        1.0 / np.sqrt(2.0)
    )


def piecelingauss(x):
    """SSJ PieceLinGauss benchmark: centered positive part of normal sum."""
    z = _normal_scores(x)
    s = z.shape[-1]
    a = ndtr(-1.0) - np.exp(-0.5) / np.sqrt(2.0 * np.pi)
    return np.maximum(np.sum(z, axis=-1) / np.sqrt(s) - 1.0, 0.0) + a


def indsumnormal(x):
    """SSJ IndSumNormal benchmark: centered indicator of a normal sum."""
    z = _normal_scores(x)
    s = z.shape[-1]
    return (np.sum(z, axis=-1) / np.sqrt(s) >= 1.0).astype(float) - ndtr(-1.0)
