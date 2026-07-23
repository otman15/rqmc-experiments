#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 23 11:29:56 2026

@author: otman
"""

import numpy as np


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