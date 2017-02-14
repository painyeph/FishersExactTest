#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Author: eph

from math import log, exp

try:

    from math import lgamma  # Python 2.7+

except ImportError:

    _cache = {0: 0., 1: 0.}
    def lgamma(z):
        if z in _cache: return _cache[z]
        x = 0.1659470187408462e-06 / (z + 7)
        x += 0.9934937113930748e-05 / (z + 6)
        x -= 0.1385710331296526 / (z + 5)
        x += 12.50734324009056 / (z + 4)
        x -= 176.6150291498386 / (z + 3)
        x += 771.3234287757674 / (z + 2)
        x -= 1259.139216722289 / (z + 1)
        x += 676.5203681218835 / z
        x += 0.9999999999995183
        _cache[z] = x = log(x) - 5.58106146679532777 - z + (z-0.5)*log(z+6.5)
        return x


def _maxn():
    l = 1; n = 2; h = float('inf')
    while l < n:
        if abs(lgamma(n+1) - lgamma(n) - log(n)) >= 1: h = n
        else: l = n
        n = (l + min(h, l * 3)) // 2
    return n

LN10 = log(10)
NINF = float('-inf')
MAXN = _maxn()


# ======================== Full Test ========================

def test1(a, b, c, d):
    result = mlnTest2(a, a+b, a+c, a+b+c+d)
    return exp(-result[0]), exp(-result[1]), exp(-result[2])

def test2(a, ab, ac, abcd):
    result = mlnTest2(a, ab, ac, abcd)
    return exp(-result[0]), exp(-result[1]), exp(-result[2])

def mlnTest1(a, b, c, d):
    return mlnTest2(a, a+b, a+c, a+b+c+d)

def mlnTest2(a, ab, ac, abcd):
    if 0 > a or a > ab or a > ac or ab > abcd or ac > abcd: raise ValueError('invalid table')
    if abcd > MAXN: raise ValueError('math domain error')
    a_min = max(0, ab+ac-abcd)
    a_max = min(ab, ac)
    if a_min == a_max: return 0., 0., 0.
    i = int(round(ab * ac / abcd))
    p0 = lgamma(ab+1) + lgamma(ac+1) + lgamma(abcd-ac+1) + lgamma(abcd-ab+1) - lgamma(abcd+1)
    pi = lgamma(i+1) + lgamma(ab-i+1) + lgamma(ac-i+1) + lgamma(abcd-ab-ac+i+1)
    pa = lgamma(a+1) + lgamma(ab-a+1) + lgamma(ac-a+1) + lgamma(abcd-ab-ac+a+1)
    pl = pi if i < a else pa
    pr = pi if i > a else pa
    pt = pa
    sl = exp(pl - pa)
    sr = exp(pr - pa)
    st = 1.
    for i in range(a_min, a_max+1):
        if i == a: continue
        pi = lgamma(i+1) + lgamma(ab-i+1) + lgamma(ac-i+1) + lgamma(abcd-ab-ac+i+1)
        if i < a: sl += exp(pl - pi)
        if i > a: sr += exp(pr - pi)
        if pi >= pt: st += exp(pt - pi)
    if pl < pr:
        return -log(1.-max(0., exp(p0-pa)*(sr-1.))), max(0., pa-p0-log(sr)), max(0., pa-p0-log(st))
    else:
        return max(0., pa-p0-log(sl)), -log(1.-max(0., exp(p0-pa)*(sl-1.))), max(0., pa-p0-log(st))

def mlog10Test1(a, b, c, d):
    result = mlnTest2(a, a+b, a+c, a+b+c+d)
    return result[0]/LN10, result[1]/LN10, result[2]/LN10

def mlog10Test2(a, ab, ac, abcd):
    result = mlnTest2(a, ab, ac, abcd)
    return result[0]/LN10, result[1]/LN10, result[2]/LN10


# ======================== Left Tail Only ========================

def test1l(a, b, c, d):
    return exp(-mlnTest2l(a, a+b, a+c, a+b+c+d))

def test2l(a, ab, ac, abcd):
    return exp(-mlnTest2l(a, ab, ac, abcd))

def mlnTest1l(a, b, c, d):
    return mlnTest2l(a, a+b, a+c, a+b+c+d)

def mlnTest2l(a, ab, ac, abcd):
    if 0 > a or a > ab or a > ac or ab > abcd or ac > abcd: raise ValueError('invalid table')
    if abcd > MAXN: raise ValueError('math domain error')
    a_min = max(0, ab+ac-abcd)
    a_max = min(ab, ac)
    if a_min == a_max: return 0.
    i = int(round(ab * ac / abcd))
    p0 = lgamma(ab+1) + lgamma(ac+1) + lgamma(abcd-ac+1) + lgamma(abcd-ab+1) - lgamma(abcd+1)
    pa = lgamma(a+1) + lgamma(ab-a+1) + lgamma(ac-a+1) + lgamma(abcd-ab-ac+a+1)
    if i < a:
        return -log(1. - max(0., exp(p0 - pa) * sum(
            exp(pa - lgamma(i+1) - lgamma(ab-i+1) - lgamma(ac-i+1) - lgamma(abcd-ab-ac+i+1))
            for i in range(a+1, a_max+1) )))
    else:
        return max(0., pa - p0 - log(1. + sum(
            exp(pa - lgamma(i+1) - lgamma(ab-i+1) - lgamma(ac-i+1) - lgamma(abcd-ab-ac+i+1))
            for i in range(a_min, a) )))

def mlog10Test1l(a, b, c, d):
    return mlnTest2l(a, a+b, a+c, a+b+c+d)/LN10

def mlog10Test2l(a, ab, ac, abcd):
    return mlnTest2l(a, ab, ac, abcd)/LN10


# ======================== Right Tail Only ========================

def test1r(a, b, c, d):
    return exp(-mlnTest2r(a, a+b, a+c, a+b+c+d))

def test2r(a, ab, ac, abcd):
    return exp(-mlnTest2r(a, ab, ac, abcd))

def mlnTest1r(a, b, c, d):
    return mlnTest2r(a, a+b, a+c, a+b+c+d)

def mlnTest2r(a, ab, ac, abcd):
    if 0 > a or a > ab or a > ac or ab > abcd or ac > abcd: raise ValueError('invalid table')
    if abcd > MAXN: raise ValueError('math domain error')
    a_min = max(0, ab+ac-abcd)
    a_max = min(ab, ac)
    if a_min == a_max: return 0.
    i = int(round(ab * ac / abcd))
    p0 = lgamma(ab+1) + lgamma(ac+1) + lgamma(abcd-ac+1) + lgamma(abcd-ab+1) - lgamma(abcd+1)
    pa = lgamma(a+1) + lgamma(ab-a+1) + lgamma(ac-a+1) + lgamma(abcd-ab-ac+a+1)
    if i > a:
        return -log(1. - max(0., exp(p0 - pa) * sum(
            exp(pa - lgamma(i+1) - lgamma(ab-i+1) - lgamma(ac-i+1) - lgamma(abcd-ab-ac+i+1))
            for i in range(a_min, a) )))
    else:
        return max(0., pa - p0 - log(1. + sum(
            exp(pa - lgamma(i+1) - lgamma(ab-i+1) - lgamma(ac-i+1) - lgamma(abcd-ab-ac+i+1))
            for i in range(a+1, a_max+1) )))

def mlog10Test1r(a, b, c, d):
    return mlnTest2r(a, a+b, a+c, a+b+c+d)/LN10

def mlog10Test2r(a, ab, ac, abcd):
    return mlnTest2r(a, ab, ac, abcd)/LN10


# ======================== Two Tails Only ========================

def test1t(a, b, c, d):
    return exp(-mlnTest2t(a, a+b, a+c, a+b+c+d))

def test2t(a, ab, ac, abcd):
    return exp(-mlnTest2t(a, ab, ac, abcd))

def mlnTest1t(a, b, c, d):
    return mlnTest2t(a, a+b, a+c, a+b+c+d)

def mlnTest2t(a, ab, ac, abcd):
    if 0 > a or a > ab or a > ac or ab > abcd or ac > abcd: raise ValueError('invalid table')
    if abcd > MAXN: raise ValueError('math domain error')
    a_min = max(0, ab+ac-abcd)
    a_max = min(ab, ac)
    if a_min == a_max: return 0.
    p0 = lgamma(ab+1) + lgamma(ac+1) + lgamma(abcd-ac+1) + lgamma(abcd-ab+1) - lgamma(abcd+1)
    pt = lgamma(a+1) + lgamma(ab-a+1) + lgamma(ac-a+1) + lgamma(abcd-ab-ac+a+1)
    st = 1.
    for i in range(a_min, a):
        pi = lgamma(i+1) + lgamma(ab-i+1) + lgamma(ac-i+1) + lgamma(abcd-ab-ac+i+1)
        if pi < pt: break
        st += exp(pt - pi)
    for i in range(a_max, a, -1):
        pi = lgamma(i+1) + lgamma(ab-i+1) + lgamma(ac-i+1) + lgamma(abcd-ab-ac+i+1)
        if pi < pt: break
        st += exp(pt - pi)
    return max(0., pt - p0 - log(st))

def mlog10Test1t(a, b, c, d):
    return mlnTest2t(a, a+b, a+c, a+b+c+d)/LN10

def mlog10Test2t(a, ab, ac, abcd):
    return mlnTest2t(a, ab, ac, abcd)/LN10
