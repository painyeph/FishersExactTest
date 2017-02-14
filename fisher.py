# eph 091204 created
# eph 100713 fixed: "return 0.0, 0.0, 0.0" in single tests
# eph 170214 using math.lgamma in Python 2.7+

from math import log, exp

LN10 = log(10)

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
    a_min = max(0, ab+ac-abcd)
    a_max = min(ab, ac)
    if a_min == a_max: return 0.0, 0.0, 0.0
    ps = []
    lpmax = -1e1000
    rpmax = -1e1000
    for i in range(a_min, a_max+1):
        p = _lnProbability(i, ab-i, ac-i, abcd-ab-ac+i)
        ps.append(p)
        if i <= a and lpmax < p: lpmax = p
        if i >= a and rpmax < p: rpmax = p
        pass
    tpmax = min(lpmax, rpmax)
    left_tail = 0
    right_tail = 0
    two_tails = 0
    c = a-a_min
    for i in range(a_max-a_min+1):
        if i <= c: left_tail += exp(ps[i]-lpmax)
        if i >= c: right_tail += exp(ps[i]-rpmax)
        if ps[i] <= ps[c]: two_tails += exp(ps[i]-tpmax)
    left_tail = max(0.0, -lpmax-log(left_tail))
    right_tail = max(0.0, -rpmax-log(right_tail))
    two_tails = max(0.0, -tpmax-log(two_tails))
    return left_tail, right_tail, two_tails

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
    a_min = max(0, ab+ac-abcd)
    a_max = min(ab, ac)
    if a_min == a_max: return 0.0
    ps = []
    lpmax = -1e1000
    for i in range(a_min, a_max+1):
        p = _lnProbability(i, ab-i, ac-i, abcd-ab-ac+i)
        ps.append(p)
        if i <= a and lpmax < p: lpmax = p
        pass
    left_tail = 0
    c = a-a_min
    for i in range(a_max-a_min+1):
        if i <= c: left_tail += exp(ps[i]-lpmax)
    left_tail = max(0.0, -lpmax-log(left_tail))
    return left_tail

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
    a_min = max(0, ab+ac-abcd)
    a_max = min(ab, ac)
    if a_min == a_max: return 0.0
    ps = []
    rpmax = -1e1000
    for i in range(a_min, a_max+1):
        p = _lnProbability(i, ab-i, ac-i, abcd-ab-ac+i)
        ps.append(p)
        if i >= a and rpmax < p: rpmax = p
        pass
    right_tail = 0
    c = a-a_min
    for i in range(a_max-a_min+1):
        if i >= c: right_tail += exp(ps[i]-rpmax)
    right_tail = max(0.0, -rpmax-log(right_tail))
    return right_tail

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
    a_min = max(0, ab+ac-abcd)
    a_max = min(ab, ac)
    if a_min == a_max: return 0.0
    ps = []
    lpmax = -1e1000
    rpmax = -1e1000
    for i in range(a_min, a_max+1):
        p = _lnProbability(i, ab-i, ac-i, abcd-ab-ac+i)
        ps.append(p)
        if i <= a and lpmax < p: lpmax = p
        if i >= a and rpmax < p: rpmax = p
        pass
    tpmax = min(lpmax, rpmax)
    two_tails = 0
    c = a-a_min
    for i in range(a_max-a_min+1):
        if ps[i] <= ps[c]: two_tails += exp(ps[i]-tpmax)
    two_tails = max(0.0, -tpmax-log(two_tails))
    return two_tails

def mlog10Test1t(a, b, c, d):
    return mlnTest2t(a, a+b, a+c, a+b+c+d)/LN10

def mlog10Test2t(a, ab, ac, abcd):
    return mlnTest2t(a, ab, ac, abcd)/LN10

# ======================== Common ========================

try:

    from math import lgamma  # Python 2.7+

    def _lnProbability (a, b, c, d):
        return (lgamma(a+b+1)+lgamma(c+d+1)+lgamma(a+c+1)+lgamma(b+d+1)
                -lgamma(a+1)-lgamma(b+1)-lgamma(c+1)-lgamma(d+1)-lgamma(a+b+c+d+1))

except ImportError:

    _cache = {0: 0, 1: 0}

    def _lnGamma (z):
        if z in _cache: return _cache[z]
        x = 0
        x += 0.1659470187408462e-06 / (z + 8)
        x += 0.9934937113930748e-05 / (z + 7)
        x -= 0.1385710331296526 / (z + 6)
        x += 12.50734324009056 / (z + 5)
        x -= 176.6150291498386 / (z + 4)
        x += 771.3234287757674 / (z + 3)
        x -= 1259.139216722289 / (z + 2)
        x += 676.5203681218835 / (z + 1)
        x += 0.9999999999995183
        _cache[z] = x = log(x) - 6.58106146679532777 - z + (z + 0.5) * log(z + 7.5)

    def _lnProbability (a, b, c, d):
        return (_lnGamma(a+b)+_lnGamma(c+d)+_lnGamma(a+c)+_lnGamma(b+d)
                -_lnGamma(a)-_lnGamma(b)-_lnGamma(c)-_lnGamma(d)-_lnGamma(a+b+c+d))
        return x
