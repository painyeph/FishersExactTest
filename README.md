## Functions

              | p-value                                                | -log( p-value )                                              | -log10( p-value )
--------------|--------------------------------------------------------|--------------------------------------------------------------|--------------------------------------------------------------------
 left-tailed  | `test1l(a, b, c, d)` or `test2l(a, a+b, a+c, a+b+c+d)` | `mlnTest1l(a, b, c, d)` or `mlnTest2l(a, a+b, a+c, a+b+c+d)` | `mlog10Test1l(a, b, c, d)` or `mlog10Test2l(a, a+b, a+c, a+b+c+d)`
 right-tailed | `test1r(a, b, c, d)` or `test2r(a, a+b, a+c, a+b+c+d)` | `mlnTest1r(a, b, c, d)` or `mlnTest2r(a, a+b, a+c, a+b+c+d)` | `mlog10Test1r(a, b, c, d)` or `mlog10Test2r(a, a+b, a+c, a+b+c+d)`
 two-sided    | `test1t(a, b, c, d)` or `test2t(a, a+b, a+c, a+b+c+d)` | `mlnTest1t(a, b, c, d)` or `mlnTest2t(a, a+b, a+c, a+b+c+d)` | `mlog10Test1t(a, b, c, d)` or `mlog10Test2t(a, a+b, a+c, a+b+c+d)`
 triple       | `test1(a, b, c, d)` or `test2(a, a+b, a+c, a+b+c+d)`   | `mlnTest1(a, b, c, d)` or `mlnTest2(a, a+b, a+c, a+b+c+d)`   | `mlog10Test1(a, b, c, d)` or `mlog10Test2(a, a+b, a+c, a+b+c+d)`

## Comparison with scipy.stats.fisher_exact

### Speed

```bash
$ python -m timeit "from scipy import stats" "stats.fisher_exact([[8, 2], [1, 5]])[1]"
1000 loops, best of 3: 931 usec per loop

$ python -m timeit "import fisher" "fisher.test1t(8, 2, 1, 5)"
100000 loops, best of 3: 4.99 usec per loop
```

### Precision

```python
>>> from scipy import stats
>>> -log10(stats.fisher_exact([[100, 1], [10, 1000]])[1])
128.93472935802367
>>> -log10(stats.fisher_exact([[100, 1], [10, 10000]])[1])
226.6210481678513
>>> -log10(stats.fisher_exact([[100, 1], [10, 100000]])[1])
Traceback (most recent call last):
  File "<stdin>", line 1, in <module>
ValueError: math domain error
```

```python
>>> import fisher
>>> fisher.mlog10Test1t(100, 1, 10, 1000)
128.93472935802276
>>> fisher.mlog10Test1t(100, 1, 10, 10000)
226.62104816784367
>>> fisher.mlog10Test1t(100, 1, 10, 100000)
326.3812661050023
>>> fisher.mlog10Test1t(100, 1, 10, 1000000)
426.3571991266363
>>> fisher.mlog10Test1t(100, 1, 10, 10000000)
526.3547915418933
>>> fisher.mlog10Test1t(100, 1, 10, 100000000)
626.354550887692
>>> fisher.mlog10Test1t(100, 1, 10, 1000000000)
726.3545256659792
>>> fisher.mlog10Test1t(100, 1, 10, 10000000000)
826.3545146294285
>>> fisher.mlog10Test1t(100, 1, 10, 100000000000)
926.354583114538
>>> fisher.mlog10Test1t(100, 1, 10, 1000000000000)
1026.3549166719595
>>> fisher.mlog10Test1t(100, 1, 10, 10000000000000)
1126.3970256263215
>>> fisher.mlog10Test1t(100, 1, 10, 100000000000000)
1226.447616894783
>>> fisher.mlog10Test1t(100, 1, 10, 1000000000000000)
Traceback (most recent call last):
  File "<stdin>", line 1, in <module>
OverflowError: the grand total of contingency table is too large
```
