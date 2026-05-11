# Proposition: minimal quotient for the SHIFT/Lugano witness

## Statement

Let `L` be the eight-outcome SHIFT label and let `C = f(L)` be the four-valued Lugano output class obtained from
`f(0)=f(1)=0` and `f(+)=f(-)=1` componentwise. For uniformly distributed computational inputs `|abc>`, the ideal SHIFT
measurement satisfies `C = w_L(a,b,c)` with probability one and `C` is uniformly distributed over four values.

Consider any coarse-grained measurement outcome `Y` with at most `m` possible values, obtained from `L` by any stochastic
channel `P(Y|L)`. Let a decoder `d(Y)` output one of the four Lugano output classes. Then the optimal success probability is

    P_success <= min(m,4)/4.

Consequently:

- any coarse-graining with `m <= 3` outcomes has `P_success <= 3/4` and cannot violate the Lugano causal bound;
- the four-outcome quotient measurement `Y=C=f(L)` achieves `P_success=1`.

## Proof

For any decoder `d`, the success probability is

    P_success = sum_y P(Y=y, C=d(y)).

For each `y`, the optimal decoder chooses the class `c_y` maximizing `P(Y=y,C=c)`, hence

    P_success = sum_y max_c P(Y=y,C=c).

Let `c_y` be a maximizing class for `y`. Group the terms by class:

    P_success = sum_c sum_{y: c_y=c} P(Y=y,C=c)
              <= sum_c P(C=c) over classes that appear among the c_y.

At most `m` distinct classes can appear among `m` outcomes `y`. Since the four Lugano classes are uniformly distributed,

    P(C=c)=1/4.

Therefore

    P_success <= m/4    for m <= 4,

and the trivial upper bound 1 applies for `m>=4`. Hence `P_success <= min(m,4)/4`.

The bound is achieved by assigning one outcome to each of `m` Lugano classes for `m<=4`. In particular, for `m=4` the
quotient `Y=C=f(L)` achieves `P_success=1`.

## Numerical verification

The Python verification found:

- cross-error `P(C != f(L)) = 0.000e+00`;
- `P(C) = [0.25 0.25 0.25 0.25]`;
- stochastic optimum = `min(m,4)/4` for `m=1,...,8`;
- among all 4140 deterministic partitions of the 8 SHIFT labels, no partition with `<=3` blocks violates `3/4`;
- perfect success occurs if and only if the partition refines the f-output quotient:
  `all_perfect_refine_f = True`,
  `all_refine_f_perfect = True`.