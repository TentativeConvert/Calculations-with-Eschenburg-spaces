## Background

"basic invariants":  `r, s, p_1`

    r, s, s_22      coincide  =>  spaces homotopy equivalent
    r, s, s_22, p_1 coincide  =>  spaces tangentially homotopy equivalent

If values of `s_2` differ, then spaces are not homeomorphic.

## Standard parametrization
  
The "standard parametrization" of positively curved Eschenburg spaces with `|r| < R` of [CEZ06] is given by quadruples `(k1, k2, l1, l2)` with

    R >=  k1 >= k2 > l1 >= l2 >= 0

-- see [CEZ06, Lemma 1.4] and [CEZ06, proof of Prop. 1.7].
These spaces are indeed positively curved -- as

    min(l1,l2,l3) = l3 = 0
    max(l1,l2,l3) = l1
    k1, k2 > l1
    k3 = l1 + l2 - k1 - k2 < 0

the conditions of [CEZ06, (1.2)] that no `ki` is equal to either the minimum or the maximums of the `li`s are automatically satisfied.

## Algorithm
Instead of iterating over `k1`, `k2`, `l1`, `l2`, the algorithm iterates over `d`, `n`, `k1`, `k2`, where

    n := k2-l1
    d := k1-l2

Note that `1 <= n <= d`.  In terms of `n` and `d`, the formula for `|r|` can be written as

    |r| = d^2 + (k1+k2)*n + l2*(d-n).

In particular,  `|r| >= d^2`.

   
