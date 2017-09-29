## Background

The "basic invariants" of an Eschenburg space `E` are:  

    |r|  - the order of H^4(E); an odd integer)
    s    - an integer in [-|r|/2, |r|/2] describing the linking form
    p_1  - an integer in [0, |r|] describing the first Pontryagin class

We say that "the basic invariants of two Eschenburg spaces `E` and `E'` agree if

     |r| = |r'|
       s =  +/- s' 
     p_1 = p_1'

(The invariant `s` changes sign under orientation-reversing homeomorphisms.)

    r, s, s_2       coincide <=>  spaces homeomorphic
    r, s, s_22      coincide <=>  spaces homotopy equivalent
    r, s, s_22, p_1 coincide <=>  spaces tangentially homotopy equivalent


## The task (high-level formulation)
Find all pairs of positively curved Eschenburg spaces `(E, E')` with `|r| < R` whose basic invariants agree, for some given positive bound `R`.
     
## The task (low-level formulation)
[CEZ06] prove that all positively curved Eschenburg spaces with `|r| < R` can be parametrized by 
quadruples `(k1, k2, l1, l2)` with

    R >=  k1 >= k2 > l1 >= l2 >= 0

-- see [CEZ06, Lemma 1.4] and [CEZ06, proof of Prop. 1.7].  Moreover, these quadruples are required to satisfy a list of coprimacy conditions [CEZ06, (1.1)], see below.  (The conditions [CEZ06, (1.2)] are automatically satisfied.)  

## The algorithm
Instead of iterating over quadruples (`k1`,`k2`,`l1`,`l2`), the algorithm iterates over quadruples (`d`,`n`,`k1`,`k2`), where

    n := k2-l1
    d := k1-l2

Thus, below, `l1` and `l2` are to be read as short hands for `k2-n` and `k1-d`, respectively.
Note that `1 <= n <= d`.  In terms of `n` and `d`, the formula for `|r|` can be written as

    |r| = d^2 + (k1+k2)*n + l2*(d-n).          (UB)

In particular,  `|r| >= d^2 + d + 1`.
The additional conditions [CEZ06, (1.1)] that the quadruples (`d`,`n`, `k1`, `k2`) need to satisfy are:

    (n,  d)    coprime
    (k2, d)    coprime
    (k1, n)    coprime
    (k2,    k1-l1)  (= (k2, k1-n))   coprime
    (k1,    k2-l2)  (= (k1, k2+d))   coprime
    (k1-l1, k2-l2)  (= (k1+k2-n,k2-k1+d))   coprime
    
Of course, there is a straight-forward way of finding all such quadruples:  simply iterate over all possible values of `k1`, `k2 `, `l1` and `l2` between `0`and `R`and check the conditions in each case.  The problem with this approach is that it is very inefficient (i.e. very slow).  The strategy we employ instead can be summarized as follows:

#### Step 1:  Find all coprime pairs `(n,d)` with `R >= d >= n => 1`.
We use the [standard algorithm for generating Farey sequences](https://en.wikipedia.org/wiki/Farey_sequence#Next_term) to find these.  (The pairs are interpreted as reduced fractions `0 < n/d <= 1`.  This explains our choice of letters: `n` for numerator and `d` for denominator.)




