#### Background

The "basic invariants" of an Eschenburg space `E` are:  

    |r|  - the order of H^4(E); an odd integer
     s   - an integer in [-|r|/2, |r|/2] describing the linking form
    p_1  - an integer in [0, |r|] describing the first Pontryagin class

We say that "the basic invariants of two Eschenburg spaces `E` and `E'` agree if

     |r| = |r'|
      s  =  +/- s' 
     p_1 = p_1'

(The invariant `s` changes sign under orientation-reversing homeomorphisms.)

    r, s, s_2       coincide  <---->  spaces homeomorphic
    r, s, s_22      coincide  <---->  spaces homotopy equivalent
    r, s, s_22, p_1 coincide  <---->  spaces tangentially homotopy equivalent
    

#### The task (high-level formulation)
Find all pairs of positively curved Eschenburg spaces `(E, E')` with `|r| <= R` whose basic invariants agree, for some given positive bound `R`.
     
#### The task (low-level formulation)
[CEZ06] prove that all positively curved Eschenburg spaces with `|r| <= R` can be parametrized by 
quadruples `(k1, k2, l1, l2)` with

    R >=  k1 >= k2 > l1 >= l2 >= 0                  (1)

-- see [CEZ06, Lemma 1.4] and [CEZ06, proof of Prop. 1.7].  Moreover, these quadruples are required to satisfy a list of coprimacy conditions [CEZ06, (1.1)], see below.  (The conditions of [CEZ06, (1.2)] for positive curvature are automatically satisfied in this parametrization.)  

#### The algorithm
Instead of iterating over quadruples (`k1`,`k2`,`l1`,`l2`), the algorithm iterates over quadruples (`d`,`n`,`k1`,`k2`), where

    n := k2-l1
    d := k1-l2

Thus, below, `l1` and `l2` are to be read as short hands for `k2-n` and `k1-d`, respectively.
Note that `d >= n > 0`.  In terms of `n` and `d`, the formula for `|r|` can be written as

    |r| = d^2 + (k1+k2)*n + l2*(d-n).                (2a)

In particular:  

    |r| >= d^2 + d + 1.                              (2b)
    |r| >= d^2 + (1+k2)*n                            (2c)
        
The additional conditions [CEZ06, (1.1)] that the quadruples (`d`,`n`, `k1`, `k2`) need to satisfy are:

    (n,  d)    coprime
    (k2, d)    coprime
    (k1, n)    coprime
    (k2,    k1-l1)  (= (k2, k1-n))   coprime          (3)
    (k1,    k2-l2)  (= (k1, k2+d))   coprime
    (k1-l1, k2-l2)  (= (k1+k2-n,k2-k1+d))   coprime
    
Of course, there is a straight-forward way of finding all such quadruples:  simply iterate over all possible values of `k1`, `k2 `, `l1` and `l2` between `0`and `R`and check the conditions in each case.  The problem with this approach is that it is very inefficient (i.e. very slow).  The strategy we employ instead can be summarized as follows:

##### Step 1:  Find all coprime pairs `(n,d)` with `max_d >= d >= n > 0`, where `max_d := sqrt(R-3/4) - 1/2`
We employ the [standard algorithm for generating Farey sequences](https://en.wikipedia.org/wiki/Farey_sequence#Next_term) to find these.  (The pairs are interpreted as reduced fractions `0 < n/d <= 1`.  This explains our choice of letters: `n` for numerator and `d` for denominator.)  The necessity of the inequalities `d >= n > 0` follows from (1).  The necessity of the inequality `R' > d` follows from (2b).

##### Step 2:  Find all `k2` such that `(k2,d)` coprime with `max_k2 >= k2 > 0`, where `max_k2 := (R-d^2)/n - 1`.
Here, the inequality `k2 > 0` follows from (1) and `max_k2 >= k2` follows from (2c). 
To find all these `k2`, proceed in two substeps:

(a) First, search only in the range `d >= k2 > 0`.  Note that we need to allow the boundary case `d = k2` because `d` may be `1`.
On the other hand, we can of course make use of the upper bound, so really we search in the range
`min(d, max_k2) >= k2 > 0`.
(b) All remaining `k2` with `(k2,d)` coprime will be of the form `k2 = k2' + i*d` for some positive `i`, where `k2'` is as in (a). 

##### Step 3:  Find all `k1 >= k2` such that `(k1, n)` coprime, `|r| < R` and the remaining conditions of (1) are satisfied.
Using (2a), the condition `|r| < R` can be rewritten as:

    (R-n*(k2+d))/d >= k1

The remaining conditions of (1) translate as follows:

    l1 >= l2  <---->   k2+d-n >= k1
    l2 >= 0   <---->       k1 >= d

So the range we need to search in is `max_k1 > k1 >= min_k1` with
 
    min_k1 := max(d, k2)
    max_k1 := min(k2+d-n, (R-n*(k2+d))/d )    
   
Again, to find these `k1`, we proceed in two substeps:

(a) First, search only in the range `min_k1 + n > k1 >= min_k1`. Again, we should also take into account our upper bound, so the actual search will be in the range  `min(min_k1 + n, max_k1) > k1 >= min_k1`. 
(b) All remaining `k1` with `(k1,d)` coprime will be of the form `k1 = k1' + i*d` for some positive `i`, where `k1'`is as in (a).

##### Step 4:  Check the remaining conditions (3).

To save memory, in the actual algorithm the steps are interlaced -- as soon as we've found a Farey pair `(n,d)`, we look for a possible value of `k2`, as soon as we've found that value, we look for `k1`, etc. until we run out of possibilities;  then we proceed to the next Farey pair.
