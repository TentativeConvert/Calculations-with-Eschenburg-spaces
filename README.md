#### Background

The "basic invariants" of an Eschenburg space `E` are:  

    |r|  - the order of H^4(E); an odd integer
     s   - an integer in (-|r|/2, |r|/2) describing the linking form
     p₁  - an integer in [0, |r|) describing the first Pontryagin class

We say that "the basic invariants of two Eschenburg spaces `E` and `E'` agree if

     |r| = |r'|
      s  = ±s' 
      p₁ = p₁'

(The invariant `s` changes sign under orientation-reversing homeomorphisms.)

    r, s, s₂₂       coincide  <---->  spaces homeomorphic
    r, s, s₂₂       coincide  <---->  spaces homotopy equivalent
    r, s, s₂₂, p₁   coincide  <---->  spaces tangentially homotopy equivalent
    

#### The task (high-level formulation)
Find all pairs of positively curved Eschenburg spaces `(E, E')` with `|r| ≤ R` whose basic invariants agree, for some given positive bound `R`.
     
#### The task (low-level formulation)
[CEZ06] prove that all positively curved Eschenburg spaces with `|r| ≤ R` can be parametrized by 
quadruples `(k₁, k₂, l₁, l₂)` with

    R ≥  k₁ ≥ k₂ > l₁ ≥ l₂ ≥ 0                  (1)

-- see [CEZ06, Lemma 1.4] and [CEZ06, proof of Prop. 1.7].  Moreover, these quadruples are required to satisfy a list of coprimacy conditions [CEZ06, (1.1)], see below.  (The conditions of [CEZ06, (1.2)] for positive curvature are automatically satisfied in this parametrization.)  

#### The algorithm
Instead of iterating over quadruples (`k₁`,`k₂`,`l₁`,`l₂`), the algorithm iterates over quadruples (`d`,`n`,`k₁`,`k₂`), where

    n := k₂-l₁
    d := k₁-l₂

Thus, below, `l₁` and `l₂` are to be read as short hands for `k₂-n` and `k₁-d`, respectively.
Note that `d ≥ n > 0`.  In terms of `n` and `d`, the formula for `|r|` can be written as

    |r| = d² + (k₁+k₂)*n + l₂*(d-n).                (2a)

In particular:  

    |r| ≥ d² + d + 1.                              (2b)
    |r| ≥ d² + (1+k₂)*n                            (2c)
        
The additional conditions [CEZ06, (1.1)] that the quadruples (`d`,`n`, `k₁`, `k₂`) need to satisfy are:

    (n,  d)    coprime
    (k₂, d)    coprime
    (k₁, n)    coprime
    (k₂,    k₁-l₁)  (= (k₂, k₁-n))   coprime          (3)
    (k₁,    k₂-l₂)  (= (k₁, k₂+d))   coprime
    (k₁-l₁, k₂-l₂)  (= (k₁+k₂-n,k₂-k₁+d))   coprime
    
Of course, there is a straight-forward way of finding all such quadruples:  simply iterate over all possible values of `k₁`, `k₂ `, `l₁` and `l₂` between `0`and `R`and check the conditions in each case.  The problem with this approach is that it is very inefficient (i.e. very slow).  The strategy we employ instead can be summarized as follows:

##### Step 1:  Find all coprime pairs `(n,d)` with `max_d ≥ d ≥ n > 0`, where `max_d := sqrt(R-3/4) - 1/2`
We employ the [standard algorithm for generating Farey sequences](https://en.wikipedia.org/wiki/Farey_sequence#Next_term) to find these.  (The pairs are interpreted as reduced fractions `0 < n/d ≤ 1`.  This explains our choice of letters: `n` for numerator and `d` for denominator.)  The necessity of the inequalities `d ≥ n > 0` follows from (1).  The necessity of the inequality `R' > d` follows from (2b).

##### Step 2:  Find all `k₂` such that `(k₂,d)` coprime with `max_k₂ ≥ k₂ > 0`, where `max_k₂ := (R-d²)/n - 1`.
Here, the inequality `k₂ > 0` follows from (1) and `max_k₂ ≥ k₂` follows from (2c). 
To find all these `k₂`, proceed in two substeps:

(a) First, search only in the range `d ≥ k₂ > 0`.  
    Note that we need to allow the boundary case `d = k₂` because `d` may be `1`.
    On the other hand, we can of course make use of the upper bound, so really we search in the range `min(d, max_k₂) ≥ k₂ > 0`.
    
(b) All remaining `k₂` with `(k₂,d)` coprime will be of the form `k₂ = k₂' + i*d` for some positive `i`, 
    where `k₂'` is as in (a). 

##### Step 3:  Find all `k₁ ≥ k₂` such that `(k₁, n)` coprime, `|r| < R` and the remaining conditions of (1) are satisfied.
Using (2a), the condition `|r| < R` can be rewritten as:

    (R-n*(k₂+d))/d ≥ k₁

The remaining conditions of (1) translate as follows:

    l₁ ≥ l₂  <---->   k₂+d-n ≥ k₁
    l₂ ≥ 0   <---->       k₁ ≥ d

So the range we need to search in is `max_k₁ > k₁ ≥ min_k₁` with
 
    min_k₁ := max(d, k₂)
    max_k₁ := min(k₂+d-n, (R-n*(k₂+d))/d )    
   
Again, to find these `k₁`, we proceed in two substeps:

(a) First, search only in the range `min_k₁ + n > k₁ ≥ min_k₁`. 
    Again, we should also take into account our upper bound, so the actual search will be in the range  `min(min_k₁ + n, max_k₁) > k₁ ≥ min_k₁`. 
    
(b) All remaining `k₁` with `(k₁,n)` coprime will be of the form `k₁ = k₁' + i*n` for some positive `i`, where `k₁'`is as in (a).

##### Step 4:  Check the remaining conditions (3).

To save memory, in the actual algorithm the steps are interlaced -- as soon as we've found a Farey pair `(n,d)`, we look for a possible value of `k₂`, as soon as we've found that value, we look for `k₁`, etc. until we run out of possibilities;  then we proceed to the next Farey pair.
