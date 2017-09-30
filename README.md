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
Lemma 1.4 of [CEZ06] shows that all positively curved Eschenburg spaces can be parametrized by quadruples `(k₁, k₂, l₁, l₂)` with

       k₁ ≥ k₂ > l₁ ≥ l₂ ≥ 0                              (1)

Moreover, these quadruples are required to satisfy a list of coprimacy conditions 

      (1.1) of [CEZ06]                                    (2)

which we will spell out below. 
(The conditions of [CEZ06, (1.2)] for positive curvature are automatically satisfied in this parametrization.) 
The task is thus to find all quadruples satisfying (1), (2) and

      |r| ≤ R                                             (3)

for some given upper bound `R`.  By [CEZ06, proof of Prop. 1.7], such an upper bound implies `R ≥ k₁`. 
So there is a straight-forward way of finding all such quadruples: 
simply iterate over all possible values of `k₁`, `k₂ `, `l₁` and `l₂` between `0`and `R`and check the conditions in each case. 
The problem with this approach is that it is very inefficient (i.e. very slow). 
We therefore use a slightly refined strategy, outlined below. 

#### The algorithm
Instead of iterating over quadruples (`k₁`,`k₂`,`l₁`,`l₂`), the algorithm iterates over quadruples (`d`,`n`,`k₁`,`k₂`), where

    n := k₂-l₁
    d := k₁-l₂

In terms of these quadrauples, conditions (1) are above are equivalent to the following conditions: 

    d ≥ n > 0                                             (1'a)
    k₁ ≥ d                                                (1'b)
    k₁ ≥ k₂ ≥ k₁+n-d                                      (1'c)

(The condition `d ≥ n` is in fact superfluous, but it will be useful to know that it holds.) 
The additional conditions [CEZ06, (1.1)] that the quadruples (`d`,`n`, `k₁`, `k₂`) need to satisfy are: 

    (n,  d)    coprime                                    (2'a)
    (k₁, n)    coprime                                    (2'b)
    (k₂, d)    coprime                                    (2'c)
    (k₂,    k₁-l₁)  (= (k₂, k₁-n))   coprime             \
    (k₁,    k₂-l₂)  (= (k₁, k₂+d))   coprime              (2'd)
    (k₁-l₁, k₂-l₂)  (= (k₁+k₂-n,k₂-k₁+d))   coprime      /


In terms of `n` and `d`, the formula for `|r|` can be written as: 

    |r| = k₁*d + k₂*n + d*n                

Thus, condition (3) is equivalent to: 

    R ≥ k₁*d + k₂*n + d*n                                  (3')

As `k₁ ≥ d`, `k₂ ≥ n ≥ 1` and `n ≥ 1`, this implies in particular that:

    R ≥ d² + d + 1                                         (3'a)
    R ≥ d² + n² + k₁*n                                     (3'b)
    
##### Step (a):  Find all coprime pairs `(n,d)` with `D ≥ d ≥ n > 0`, where `D := sqrt(R-3/4) - 1/2`
We employ the [standard algorithm for generating Farey sequences](https://en.wikipedia.org/wiki/Farey_sequence#Next_term) to find these. 
(The pairs are interpreted as reduced fractions `0 < n/d ≤ 1`. 
This explains our choice of letters: `n` for numerator and `d` for denominator.) 
The value of the upper bound `D` follows from (3'a). 

##### Step (b):  Find all `k₁` such that `(k₁,n)` coprime with `K₁ ≥ k₁ ≥ d`, where `K₁ := (R-d²-n²)/n`. 
The value of `K₁` follows from (3'b).  To find all these `k₁`, proceed in two substeps: 

(b.1) First, search only in the range `d+n > k₁ ≥ d`. 
      We can of course make use of the upper bound, so really we search in the range `min(d+n-1, K₁) ≥ k₁ ≥ d`. 
    
(b.2) All remaining `k₁` with `(k₁,d)` coprime will be of the form `k₁ = k₁' + i*d` for some positive `i`,
      where `k₁'` is as in (b.1).


##### Step (c):  Find all  `k₂` such that `(k₂, d)` coprime such that `k₁ ≥ k₂ ≥ k₁+n-d` and such that (3') is satisfied.
Condition (3') is equivalent to `(R-d*(k₁+n))/n ≥ k₁`.  So the range we need to search in is `K₂ > k₂ ≥ k₁+n-d` with
 
    K₂ := max((R-d*(k₁+n))/n, k₁).
   
Again, to find these `k₂`, we proceed in two substeps:

(a) First, search only in the range `k₁+n > k₂ ≥ k₁+n-d`.
    Again, we should also take into account our upper bound, so the actual search will be in the range  `min(k₁+n-1, K₂) ≥ k₁ ≥ k₁+n-d`.
    
(b) All remaining `k₂` with `(k₂,n)` coprime will be of the form `k₂ = k₂' + i*n` for some positive `i`, where `k₂'`is as in (a).

##### Step (d):  Check the remaining conditions (2d).

To save memory, in the actual algorithm the steps are interlaced -- as soon as we've found a Farey pair `(n,d)`, we look for a possible value of `k₂`, as soon as we've found that value, we look for `k₁`, etc. until we run out of possibilities;  then we proceed to the next Farey pair.