## Aim
The aim of this little snippet of code is to find all pairs of Eschenburg spaces `E` with bounded `|r| := |H⁴(E)|` whose basic polynomial invariants agree.  This task is subdivided into two (sub)tasks:

**Task 1:** Find all Eschenburg spaces with `|r| ≤ R` for some given bound `R`.  
**Task 2:** Given a list of all such spaces, find all pairs on this list whose invariants `r` and `s` agree.

The implementation of Task 1 is described in detail below.  Task 2 is fairly straight-forward.


## Background

The *basic polynomial invariants* of an Eschenburg space `E` are:  

    |r|  - the order of H⁴(E); an odd integer
     s   - an integer in (-|r|/2, |r|/2) describing the linking form
     p₁  - an integer in [0, |r|) describing the first Pontryagin class

We say that the basic polynomial invariants of two Eschenburg spaces `E` and `E'` *agree* if

     |r| = |r'|
      s  = ±s' 
      p₁ = p₁'

(The invariant `s` changes sign under orientation-reversing homeomorphisms.)  The agreement of these basic invariants does not itself have any geometric interpretation.  However, if in addition a certain invariant s₂₂ agrees, then the spaces are tangentially homotopy equivalent; if a certain invariant s₂ agrees, they are homeomorphic [[CEZ07](https://doi.org/10.1007/s00208-007-0102-6), Thm 2.3].  (Again, signs have to be taken into account in the definition of "agree".) 

    r, s,   s₂₂    agree  ⇔  spaces homotopy equivalent
    r, s, p₁, s₂₂  agree  ⇔  spaces tangentially homotopy equivalent
    r, s, p₁, s₂   agree  ⇔  spaces homeomorphic
    

## Implementation of Task 1
The task is to find all positively curved Eschenburg spaces `E` with `|r| ≤ R`, for some given positive bound `R`.
  
Lemma 1.4 of [[CEZ07](https://doi.org/10.1007/s00208-007-0102-6)] shows that all positively curved Eschenburg spaces can be parametrized by quadruples `(k₁, k₂, l₁, l₂)` with

       k₁ ≥ k₂ > l₁ ≥ l₂ ≥ 0                              (1)

Moreover, these quadruples are required to satisfy a list of coprimacy conditions 

      (1.1) of [[CEZ07](https://doi.org/10.1007/s00208-007-0102-6)]                                    (2)

which we will spell out below. 
(The conditions of [[CEZ07](https://doi.org/10.1007/s00208-007-0102-6), (1.2)] for positive curvature are automatically satisfied in this parametrization.) 
The task is thus to find all quadruples satisfying (1), (2) and

      |r| ≤ R                                             (3)

for some given upper bound `R`.  By [[CEZ07](https://doi.org/10.1007/s00208-007-0102-6), proof of Prop. 1.7], such an upper bound implies `R ≥ k₁`. 
So there is a straight-forward way of finding all such quadruples: 
simply iterate over all possible values of `k₁`, `k₂ `, `l₁` and `l₂` between `0`and `R`and check the conditions in each case. 
The problem with this approach is that it is very inefficient (i.e. very slow). 
We therefore use a slightly refined strategy, outlined below. 

### The algorithm
Instead of iterating over quadruples (`k₁`,`k₂`,`l₁`,`l₂`), the algorithm iterates over quadruples (`d`,`n`,`k₁`,`k₂`), where

    n := k₂-l₁
    d := k₁-l₂

In terms of these quadrauples, conditions (1) are above are equivalent to the following conditions: 

    d ≥ n > 0                                             (1'a)
    k₁ ≥ d                                                (1'b)
    k₁ ≥ k₂ ≥ k₁+n-d                                      (1'c)

(The condition `d ≥ n` is in fact superfluous, but it will be useful to know that it holds.) 
The additional conditions [[CEZ07](https://doi.org/10.1007/s00208-007-0102-6), (1.1)] that the quadruples (`d`,`n`, `k₁`, `k₂`) need to satisfy are: 

    (n,  d)    coprime                                    (2'a)
    (k₁, n)    coprime                                    (2'b)
    (k₂, d)    coprime                                    (2'c)
    (k₂,    k₁-l₁)  (= (k₂, k₁+n))   coprime             \
    (k₁,    k₂-l₂)  (= (k₁, k₂+d))   coprime              (2'd)
    (k₁-l₁, k₂-l₂)  (= (k₁-k₂+n,k₂-k₁+d))   coprime      /


In terms of `n` and `d`, the formula for `|r|` can be written as: 

    |r| = k₁d + nd + nk₂

Thus, condition (3) is equivalent to: 

    R ≥ k₁d + nd + nk₂                                     (3')

As `k₁ ≥ d` and `k₂ ≥ n ≥ 1`, this implies in particular that:

    R ≥ d² + d + 1                                         (3'a)
    R ≥ k₁d + nd + n²                                      (3'b)
    
##### Step (a):  Find all coprime pairs `(n,d)` with `D ≥ d ≥ n > 0`, where `D := sqrt(R-3/4) - 1/2`
We employ the [standard algorithm for generating Farey sequences](https://en.wikipedia.org/wiki/Farey_sequence#Next_term) to find these. 
(The pairs are interpreted as reduced fractions `0 < n/d ≤ 1`. 
This explains our choice of letters: `n` for numerator and `d` for denominator.) 
The value of the upper bound `D` follows from (3'a). 

##### Step (b):  Find all `k₁` such that `(k₁,n)` coprime with `K₁ ≥ k₁ ≥ d`, where `K₁ := (R-n²)/d - n`. 
The value of `K₁` follows from (3'b).  To find all these `k₁`, proceed in two substeps: 

(b.1) First, search only in the range `d+n > k₁ ≥ d`. 
      We can of course make use of the upper bound, so really we search in the range `min(d+n-1, K₁) ≥ k₁ ≥ d`. 
    
(b.2) All remaining `k₁` with `(k₁,n)` coprime will be of the form `k₁ = k₁' + in` for some positive integer `i`,
      where `k₁'` is one of the values found in (b.1).


##### Step (c):  Find all  `k₂` such that `(k₂, d)` coprime such that `k₁ ≥ k₂ ≥ k₁+n-d` and such that (3') is satisfied.
Condition (3') is equivalent to `(R-k₁d)/n - d ≥ k₂`.  So the range we need to search in is `K₂ ≥ k₂ ≥ k₁+n-d` with
 
    K₂ := min((R-k₁d)/n - d, k₁).
   
Again, to find these `k₂`, we proceed in two substeps:

(c.1) First, search only in the range `k₁+n > k₂ ≥ k₁+n-d`.
      Again, we should also take into account our upper bound, so the actual search will be in the range
      `min(k₁+n-1, K₂) ≥ k₂ ≥ k₁+n-d`.

(c.2) All remaining `k₂` with `(k₂,d)` coprime will be of the form `k₂ = k₂' + id` for some positive integer `i`, 
      where `k₂'`is one of the values found in (c.1).

##### Step (d):  Check the remaining conditions (2'd).

To save memory, in the actual algorithm the steps are interlaced -- as soon as we've found a Farey pair `(n,d)`, we look for a possible value of `k₁`, as soon as we've found that value, we look for `k₂`, etc. until we run out of possibilities;  then we proceed to the next Farey pair.

## Limits
All integers are implemented using the data type `long`, which can store values up to ±2³¹ (so more than ±10⁹). The biggest value occuring is the upper bound `R` on `|r|`.  So in theory, the programme can find all Eschenburg spaces with `|r| ≤ 2³¹`.  Calculations for `|r| ≤ 100.000` should complete within a few minutes on standard machines.

## Reference
[[CEZ07](https://doi.org/10.1007/s00208-007-0102-6)] T. Chinburg, C. Escher and W. Ziller: *Topological properties of Eschenburg spaces and 3-Sasakian manifolds.* Math. Ann. **339** (2007), no. 3, pp. 3–20
