[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1156300.svg)](https://doi.org/10.5281/zenodo.1156300)

# Calculations with Eschenburg spaces

This C++ program implements and combines a subset of the features of the following two pieces of code:

-  The (unpublished) C program described by T. Chinburg, C. Escher and W. Ziller in *Topological properties of Eschenburg spaces and 3-Sasakian manifolds* ([\[CEZ07\]](#references)), which generates lists of pairs of Eschenburg spaces whose "basic polynomial invariants" (`|r|`, `|s|` and `p₁`) agree. 
-  The Maple program for computing invariants of Eschenburg spaces, also mentioned in [\[CEZ07\]](#references) and available on [W. Ziller's homepage](https://www.math.upenn.edu/~wziller/research.html).

The program was used to find the pairs of positively curved, tangentially homotopy euqivalent but non-homeomorphic pairs of positively curved Eschenburg spaces in ([\[GZ17\]](#references).

## Installation

Two precompiled binaries are available with the lastest [release](https://github.com/TentativeConvert/Eschenburg-calculator/releases/latest): 

    esch        -- an executable file for 64 bit Linux systems
    esch.exe    -- an executable file for 64 bit Windows systems
    
Both have been compiled on Ubuntu 14.04, using `gcc` and `mingw`, respectively.  With a bit of luck, one of these will run on your system.  In this case, simply download the respective file to a folder of your choice.

If these binaries do not work for you, you will need to download/clone the complete repository and compile from scratch.  Call `make win` or `make nix` to compile using the supplied `Makefile`s.  You will likely need to adapt the `Makefile`s to suit your system.   Your may also need to download and install some boost libraries.
There appear to be precompiled versions of MinGW, the GNU C compiler for Windows, in which the boost libraries are already included:
[nuwen.net/mingw.html](https://nuwen.net/mingw.html)  (no warranty).


## Usage

Once [installed](#installation), the program can be run from the command line:  simply navigate to the directory in which the program is installed and call `./esch` (on Linux) or `esch.exe` (on Windows).  Starting the program like this, without any additional parameters, will display some short usage instructions.  

### Analyse single space
To analyse the Eschenburg space E defined by the parameters (k,l) = (k₁,k₂,k₃,l₁,l₂,l₃), enter one of
								
      ./esch [k₁,k₂,k₃,l₁,l₂,l₃]					
      ./esch "[k₁, k₂, k₃, l₁, l₂, l₃]"				

For details on the invariants computed, see the [Mathematical background](#mathematical-background) below.

With the default configuration, the output of the program should be reliable for parameters kᵢ and lᵢ of absolute values up to 1500 (see [docs/limits.pdf](docs/limits.pdf) and [Configuration](#configuration) below).

### Count 'isomorphism' classes in a range
To count the number of various 'isomorphism' classes of positively curved Eschenburg spaces in a certain range, say with |r| < 5000, enter:						
								
     ./esch r=5000						

(or `esch.exe r=5000` on Windows).  Output will be written to the following files in the same directory:

|                       |                             |
| --------------------- | --------------------------- | 
|  `list1-he.txt`       | homotopy classes            |
|  `list2-the.txt`      | tangential homotopy classes |
|  `list3-homeo.txt`    | homeomorphism classes       |
     
In addition to the counting statistics, these files will contain lists of tuples that specify the same isomorphism class.  
See the examples files in the folder [bin_nix64](bin_nix64).  Note that the files will be overwritten the next time the program is run.  If you want to keep the results, make a copy of these files.  To limit the maximum file size, the maximum number of tuples listed in each file can be controlled with the command-line option `print = XXX`, e.g.

     ./esch r=5000 print=10000						

The default value is controlled by the [configuration variable](#configuration) `DEFAULT_MAX_TUPLES_PER_TUPLE_SIZE_PER_FILE`, which can only be set at compile time (current value is 1010).

With the default configuration for the data types used, results should be reliable up to values of `|r|` ≤ 600.000 (see [docs/limits.pdf](docs/limits.pdf) and [Configuration](#configuration) below).


## Mathematical background

For background on Eschenburg spaces, see [\[GZ\]](#references) or [\[CEZ07]\](#references).  Here we just summarize some of the key points.  The spaces were first introduced and studied in [\[Esch82\]](#references).  Each Eschenburg space is a (bi)quotient of SU(3) by an action of S¹.  Following [\[CEZ07\]](#references), we specify this action of S¹ and the resulting Eschenburg space E = E(k,l) by a six-tuple of integer parameters (k,l) = (k₁,k₂,k₃,l₁,l₂,l₃), which we refer to as the **parameter vector** of E.  The parameters need to satisfy k₁+k₂+k₃ = l₁+l₂+l₃, and some further conditions that ensure that the S¹-action is free, see [\[CEZ07, (1.1)\]](#references).  

### Positive curvature

Eschenburg defines a torus-invariant metric on a compact Lie group to be a left-invariant metric which is bi-invariant with respect to some maximal torus [\[Esch84, §34\]](#references).  Let us define a **left Eschenburg metric** on an Eschenburg space E(k,l) to be a quotient metric of a metric on SU(3) that is torus-invariant and invariant with respect to the S¹-action specified by the parameter vector (k,l).  An Eschenburg space is said to be **positively curved** if it has positive sectional curvature with respect to *some* left Eschenburg metric.  This positive curvature condition translates into a simple condition on the parameters, see [\[Esch84, Satz 414\]](#references) or [\[CEZ07, (1.2)\]](#references). 

### Normalized parameter vectors

The diffeomorphism class of an Eschenburg space E = E(k,l) is invariant under the following transformations of its parameter vector (cf. [\[CEZ07, above Lemma 1.4\]](#references)):

- Permute the parameters (k₁,k₂,k₃), or permute the parameters (l₁,l₂,l₃).
- Simultaneously switch the signs of all parameters.
- Simultaneously add a fixed integer n to each parameter.

Moreover, the positive curvature condition is also unaffected by these changes.  It follows that any positively curved Eschenburg space is diffeomorphic to an Eschenburg space E = E(k,l) with 
   
    k₁ ≥ k₂ > l₁ ≥ l₂ ≥ l₃ = 0     (*)
 
 [\[CEZ07, Lemma 1.4\]](#references).  More precisely, whenever E(k',l') is positively curved, we can find a parameter vector (k,l) satisfying (\*) and left Eschenburg metrics g and g' such that (E(k',l'),g') and (E(k,l),g) are isometric and have positive sectional curvature.  We call parameter vectors satisfying (\*) **normalized**.  All parameter vectors generated by the program when in “counting mode” are normalized in this sense.

<!--Whenever E(k,l) and E(k',l') are equal in this sense, they are defined by equivalent actions: there exists an isometry SU(3) → SU(3), i.e. a diffeomorphism that respects the above left-invariant metric on SU(3), which is equivariant with respect to the action of S¹ defined by (k,l) on the source and the action of S¹ defined by (k',l') on the target.

<!-- In the implementation notes [docs/esch_generate.md](docs/esch_generate.md), we refer to such a choice of parameters for E as **standard parametrization**. -->

### Invariants 

The invariants computed by the program are:

|         | range                       | [\[CEZ07\]](#references)  | [\[Mil00\]](#references)    | definition/interpretation                    |
| ------- | --------------------------- | ------------------------- | ---- | ------------------------------------------------------------------- |
| `r`     | ∈ ℕ                        | r = \|r(k,l)\|            | L₂   | \|σ₂(k) - σ₂(l)\| = order of H⁴(E)	           |			   
| `s`	  | ∈ (-`r`/2, `r`/2) ∩ ℤ | s			    | L₃   | representative of σ₃(k) - σ₃(l) ∈ ℤ/`r` (determines linking form)  |   
| `Σ`	  | ∈ {-1, 0, 1}		| –			    | r	   | representative of σ₁(l) ∈ ℤ/3        | 
| `p₁`	  | ∈ {0, ..., `r`-1}       | p₁			    | –	   | representative of first Pontryagin class ∈	 H⁴(E) = ℤ/`r`  | 
| `s₂`	  | ∈ (-1/2, 1/2]		| s₂			    | –	   | representative of KS-invariant s₂ ∈ ℚ/ℤ	 | 
| `s₂₂`	  | ∈ (-1/2, 1/2]		| s₂₂			    | –	   | representative of KS-invariant s₂₂ ∈ ℚ/ℤ	 | 

Here, the left column denotes the notation for the various invariants used by the program; 
the third and fourth columns indicate the notation used in the literature. 
Evidently, we are trying to follow the notation of [\[CEZ07\]](#references)  as closely as possible. 
Note that the invariants `s₂` and `s₂₂` can only be computed when the parameters  `(k₁,k₂,k₃,l₁,l₂,l₃)`  satisfy a certain (weak) condition called "condition C" [\[CEZ07, §2\]](#references).  For a translation between the parametrization of Eschenburg spaces used in [\[Mil00\]](#references) and the parametrization used here, see [\[Sha02, end of §2\]](#references).

### Classification

The table below summarizes which invariants need to agree in order for two Eschenburg spaces satisfying condition C to have the same homotopy type/homeomorphism type etc. 
For example, the first line says that, according [\[Kru98\]](#references), two Eschenburg spaces are homotopy equivalent through an orientation-preserving homotopy equivalence if and only if their invariants `r`, `s`, and `s₂₂` agree.  

| invariants … agree                                                 |⇔|   spaces agree up to …                   | References                                                          |
| ------------------------------------------------------------------ | --- | ------------------------------------- | -------------------------------------------------------------------- | 
| `r`, `s`, `s₂₂`       |⇔| oriented homotopy equivalence            | [\[Kru98\]](#references) |
| `r`, `s`, `s₂₂`, `p₁` |⇔| oriented tangential homotopy equivalence | [\[GZ\]](#references)    |
| `r`, `s`, `s₂`, `p₁`  |⇔| oriented homeomorphism                   | [\[Kru05\]](#references) <br>[\[CEZ07, Thm&nbsp;2.3\]](#references)  | 
|   |   |   |    |
| `r`, `\|s\|`, `\|s₂₂\|`, `sign(s)·sign(s₂₂)`          |⇔| homotopy equivalence            | [\[Mil00\]](#references) <br>[\[Kru98\]](#references)                |
| `r`, `\|s\|`, `\|s₂₂\|`, `sign(s)·sign(s₂₂)`, `p₁`    |⇔| tangential homotopy equivalence | [\[GZ\]](#references) |
| `\|r\|`, `\|s\|`, `\|s₂\|`,  `sign(s)·sign(s₂)`, `p₁` |⇔| homeomorphism                   | [\[Kru05\]](#references) <br>[\[CEZ07, Thm&nbsp;2.3\]](#references)  | 

In the second half of the table, the `sign` of an invariant is defined as one of the values `+1, 0, -1` in the evident way.  The sign of `s₂` is `0` if and only if `s₂ = 0` or `s₂ = 1/2` (since 1/2 = -1/2 in ℚ/ℤ), and likewise for `s₂₂`.   
    

## Code base
The code is structured as follows:

      esch_space.h            implement       class Space
      esch_space.cpp 
 
All code for computing the above invariants is contained in this class.  The class also contains comparison functions that determine whether two spaces are homotopy equivalent/tangentially homotopy equivalent/homeomorphic. 
 
      esch_tuples.h           implement       class SpaceTuple
      esch_tuples.cpp                         class SpaceTupleList
      esch_generate.cpp

The class `SpaceTuple` is a simple wrapper around `std::deque< Space >`.  The class `SpaceTupleList` is a wrapper around `std::deque < SpaceTuple >` with two interesting constructors:

One constructor, implemented separately in `esch_generate.cpp`, first generates a list of all positively curved Eschenburg spaces with `|r|` bounded by a given integer.  (More precisely, it generates a list of parameter values (k₁,k₂,k₃,l₁,l₂,l₃) that specify an Eschenburg space with `|r|` bounded by this integer.)  It then uses Milgram's homotopy invariants to find tuples of homotopy equivalent spaces on this list.  See [docs/esch_generate.md](docs/esch_generate.md) for further details.

The other constructor takes an existing list of tuples and a "filter" as input.  Possible "filters" are "tangential homotopy class" or  "homeomorphism class".  The constructor looks for (sub-)tuples of spaces in the given list of tuples that fall into the same isomorphism class according to the "filter".

The interface of these classes is demonstrated in `esch.cpp`.

#### Configuration

When compiling from scratch, the data types used in the computations and a few other options can be set in `config.h`.  See [docs/limits.pdf](docs/limits.pdf).


## References

|            |      |
| ---------- | ---- | 
| \[CEZ07\]  | [T. Chinburg, C. Escher and W. Ziller: *Topological properties of Eschenburg spaces and 3-Sasakian manifolds.*](https://doi.org/10.1007/s00208-007-0102-6) Math. Ann. **339** (2007), no. 3, 3–20. |
| \[Esch82\] | [J.-H. Eschenburg: *New examples of manifolds with strictly positive curvature.*](https://doi.org/10.1007/BF01389224) Invent. Math. **66** (1982), no. 3, 469–480. |
| \[Esch84\] | [J.-H. Eschenburg: *Freie isometrische Aktionen auf kompakten Lie-Gruppen mit positiv gekrümmten Orbiträumen.*](https://opus.bibliothek.uni-augsburg.de/opus4/frontdoor/index/index/docId/4219) Schriftenreihe des Mathematischen Instituts der Universität Münster, 2. Serie, Band 32, 1984. |
| \[Kru98\]  | [B. Kruggel: *Kreck-Stolz invariants, normal invariants and the homotopy classification of generalised Wallach spaces.*](https://doi.org/10.1093/qmathj/49.4.469) Quart. J. Math. Oxford Ser. (2) **49** (1998), no. 4, 469–485. |
| \[Kru05\]  | [B. Kruggel: *Homeomorphism and diffeomorphism classification of Eschenburg spaces.*](https://doi.org/10.1093/qmath/hah031) Quart. J. Math. Oxford Ser. (2) **56** (2005), no. 4, 553–577. |
| \[Mil00\]  | [R. J. Milgram: *The classification of Aloff-Wallach manifolds and their generalizations.*](https://mathscinet.ams.org/mathscinet-getitem?mr=1747543) Surveys on surgery theory, Vol. 1, 379–407, Ann. of Math. Stud., 145, Princeton Univ. Press, Princeton, NJ, 2000. |
| \[Sha02\]  | [K. Shankar: *Strong inhomogeneity of Eschenburg spaces.*](https://doi.org/10.1307/mmj/1022636754) Michigan Math. J. **50** (2002), no. 1, 125–141. |
| \[GZ\] | D. González-Álvaro, M. Zibrowius, *Open manifolds with positively curved souls* (in preparation). |

