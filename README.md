# Calculations with Eschenburg spaces

This C++ program implements and combines a subset of the features of the following two pieces of code:

-  The (unpublished) C program described by T. Chinburg, C. Escher and W. Ziller in *Topological properties of Eschenburg spaces and 3-Sasakian manifolds* ([\[CEZ07\]](#references)), which generates lists of pairs of Eschenburg spaces whose "basic polynomial invariants" (`|r|`, `|s|` and `p₁`) agree. 
-  The Maple program for computing invariants of Eschenburg spaces, also mentioned in [\[CEZ07\]](#references) and available on [W. Ziller's homepage](https://www.math.upenn.edu/~wziller/research.html).


## Installation

Two precompiled binaries are available in the folders [bin_nix64](bin_nix64) and [bin_win64](bin_win64):

    bin_nix64/esch        -- for 64 bit Linux systems
    bin_win64/esch.exe    -- for 64 bit Windows systems
    
(Both have been compiled on Ubuntu 14.04, using `gcc` and `mingw`, respectively.)  With a bit of luck, one of these will run on your system.  In this case, simply download the respective file to a folder of your choice.

If the above binaries are not appropriate for your system, download/clone the complete repository and compile from scratch.  Call `make win` or `make nix` to compile using the supplied `Makefile`s.  You will likely need to adapt the `Makefile`s to suit your system.   Your may also need to download and install some boost libraries.
There appear to be precompiled versions of MinGW, the GNU C compiler for Windows, in which the boost libraries are already included:
[nuwen.net/mingw.html](https://nuwen.net/mingw.html)  (no warranty).


## Usage

Once [installed](#installation), the program can be run from the command line:  simply navigate to the directory in which the program is installed and call `./esch` (on Linux) or `esch.exe` (on Windows).  Starting the program like this, without any additional parameters, will display some short usage instructions.  

#### Analyse single space
To analyse the Eschenburg space `E` described by parameters `(k₁,k₂,k₃,l₁,l₂,l₃)`, enter one of
								
      ./esch [k₁,k₂,k₃,l₁,l₂,l₃]					
      ./esch "[k₁, k₂, k₃, l₁, l₂, l₃]"				

(or `esch.exe [k₁,k₂,k₃,l₁,l₂,l₃]` on Windows systems).  The invariants computed are:

     r   := σ₂(k) - σ₂(l) ∈ ℤ
         (|r| is the order of H⁴(E))
     s   := σ₃(k) - σ₃(l) ∈ ℤ/|r|, normalized to lie in the range (-|r|/2, |r|/2)
         (s determines the linking form)
     p₁  ∈ ℤ/|r|, normalized to be an integer in [0, |r|) 
         (p₁ is the first Pontryagin class)
     s₂  ∈ ℚ/ℤ normalized to lie in the range (-1/2, 1/2]
     s₂₂ ∈ ℚ/ℤ normalized to lie in the range (-1/2, 1/2]

For more details on these invariants, consult  [\[CEZ07\]](#references).  By  [\[CEZ07, Thm 2.3\]](#references), they classify positively curved Eschenburg spaces up to homotopy equivalence and homeomorphism as follows:

    |r|, |s|, |s₂₂|, sign(s·s₂₂)      agree   ⇔  spaces are homotopy equivalent
    |r|, |s|, |s₂₂|, sign(s·s₂₂), p₁  agree   ⇔  spaces are tangentially homotopy equivalent
    |r|, |s|, |s₂|,  sign(s·s₂),  p₁  agree   ⇔  spaces are homeomorphic

Note, however, that formulas for computing the invariants `s₂` and `s₂₂` are only known when the parameters `(k₁,k₂,k₃,l₁,l₂,l₃)` satisfy a certain `Condition C` [\[CEZ07, §2\]](#references).  'Most' Eschenburg spaces satisfy this condition.  For Eschenburg spaces that do no satisfy this condition, no homeomorphism classification seems to be known.
    
With the default configuration, the output of the program should be reliable for parameters `kᵢ` and `lᵢ` of absolute values up to `1500` (see [docs/limits.pdf](docs/limits.pdf) and [Configuration](#configuration) below).


#### Count 'isomorphism' classes in a range
To count the number of various 'isomorphism' classes of positively curved Eschenburg spaces in a certain range, say with `|r| < 5000`, enter:						
								
     ./esch r=5000						

(or `esch.exe r=5000` on Windows).  Output will be written to the following files in the same directory:

     list1-he.txt         (homotopy classes)							
     list2-the.txt        (tangential homotopy classes)
     list3-homeo.txt      (homeomorphism classes)

In addition to the counting statistics, these files will contain lists of tuples that specify the same isomorphism class.  
See the examples files in the folder [bin_nix64](bin_nix64).  Note that the files will be overwritten the next time the program is run.  If you want to keep the results, make a copy of these files.  To limit the maximum file size, the maximum number of tuples listed in each file can be controlled with the command-line option `print = XXX`, e.g.

     ./esch r=5000 print=10000						

The default value is controlled by the [configuration variable](#configuration) `DEFAULT_MAX_TUPLES_PER_TUPLE_SIZE_PER_FILE`, which can only be set at compile time (current value is `1010`).

With the default configuration for the data types used, results should be reliable up to values of `|r| ≤ 600.000` (see [docs/limits.pdf](docs/limits.pdf) and [Configuration](#configuration) below).


## Code base
The code is structured as follows:

      esch_space.h            implement       class Space
      esch_space.cpp 
 
All code for computing the above invariants is contained in this class.  The class also contains comparison functions that determine whether two spaces are homotopy equivalent/tangentially homotopy equivalent/homeomorphic. 
 
      esch_tuples.h           implement       class SpaceTuple
      esch_tuples.cpp                         class SpaceTupleList
      esch_generate.cpp

The class `SpaceTuple` is a simple wrapper around `std::deque< Space >`.  The class `SpaceTupleList` is a wrapper around `std::deque < SpaceTuple >` with two interesting constructors:

One constructor, implemented separately in `esch_generate.cpp`, first generates a list of all positively curved Eschenburg spaces with `|r|` bounded by a given integer.  (More precisely, it generates a list of parameter values `(k₁,k₂,k₃,l₁,l₂,l₃)` that specify an Eschenburg space `|r|` bounded by this integer.)  It then looks for spaces on this list whose parameters `|r|` and `|s|` agree and saves a list of tuples of such spaces.  See [docs/esch_generate.md](docs/esch_generate.md) for further details.

The other constructor takes an existing list of tuples and a "filter" as input.  Possible "filters" are "homotopy class" or  "homeomorphism class", for example.  The constructor looks for (sub-)tuples of spaces in the given list of tuples that fall into the same isomorphism class according to the "filter".

The interface of these classes is demonstrated in `esch.cpp`.

#### Configuration

When compiling from scratch, the data types used in the computations and a few other options can be set in `config.h`.  See [docs/limits.pdf](docs/limits.pdf).


## References
\[CEZ07\] [T. Chinburg, C. Escher and W. Ziller: *Topological properties of Eschenburg spaces and 3-Sasakian manifolds.*](https://doi.org/10.1007/s00208-007-0102-6)  Math. Ann. **339** (2007), no. 3, pp. 3–20
