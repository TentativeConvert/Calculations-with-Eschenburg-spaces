
# Mathematical background

See [\[CEZ07\]](#references) ...

Invariants of an Eschenburg space `E` with parameters (k₁,k₂,k₃,l₁,l₂,l₃) are:



# Usage

Once [installed](#installation), the program can be run from the command line:  simply navigate to the directory in which the program is installed and call `./esch` (on linux) or `esch.exe` (on windows).  Starting the program like this, without any additional parameters, will display some short usage instructions.  

### Analyse single space
To analyse the Eschenburg space `E` described by parameters `(k₁,k₂,k₃,l₁,l₂,l₃)`, enter one of
								
      ./esch [k₁,k₂,k₃,l₁,l₂,l₃]					
      ./esch "[k₁, k₂, k₃, l₁, l₂, l₃]"				

(or `esch.exe [k₁,k₂,k₃,l₁,l₂,l₃]` on windows systems).  The invariants computed are:

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
    
With the default configuration, results should be reliable for parameters `kᵢ` and `lᵢ` of absolute value at most ...   ([see details below](#configuration)).


### Count 'isomorphism' classes in a range
To count the number of various 'isomorphism' classes of positively curved Eschenburg spaces in a certain range, say with |r| < 5000, enter:						
								
     ./esch r=5000						

(or `esch.exe r=5000` on windows).  Output will be written to the following files in the same directory:

     list1-he.txt	      (homotopy classes)							
     list2-the.txt        (tangential homotopy classes)
     list3-homeo.txt      (homeomorphism classes)

In addition to the counting statistics, these files will contain lists of tuples that specify the same isomorphism class.  To limit the maximum file size, the maximum number of tuples listed is controlled by the [configuration variable](#configuration) `MAX_TUPLES_PER_TUPLESIZE_PER_FILE`, which can only be set at compile time.  

Note that the files will be overwritten the next time the program is run.  If you want to keep the results, make a copy of these files.

With the default configuration, results should be reliable up to parameters of ...   ([see details below](#configuration)).


# Installation

Two precompiled binaries are available in the repository:

    bin_nix64/esch        -- for 64 bit linux systems
    bin_win64/esch.exe    -- for 64 bit windows systems
    
(Both have been compiled on Ubuntu 14.04, using `gcc` and `mingw`, respectively.)  With a bit of luck, one of these will run on your system.  In this case, simply download the respective file to a folder of your choice.

If the above binaries are not appropriate for your system, download/clone the complete repository and compile from scratch.  Call `make win` or `make nix` to compile using the supplied `Makefile`s.  You will likely need to adapt the `Makefile`s to suit your system.   Your may also need to download and install some boost libraries.


# Configuration

When compiling from scratch, some options ....


# Implementation
The code is structured as follows:

      esch_space.*
    implement 
      class Space()
             
`esch_space.*` implements ....    

      esch_lists.*
      esch_generate.*
     implements 
       class SpaceTuple()
       class SpaceTupleList()

These are simple wrapper classes with interesting constructors:

One constructor, implemented separately in `esch_generate.cpp`, generates .... .  See [docs/esch_generate.md](docs/esch_generate.md), for details how this algorithm is implemented.

The other constructor takes an existing list and ....


## References
\[CEZ07\] [T. Chinburg, C. Escher and W. Ziller: *Topological properties of Eschenburg spaces and 3-Sasakian manifolds.*](https://doi.org/10.1007/s00208-007-0102-6)  Math. Ann. **339** (2007), no. 3, pp. 3–20
