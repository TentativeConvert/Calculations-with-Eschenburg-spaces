
## Mathematical background

See [\[CEZ07\]](#references) ...

Invariants of an Eschenburg space `E` with parameters (k₁,k₂,k₃,l₁,l₂,l₃) are:

     r   := ... - ....;
         |r| is the order of H⁴(E)   
     s   := ... - ... modulo |r|, normalized to be an integer in (-|r|/2, |r|/2);
         s determines the linking form
     p₁  := first Pontryagin class, normalized to be an integer in [0, |r|) 
     s₂  :=  ...
     s₂₂ :=  ...

## Usage

Once [installed](#Installation), the program can be run from the command line:  simply navigate to the directory in which the program is installed and call `./esch` (on linux) or `esch.exe` (on windows).  Starting the program like this, without any additional parameters, will display some short usage instructions.  

### Analyse single space
To analyse an Eschenburg space described by certain parameters, enter one of
								
      ./esch [75,54,-51,46,32,0]					
      ./esch "\[75, 54, -51, 46, 32, 0\]"				

(or `esch.exe [75,54,-51,46,32,0]` on windows systems).  
With the default configuration, results should be reliable up to parameters of ...   ([see details below](#configuration)).

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

## Installation

On Linux ....
(see also configuration below)

On Windows ...

## Configuration

When compiling from scratch, some options ....

## Implementation
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
