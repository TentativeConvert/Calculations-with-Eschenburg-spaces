## Notes

Version 2 implemented a completely new algorithm for finding the list of spaces (the algorithm described in the current main README.)

Version 3, differs from Version 2 mainly in the use of different containers:  the list of spaces is now realized as deque of deques, and deques are converted to vectors only before sorting & looking for pairs.  This change was institigated by problems with memory management and followed a series of trials/expriments sketched in the "timings" section below.

The crucial difference between std::vector and std::deque (for the purposes of this program) seems to be that vectors require a single contiguous block of memory, and thus may need to be moved around when they grow.  On the other hand, it seems that std::sort can operate much faster on a vector than it can operate on a deque (as the below timings show, despite conflicting reports on the web).

Version 4, i.e. the current version in the main folder, is the first version that is split up into several files.  Moreover, it has a more efficient (and more flexible) two-stage sorting algorithm.

## Timings

all times in seconds, in the format t₁/t₂, where  
t₁ = time for finding list of spaces  
t₂ = total time  

| version | R ≤ 50.000, xps | R ≤ 50.000, reh | R ≤ 100.000, reh | R ≤ 150.000, reh | R ≤ 200.000, reh |
| --- | :---: | :---: | :---: | :---: | :---: |
| v2 | 15/17 | 12/14 | 54/62 | 550/2.247 | [killed] | 
| v2_ad | 15/19 | | | | 
| v2_dd | 15/22 | 12/17 | 51/75 | 128/3.067 |  |
| v2_p1 | 16/18 | | | |      
| v2_sd  | 14/26 | | | | 
| v3 | 14/19 | 10/15 | 45/70 | 107/181 | 196/367 |
| v4 | 14/17 | 10/13 | 45/59 | 107/138 | 196/ ?  |

Short explanation of the version names:

| version | file-name | implemantion |
| --- | --- | --- |
| v2 | `esch-gen-2.c` |   list of spaces = array of vectors, one vector for each value of `\|r\|` |
| v2_ad | – | list of spaces = array of deques, one deque for each value of `\|r\|` |
| v2_dd | – | list of spaces = deque of deques, one deque for each value of `\|r\|` |
| v2_p1 | `esch-gen-2-p1.c` | list of spaces = array of vectors, one vector for each value of `p₁` |
| v2_sd | `esch-gen-2-single-deque.c` | list of spaces = single deque |
| v3 | `esch-gen-3.c` | list of spaces = deque of deques, conversion to vector before sorting |
| v4 | [current]  | list of all spaces as in v3, sorting in two stages |

