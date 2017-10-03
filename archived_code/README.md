## Notes on old versions

v0:  contains computations of s₂ & s₂₂ and test for "condition C" not contained in later versions

v2:  completely new algorithm for finding the list of spaces

v3 = current version (see main folder)  
list of spaces no realized as deque of deques  
conversion to vector before sorting & looking for pairs


## Timings

all times in seconds, in the format t1/t2, where  
t1 = time for finding list of spaces  
t2 = total time  

| version | R ≤ 50.000, xps | R ≤ 50.000, reh | R ≤ 100.000, reh | R ≤ 150.000, reh | R ≤ 200.000, reh|
|---|
| v2 | 15/17 | 12/14 | 54/62 | 550/2.247 | [killed] | 
| v2_ad | 15/19 | | | | 
| v2_dd | 15/22 | 12/17 | 51/75 | 128/3.067 |  |
| v2_p1 | 16/18 | | | |      
| v2_sd  | 14/26 | | | | 
| v3 | 14/19 | 10/15 | 45/70 | 107/181 | 196/367 |

Explanation of the version names:

| version | file-name | implemantion |
|---|
| v2 | `esch-gen-2.c` |   list of spaces = array of vectors, one vector for each value of `|r|` |
| v2_ad| – | list of spaces = array of deques, one deque for each value of `|r|` |
| v2_dd| – list of spaces = deque of deques, one deque for each value of `|r|` |
| v2_p1| `esch-gen-2-p1.c` | list of spaces = array of vectors, one vector for each value of `p₁` |
| v2_sd| `esch-gen-2-single-deque.c` | list of spaces = single deque |
| v3| [current] | list of spaces = deque of deques, conversion to vector before sorting |


