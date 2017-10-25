/* compile with: 
   make
*/
#include <functional>
// eschenburg classes:
#include "esch_space.h"     // class Space
#include "esch_tuples.h"    // class SpaceTuple
                            // = wrapper for deque< Space >
                            // class SpaceTupleList
                            // = wrapper for deque< Space_tuples >

//////////////////////////////////////////////////
// Main routine:
main(){
  long R;
  printf("\nMaximum value of |r|: ");
  if(scanf("%ld",&R) != 1 || R <= 0) return -1; 

  SpaceTupleList tuples_rs(R);
  tuples_rs.print("list_rs.txt", "Found the following numbers of tuples of spaces whose invariants r and s agree:");
  tuples_rs.compute_KS_invariants();

  SpaceTupleList tuples_he(tuples_rs,Space::compareHomotopyType, "Looking for homotopy equivalent tuples ...");
  tuples_he.print("list_he.txt", "Found the following numbers of tuples of homotopy equivalent spaces:");

  return 0;
}



