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

#include <iostream>
#include <string>
#include <vector>

//////////////////////////////////////////////////

int show_usage(const char* name);
int generate_lists(const long& R);
int analyse_space(const long& k1, const long& k2, const long& k3, 
		  const long &l1, const long& l2, const long& l3);
void RemoveSpaces(char* source);

int main(int argc, char* argv[])
{
  long R;
  long k1, k2, k3, l1, l2, l3;
  if (argc == 2) {
    if (sscanf(argv[1],"r=%ld",&R) == 1)
      return generate_lists(R);
    else if (sscanf(argv[1],"[%ld,%ld,%ld,%ld,%ld,%ld]",&k1,&k2,&k3,&l1,&l2,&l3) == 6)
      return analyse_space(k1,k2,k3,l1,l2,l3);
  }
  return show_usage(argv[0]);
}


int analyse_space(const long& k1, const long& k2, const long& k3, 
		  const long &l1, const long& l2, const long& l3)
{
  Space E({k1,k2,k3},{l1,l2,l3});
  if (E.is_positively_curved_space())
    {
      E.compute_KS_invariants();
      E.print();
      return 0;
    }
  else
    {
      E.print();
      return 1;
    }
}


int generate_lists(const long& R){
  SpaceTupleList tuples_rs(R);
  tuples_rs.print("list_rs.txt", "Found the following numbers of tuples of spaces whose invariants r and s agree:");
  tuples_rs.compute_KS_invariants();
  
  SpaceTupleList tuples_he(tuples_rs,Space::compareHomotopyType, "Looking for homotopy equivalent tuples ...");
  tuples_he.print("list_he.txt", "Found the following numbers of tuples of homotopy equivalent spaces:");
  
  return 0;
}

int show_usage(const char* name)
{
  printf("\
\n -------   PROGRAM FOR CALCULATIONS WITH ESCHENBURG SPACES    ------- \
\n To analyses Eschenburg space described by certain parameters, enter:	\
\n									\
\n       %1$s [50,49,48,16,3,-2]					\
\n or    %1$s \"[50, 49, 48, 16, 6, -2]\"				\
\n									\
\n To generates various lists of positively curved Eschenburg spaces	\
\n with |r| < 5000, enter:						\
\n									\
\n        %1$s r=5000							\
\n									\
\n",name);
  return 1;
}
