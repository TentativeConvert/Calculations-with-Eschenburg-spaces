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
int generate_lists(const INT_R& R);
int analyse_space(const INT_p& k1, const INT_p& k2, const INT_p& k3, 
		  const INT_p &l1, const INT_p& l2, const INT_p& l3);
//void RemoveSpaces(char* source);

int main(int argc, char* argv[])
{
  long long R;
  long long k1, k2, k3, l1, l2, l3;
  if (argc == 2) {
    if (sscanf(argv[1],"r=%lld",&R) == 1)
      return generate_lists((INT_R)R);
    else if (sscanf(argv[1],"[%lld,%lld,%lld,%lld,%lld,%lld]",
		    &k1,&k2,&k3,&l1,&l2,&l3) == 6)
      return analyse_space((INT_p)k1,(INT_p)k2,(INT_p)k3,(INT_p)l1,(INT_p)l2,(INT_p)l3);
  }
  return show_usage(argv[0]);
}


int analyse_space(const INT_p& k1, const INT_p& k2, const INT_p& k3, 
		  const INT_p &l1, const INT_p& l2, const INT_p& l3)
{
  Space E({k1,k2,k3},{l1,l2,l3});
  if (E.is_space())
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


int generate_lists(const INT_R& R){
  SpaceTupleList tuples_rs(R);
  tuples_rs.compute_KS_invariants();
  
  SpaceTupleList tuples_he(tuples_rs,Space::compareHomotopyType, "Looking for tuples of homotopy equivalent spaces ...");
  tuples_he.print("list1-he.txt", "Found the following numbers of tuples of homotopy equivalent spaces:");

  SpaceTupleList tuples_the(tuples_he,Space::compareTangentialHomotopyType, "Looking for tuples of tangentially homotopy equivalent spaces ...");
  tuples_the.print("list2-the.txt", "Found the following numbers of tuples of tangentially homotopy equivalent spaces:");
  
  SpaceTupleList tuples_homeo(tuples_the,Space::compareHomeomorphismType, "Looking for tuples of homeomorphic spaces ...");
  tuples_homeo.print("list3-homeo.txt", "Found the following numbers of tuples of homeomorphic spaces:");
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
