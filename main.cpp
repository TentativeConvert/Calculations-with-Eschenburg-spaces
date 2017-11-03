/* compile with: 
   make
*/
#include <functional>
#include <vector>
#include <limits>
#include <typeinfo>
#include "config.h"
#include "esch_space.h"     // class Space
#include "esch_tuples.h"    // class SpaceTuple
                            // = wrapper for deque< Space >
                            // class SpaceTupleList
                            // = wrapper for deque< Space_tuples >
int show_usage(const char* name);
int generate_lists(const INT_R& R);
int analyse_space(const INT_P& k1, const INT_P& k2, const INT_P& k3, 
		  const INT_P &l1, const INT_P& l2, const INT_P& l3);


int main(int argc, char* argv[])
{
  long long R;
  long long k1, k2, k3, l1, l2, l3;
  if (argc == 2) {
    if (sscanf(argv[1],"r=%lld",&R) == 1)
      return generate_lists((INT_R)R);
    else if (sscanf(argv[1],"[%lld,%lld,%lld,%lld,%lld,%lld]",
		    &k1,&k2,&k3,&l1,&l2,&l3) == 6)
      return analyse_space((INT_P)k1,(INT_P)k2,(INT_P)k3,(INT_P)l1,(INT_P)l2,(INT_P)l3);
  }
  return show_usage(argv[0]);
}


int analyse_space(const INT_P& k1, const INT_P& k2, const INT_P& k3, 
		  const INT_P &l1, const INT_P& l2, const INT_P& l3)
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
\n To analyse Eschenburg space described by certain parameters, enter:	\
\n									\
\n       %1$s [75,54,-51,46,32,0]					\
\n or    %1$s \"[75, 54, -51, 46, 32, 0]\"				\
\n									\
\n To generate various lists of positively curved Eschenburg spaces	\
\n with |r| < 5000, enter:						\
\n									\
\n        %1$s r=5000							\
\n									\
\n",name);
  printf(" -------------------------------------------------------------------- ");
  printf("\n Size of data types (see config.h & docs/limits.pdf): \n");
  printf("   INT_P:    %3ld bit\n",sizeof(INT_P)*8);
  printf("   INT_R:    %3ld bit\n",sizeof(INT_R)*8);
  printf("   INT_KS:   %3ld bit\n",sizeof(INT_KS)*8);
  printf("   FLOAT_KS: %3d bit significand\n",std::numeric_limits<FLOAT_KS>::digits);
  printf(" -------------------------------------------------------------------- \n");
  return 1;
}
