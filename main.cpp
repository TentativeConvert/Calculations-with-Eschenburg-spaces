/* compile with: 
   make
*/
#include <functional>
#include <vector>
#include <limits>
#include <typeinfo>
#include <string>
#include "config.h"
#include "esch_space.h"     // class Space
#include "esch_tuples.h"    // class SpaceTuple
                            // = wrapper for deque< Space >
                            // class SpaceTupleList
                            // = wrapper for deque< Space_tuples >
int show_usage(std::string myname);
int generate_lists(const INT_R& R, const size_t& max_tuples = DEFAULT_MAX_TUPLES_PER_TUPLE_SIZE_PER_FILE);
int analyse_space(const INT_P& k1, const INT_P& k2, const INT_P& k3, 
		  const INT_P &l1, const INT_P& l2, const INT_P& l3);


int main(int argc, char* argv[])
{
  long long R;
  long long k1, k2, k3, l1, l2, l3;
  long long max_tuples;
  if (argc >= 2) {
    if (sscanf(argv[1],"r=%lld",&R) == 1 && argc == 2)
      return generate_lists((INT_R)R);
    else if (sscanf(argv[1],"r=%lld",&R) == 1 && sscanf(argv[2],"print=%lld",&max_tuples) == 1 && argc == 3)
      return generate_lists((INT_R)R,(size_t)max_tuples);	
    else if (sscanf(argv[1],"[%lld,%lld,%lld,%lld,%lld,%lld]",
		    &k1,&k2,&k3,&l1,&l2,&l3) == 6 && argc == 2)
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

int generate_lists(const INT_R& R, const size_t& max_tuples){
  SpaceTupleList tuples_rs(R);
  tuples_rs.compute_KS_invariants();
  
  SpaceTupleList tuples_he(tuples_rs,Space::compareHomotopyType, "homotopy classes");
  tuples_he.print("list1-he.txt", max_tuples);
  
  SpaceTupleList tuples_the(tuples_he,Space::compareTangentialHomotopyType, "tangential homotopy classes");
  tuples_the.print("list2-the.txt", max_tuples);
  
  SpaceTupleList tuples_homeo(tuples_the,Space::compareHomeomorphismType, "homeomorphism classes");
  tuples_homeo.print("list3-homeo.txt", max_tuples);
  return 0;
}


int show_usage(std::string myname)
{
  int len = myname.length();
  if (len > 30)
    myname.replace(8,len-26,"...");
  const char* name = myname.c_str();
  printf("\
\n -------   PROGRAM FOR CALCULATIONS WITH ESCHENBURG SPACES    ------- \
\n To analyse Eschenburg space described by certain parameters, enter:	\
\n									\
\n       %s [75,54,-51,46,32,0]					\
\n or    %s \"[75, 54, -51, 46, 32, 0]\"				\
\n									\
\n To generate various lists of tuples of positively curved Eschenburg	\
\n spaces with |r| < 5000, enter:					\
\n									\
\n        %s r=5000							\
\n									\
\n To limit the size of the output files, only a certain number of	\
\n tuples is written to each file by default.  This number can be	\
\n changed with the 'print=XXXX' option, e.g.:				\
\n									\
\n        %s r=5000 print=10000						\
\n									\
\n									\
\n",name,name,name,name); // using %1$s instead of %s does not work on Windows!
  printf(" -------------------------------------------------------------------- ");
  printf("\n Size of data types (see config.h & docs/limits.pdf): \n");
  printf("   INT_P:    %3ld bit\n",sizeof(INT_P)*8);
  printf("   INT_R:    %3ld bit\n",sizeof(INT_R)*8);
  printf("   INT_KS:   %3ld bit\n",sizeof(INT_KS)*8);
  printf("   FLOAT_KS: %3d bit significand\n",std::numeric_limits<FLOAT_KS>::digits);
  printf(" -------------------------------------------------------------------- \n");
  return 1;
}
