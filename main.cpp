/* compile with: 
   make
*/
#include <algorithm>
using std::sort;
#include <chrono>
using namespace std::chrono;
//////////////////////////////////////////////////
#include "aux_math.h"  // includes <cmath>
using std::abs;
using std::max;        // from <algorithm>
#include <boost/rational.hpp>
using boost::rational;
//////////////////////////////////////////////////
// eschenburg classes:
#include "esch_space.h"     // class Space
#include "esch_families.h"  // class Space_family
                            // = wrapper for deque< Space >
                            // class Deque_of_Space_families
                            // = wrapper for deque< Space_families >
#include "esch_generate.h"
//////////////////////////////////////////////////
// Tests (implementation further below)
void test_rationals(void);
void test_space(void);

//////////////////////////////////////////////////
// Main routine:

main(){
  //test_rationals();
  //test_Space();
  //std::string s = "abcd";

  long R;
  printf("\nMaximum value of |r|: ");
  if(scanf("%ld",&R) != 1 || R <= 0) return -1; 

  system_clock::time_point t1 = system_clock::now();/*timer*/
  Deque_of_Space_families families_rs;  
  generate_rs_families(families_rs, R, t1);

  //////////////////////////////////////////////////
  // We now have a deque of families (families_rs) whose invariants r & s agree.
  // Now look for those families for which also p1 agrees!

  printf("Looking for pairs whose invariants (r,s,s22) agree ...\n");

  long c_pairs = 0;
  long c_triples = 0;
  Deque_of_Space_families families_rss22;  

  long c_rsp_spaces = 0;
  long feedback_step_size = max((long)(families_rs.size()/100),(long)1);
  for(std::size_t i = 0; i < families_rs.size(); ++i)
    {
      if(i % feedback_step_size == 0) // feedback for user
	{
	  printf(" %3d%%\r",(int)((float)i/feedback_step_size));
	  fflush(stdout);
	}
      struct mycomparator {
	//bool operator() (const Space& E1, const Space& E2) { return (E1.p1() < E2.p1());}
	bool operator() (Space E1, Space E2) { return (E1.s22() < E2.s22());}  // <--- need to pass by value, not by reference for s22()
      } comparator;

      Space_family& F = families_rs[i];
      sort(F.begin(),F.end(),comparator);

      for(std::size_t i1 = 0; i1 < F.size(); )// i1 is incremented indirectly via i2
      {
	std::size_t i2 = i1+1;
	while (i2 < F.size() 
	       && (
		   ((F[i1].s() == F[i2].s() && F[i1].s22() == F[i2].s22()) 
		    || (F[i1].s() == -F[i2].s() && F[i1].s22() == -F[i2].s22())
		    || (abs(F[i1].s()) == abs(F[i2].s()) && F[i1].s22() == rational<long long>(1,2) && F[i2].s22() == rational<long long>(1,2) )
		    )
		   )
	       )
	  ++i2;
	if(i2 > i1+1)
	  {
	    //all_spaces with indexes i1,...,i2 define spaces with the same invariants
            Space_family new_family;
	    for(std::size_t j = i1; j < i2; ++j)
	      {
		new_family.push_back(F[j]);
	      }
	    families_rss22.push_back(new_family);
	  }
	i1 = i2;
      }
    }
  printf(" 100%%");
  // print first list of families to file:
  // (do this AFTER s22-values have been computed)
  families_rs.sort_and_count_families(); //stable_sort(families_rs.begin(),families_rs.end());
  families_rs.print("list_rs_human.txt");

  // print second family to file:
  families_rss22.sort_and_count_families();
  families_rss22.print("list_rss22_human.txt");

  system_clock::time_point t2 = system_clock::now();/*timer*/
  float duration = (float)duration_cast<milliseconds>( t2 - t1 ).count()/1000;/*timer*/
  //printf("\nFound %ld pairs and %ld triples in this range (total time: %.2f s).\n\n", c_pairs, c_triples,duration);
  return 0;
}


////////////////////////////////////////////////////////////////////////////////////////////////////

void test_rationals(void)
{ 
  rational<long long>a(10,2);
  if( a == rational<long long>(5,3))
    printf("schlecht\n");
  if (a == rational<long long>(5))
    printf("gut 1\n");
  if (a == rational<long long>(5,1))
    printf("gut 2\n");
}

void test_Space(void)
{  
  class Space E1;
  E1.setParameters({0,0,0},{0,0,0});
} 

