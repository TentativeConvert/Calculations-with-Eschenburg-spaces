/* compile with: 
   make
*/
#include <algorithm>
using std::sort;

#include <functional>

#include "aux_feedback.h"
#include "aux_math.h"       // includes <cmath> & rationals

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

  Feedback feedback;

  long R;
  printf("\nMaximum value of |r|: ");
  if(scanf("%ld",&R) != 1 || R <= 0) return -1; 

  Deque_of_Space_families families_rs;  
  generate_rs_families(families_rs, R);

  //////////////////////////////////////////////////
  // We now have a deque of families (families_rs) whose invariants r & s agree.
  // Now look for those families for which also p1 agrees!

  families_rs.compute_KS_invariants();
  Deque_of_Space_families families_he;
  Deque_of_Space_families::filter(families_rs,families_he,Space::compareHomotopyType);

  /*
  Deque_of_Space_families families_rsp;  
  printf("Looking for pairs whose invariants (r,s,p) agree ...\n");
  Deque_of_Space_families::filter(families_rs,families_rsp,Space::compareBasicType);l
  families_rsp.sort_and_count_families();
  families_rsp.print("list_rsp_human.txt");
  */

  // print first list of families to file:
  // (do this AFTER s22-values have been computed)
  //families_rs.sort_and_count_families(); //stable_sort(families_rs.begin(),families_rs.end());
  //families_rs.print("list_rs_human.txt");

  //print second family to file:
  families_he.sort_and_count_families();
  families_he.print("list_he_human.txt");

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

void test_sign(void)
{
  printf("%d\n",sign(rational<long long>(1,8)));
  printf("%d\n",sign(rational<long long>(-1,4)));
  printf("%d\n",sign(rational<long long>(4,2)));
  printf("%d\n",sign(rational<long long>(-5,2)));
  printf("%d\n",sign(rational<long long>(1,2)));
  printf("%d\n",sign(+3));
  printf("%d\n",sign(0));
}
