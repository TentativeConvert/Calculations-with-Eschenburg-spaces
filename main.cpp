/* compile with: 
   make
*/
#include <cstdio>
#include <cstring>
using std::strcat;
#include <vector>
using std::vector;
#include <deque>
using std::deque;
#include <array>
using std::array;
#include <algorithm>
using std::sort;
#include <chrono>
using namespace std::chrono;
//////////////////////////////////////////////////
// Auxiliary mathematics:
#include <cmath>
using std::min;
using std::max;
using std::sqrt;
using std::abs;
#include <boost/math/common_factor.hpp>
using boost::math::gcd;
#include "aux_math.h"
#include "eschenburg.h"
//#include "space_mini.h"

//////////////////////////////////////////////////
// Tests (implementation further below)
void test_rationals(void);

//////////////////////////////////////////////////
// Main routine:

main(){
  //test_rationals();

  /*class Space E1;
  E1.setParameters({0,0,0},{0,0,0});
  deque<int> my_deque;
  int a = 5;
  my_deque.push_back(a);
  my_deque.push_back(5);*/

  double epsilon = 0.1; // used as "safety buffer" against rounding erros
  
  long R;
  printf("\nMaximum value of |r|: ");
  if(scanf("%ld",&R) != 1 || R <= 0) return -1; 
  printf("Looking for Eschenburg spaces with |r| <= %ld ... \n", R);

  system_clock::time_point t1 = system_clock::now();/*timer*/

  long c_spaces = 0; // counter
  struct miniSpace {
    long d;
    long n;
    long k1;
    long k2;
    long r() const { return -(k1*d + n*d + n*k2); };
    long unreduced_s() const { return -(k1*k2*(n+d)); };
  };
  deque< deque< struct miniSpace > > all_spaces((long)((R+1)/2)); 
  // array of lists ("vectors") of spaces that we find, one list for each value of |r| <= R
  // allocated dynamically, so size is only limited by OS/hardware 
  // (see https://stackoverflow.com/a/216731/3611932)  

  // Step (a)
  long D = (long)(sqrt(R-3/4) - 1/2 + epsilon);
  long old_n = 0;       // first pair in Farey sequence is 0/1
  long old_d = 1;       // -- only used to start the algorithm
  long n = 1;           // second pair in Farey sequence is 1/D
  long d = D;           //
  while (n <= D)
    {
      if(d == 101) // progress feedback for user
	{
	  printf(" %3d%%\r", (int)(n*100/101));  // will display 1%  2%  3%  ...
	  fflush(stdout);
	}
      // Step (b.1) 
      long K1  = (long)((R-n*n)/d - n + epsilon);
      long K1_ = min(d+n-1, K1);
      for(long k1_ = d; k1_ <= K1_; ++k1_)
	{
	  if(gcd(k1_,n) != 1) continue; // i.e. try next k1_
	  
	  // Step (b.2)
	  for(long k1 = k1_; k1 <= K1; k1+=n)
	    {
	      // Step (c.1)
	      long K2 = min((long)((R-k1*d)/n - d + epsilon), k1);
	      long K2_ = min(k1+n-1, K2);
	      for(long k2_ = k1+n-d; k2_ <= K2_; ++k2_)
		{
		  if(gcd(k2_,d) != 1) continue; // i.e. try next k2_
		  
		  // Step (c.2)
		  for(long k2 = k2_; k2 <= K2; k2+=d)
		    {
		      // Step (d)
		      long l1 = k2 - n;
		      long l2 = k1 - d;

		      // Check conditions (2'd):
		      if (gcd(k2,    k1-l1) != 1) continue;
		      if (gcd(k1,    k2-l2) != 1) continue;
		      if (gcd(k1-l1, k2-l2) != 1) continue;

		      long mr = k1*d + n*d + n*k2; // "mr" = "minus r" = |r|, since r always negative in our parametrization

		      // Check that indeed |r| <= R  
		      // -- this should be automatic except for rounding errors in boundary cases 
		      if (mr > R) break; // note that r is monotonous in k2

		      // If we get this far, we've found a positively curved Eschenburg space with |r| < R.
		      // Add it to our list:
		      ++c_spaces;
		      miniSpace new_mini = {d, n, k1, k2};
		      //  0 < mr <= R   and   mr is odd
		      //  0 <= (mr-1)/2 <= (R-1)/2
		      all_spaces[(long)((mr-1)/2)].push_back(new_mini);
		    }
		}
	    }
	}
      // Step (a) continued -- generate next Farey pair:
      long k = (long)((D + old_d) / d);
      long new_n = k*n - old_n;
      long new_d = k*d - old_d;
      old_n = n;  
      old_d = d;
      n = new_n;  
      d = new_d;
    }
  printf(" 100%%");
  system_clock::time_point t2 = system_clock::now();/*timer*/
  float duration = (float)duration_cast<milliseconds>( t2 - t1 ).count()/1000;/*timer*/
  printf("\nFound %ld spaces in this range (time: %.2f s).\n", c_spaces,duration);

  //////////////////////////////////////////////////
  // List of spaces (all_spaces) is now complete.
  // Now look for families of spaces whose invariants r & s agree!
  printf("\nLooking for candidate spaces whose invariants r & s agree ...\n");

  deque< deque< class Space > > families_rs;
  // a deque of families (deques) of spaces whose invariants r & s agree
 
  long c_rs_spaces = 0;
  long feedback_step_size = max((long)((R+1)/2/100),(long)1);

  for(long hmr = 0; hmr < (R+1)/2; ++hmr){  //hmr = "half minus r" (abgerundet)
    if(hmr % feedback_step_size == 0) // feedback for user
      {
	printf(" %3d%%\r",(int)(100*hmr*2/(R+1)));
	fflush(stdout);
      }
    //------------------------------------------------
    // Sort spaces in all_spaces[hmr] by their s-invariant.
    //
    // To speed up sorting, we make a list (vector) of pairs (index, s),
    // where s is the s-invariant of the space all_spaces[hmr][index].
    // Then we only sort this list of pairs.
    struct i_s_pair {
      long i; long abs_s;  
      bool operator<(const i_s_pair& otherpair) const {return (abs_s < otherpair.abs_s);}
    };
    vector< struct i_s_pair > i_s_pairs(all_spaces[hmr].size());
    long mr = 2*hmr+1;    
    for(long i = 0; i < all_spaces[hmr].size(); ++i){
      i_s_pairs[i].i = i;
      i_s_pairs[i].abs_s = abs(signed_mod(all_spaces[hmr][i].unreduced_s(), mr));
    }
    sort(i_s_pairs.begin(),i_s_pairs.end());
    // Our list of paris (i_s_pairs) is now sorted.
    // Now find spaces where s-values match.
    for(long i1 = 0; i1 < i_s_pairs.size(); )// i1 is incremented indirectly via i2
      {
	long i2 = i1+1;
	while (i2 < i_s_pairs.size() 
	       && i_s_pairs[i1].abs_s == i_s_pairs[i2].abs_s)
	  ++i2;
	if(i2 > i1+1)
	  {
	    //all_spaces with indexes i1,...,i2 define spaces with the same invariants
            deque<class Space> new_deque;
	    for(long j = i1; j < i2; ++j)
	      {
                Space E;
		miniSpace& e = all_spaces[hmr][i_s_pairs[j].i];
		E.setParameters({e.k1, e.k2, -e.n-e.d},{e.k2-e.n, e.k1-e.d, 0});
		new_deque.push_back(E);
		++c_rs_spaces;
	      }
	    families_rs.push_back(new_deque);
	  }
	i1 = i2;
      }
    all_spaces[hmr].clear();  // free up memory space!
  }
  printf(" 100%%");
  t2 = system_clock::now();/*timer*/
  duration = (float)duration_cast<milliseconds>( t2 - t1 ).count()/1000;/*timer*/
  printf("\nFound %ld spaces with non-unique values of (r,s) in this range (total time: %.2f s).\n\n", c_rs_spaces,duration);

  //////////////////////////////////////////////////
  // We now have a deque of families (families_rs) whose invariants r & s agree.
  // Now look for those families for which also p1 agrees!

  printf("Looking for pairs whose invariants (r,s,p1) agree ...\n");
  FILE *file_maple = fopen("list_basic_pairs_maple.txt", "w");  // r, s, p1 agree
  FILE *file_human = fopen("list_basic_pairs_human.txt", "w");  // r, s, p1 agree
  FILE *file_maple_3 = fopen("list_basic_triples_maple.txt", "w");  
  FILE *file_human_3 = fopen("list_basic_triples_human.txt", "w");  
  if (file_maple == NULL || file_human == NULL || file_human_3 == NULL || file_maple_3 == NULL)
    printf("Error opening file!\n");

  long c_pairs = 0;
  long c_triples = 0;
  deque< deque< class Space > > families_rsp;  

  long c_rsp_spaces = 0;
  feedback_step_size = max((long)(families_rs.size()/100),(long)1);
  for(long i = 0; i < families_rs.size(); ++i)
    {
      if(i % feedback_step_size == 0) // feedback for user
	{
	  printf(" %3d%%\r",(int)((float)i/feedback_step_size));
	  fflush(stdout);
	}
      struct mycomparator {
	bool operator() (const Space& E1, const Space& E2) { return (E1.p1() < E2.p1());}
	//bool operator() (Space E1, Space E2) { return (E1.s22() < E2.s22());}  // <--- will need to pass by value, not by reference for s22()
      } p1_comparator;
      sort(families_rs[i].begin(),families_rs[i].end(),p1_comparator);
      for(long i1 = 0; i1 < families_rs[i].size(); ++i1)
      {
	for(long i2 = 1; i1+i2 < families_rs[i].size(); ++i2)
	  {
	    Space E1 = families_rs[i][i1];
	    Space E2 = families_rs[i][i1+i2];
	    if ( E1.p1() == E2.p1() )
	      {
		++c_pairs;
		// print for maple:
		fprintf(file_maple,"[");
		E1.print_for_maple(file_maple);
		fprintf(file_maple,",");
		E2.print_for_maple(file_maple);
		fprintf(file_maple,"]\n");
		// print for human:
		fprintf(file_human,"\n Pair %ld:\n ", c_pairs);
		E1.print_for_human(file_human);
		fprintf(file_human,"\n ");
		E2.print_for_human(file_human);
		fprintf(file_human,"\n");
		if(i2 > 1) {
		  Space E3 = families_rs[i][i1+1];
		  ++c_triples;	      
		  // print for maple:
		  fprintf(file_maple_3,"[");
		  E1.print_for_maple(file_maple_3);
		  fprintf(file_maple_3,",");
		  E3.print_for_maple(file_maple_3);
		  fprintf(file_maple_3,",");
		  E2.print_for_maple(file_maple_3);
		  fprintf(file_maple_3,"]\n");
		  // print for human:
		  fprintf(file_human_3,"\n Triple %ld:\n ", c_triples);
		  fprintf(file_human_3,"(i1 = %ld, i2 = %ld)\n ",i1,i2);
		  E1.print_for_human(file_human_3);
		  fprintf(file_human_3,"\n ");
		  E3.print_for_human(file_human_3);
		  fprintf(file_human_3,"\n ");
		  E2.print_for_human(file_human_3);
		  fprintf(file_human_3,"\n");
		}
	      }
	    else break;  // can break here since list is sorted
	}
      }
    }
  printf(" 100%%");
  t2 = system_clock::now();/*timer*/
  duration = (float)duration_cast<milliseconds>( t2 - t1 ).count()/1000;/*timer*/
  printf("\nFound %ld pairs and %ld triples in this range (total time: %.2f s).\n\n", c_pairs, c_triples,duration);

  fclose(file_maple);
  fclose(file_human);
  fclose(file_maple_3);
  fclose(file_human_3);
  return 0;
}


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
