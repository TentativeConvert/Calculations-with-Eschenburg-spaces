#include "esch_space.h"
#include "esch_tuples.h"
using std::deque;
#include <vector>
using std::vector;
#include "aux_feedback.h"
#include "aux_math.h"
#include <algorithm>
#include <cstddef> // for std::size_t


SpaceTupleList::SpaceTupleList(const INT_R& R)
{
  printf("\nLooking for Eschenburg spaces with |r| <= %ld ... \n", (long)R);
  this->indeterminacies=0;
  this->description="homotopy classes";

  Feedback feedback;
  feedback.start(100);

  std::size_t c_spaces = 0; // counter
  struct tinySpace {
    INT_P d;
    INT_P n;
    INT_P k1;
    INT_P k2;
    //INT_R r() const { return -(k1*d + n*d + n*k2); };
    INT_R unreduced_s() const { return -(k1*k2*(n+d)); };
    int_least8_t M1() const { return (int_least8_t)signed_mod(-k1-k2+n+d, 3); };
    int_least8_t M2() const { return (int_least8_t)absolute_mod(k1 + k2 -n-d + k1*k2 -k1*n-k1*d -k2*n-k2*d, 2); };
  };
  deque< deque< struct tinySpace > > all_spaces((INT_R)((R+1)/2)); 
  // deque of deque of spaces that we find: one deque for each value of |r| <= R
  // allocated dynamically, so size is only limited by OS/hardware 
  // (see https://stackoverflow.com/a/216731/3611932)  
  
  // Step (a)
  INT_P D = (INT_P)(sqrt(R-3/4) - 1/2) + 1; 
                         // + 1 because result of sqrt is rounded down
                         // sqrt(integer) returns double, which has 53 bit significand
  INT_P old_n = 0;	 // first pair in Farey sequence is 0/1
  INT_P old_d = 1;	 // -- only used to start the algorithm
  INT_P n = 1;		 // second pair in Farey sequence is 1/D
  INT_P d = D;		 //
  while (n <= D)
    {
      if(d == 101) feedback.update_percent((int)(n*100/101));
      // Step (b.1) 
      INT_P K1	= (INT_P)((R-n*n)/d - n) + 1 ; // + 1 to guard against rounding errors
      INT_P K1_ = min(d+n-1, K1);
      for(INT_P k1_ = d; k1_ <= K1_; ++k1_)
	{
	  if(gcd(k1_,n) != 1) continue; // i.e. try next k1_
	  
	  // Step (b.2)
	  for(INT_P k1 = k1_; k1 <= K1; k1+=n)
	    {
	      // Step (c.1)
	      INT_P K2 = min((INT_P)((R-k1*d)/n - d + 1), k1); // + 1 to guard against rounding errors 
	      INT_P K2_ = min(k1+n-1, K2);
	      for(INT_P k2_ = k1+n-d; k2_ <= K2_; ++k2_)
		{
		  if(gcd(k2_,d) != 1) continue; // i.e. try next k2_
		  
		  // Step (c.2)
		  for(INT_P k2 = k2_; k2 <= K2; k2+=d)
		    {
		      // Step (d)
		      INT_P l1 = k2 - n;
		      INT_P l2 = k1 - d;

		      // Check conditions (2'd):
		      if (gcd(k2,    k1-l1) != 1) continue;
		      if (gcd(k1,    k2-l2) != 1) continue;
		      if (gcd(k1-l1, k2-l2) != 1) continue;

		      INT_R mr = k1*d + n*d + n*k2; // "mr" = "minus r" = |r|, since r always negative in our parametrization

		      // Check that indeed |r| <= R  
		      // -- this should be automatic except for rounding errors in boundary cases 
		      if (mr > R) break; // note that r is monotonous in k2

		      // If we get this far, we've found a positively curved Eschenburg space with |r| < R.
		      // Add it to our list:
		      ++c_spaces;
		      tinySpace new_mini = {d, n, k1, k2};
		      //  0 < mr <= R   and   mr is odd
		      //  0 <= (mr-1)/2 <= (R-1)/2
		      all_spaces[(INT_R)((mr-1)/2)].push_back(new_mini);
		    }
		}
	    }
	}
      // Step (a) continued -- generate next Farey pair:
      INT_P k = (INT_P)((D + old_d) / d);
      INT_P new_n = k*n - old_n;
      INT_P new_d = k*d - old_d;
      old_n = n;  
      old_d = d;
      n = new_n;  
      d = new_d;
    }
  feedback.finish();
  printf(">> %9lld spaces in this range.\n", (long long)c_spaces);
 
  //////////////////////////////////////////////////
  // List of spaces (all_spaces) is now complete.
  // Now look for tuples of spaces whose invariants |r| & |s| agree!
  printf("\nLooking for spaces of same homotopy type ...\n");
 
  std::size_t counter_distinct_rs_values = 0;
  // std::size_t counter_singletons = 0;
  singletons = 0;
  feedback.start((std::size_t)(R+1)/2);

  for(INT_R hmr = 0; hmr < (R+1)/2; ++hmr){  //hmr = "half minus r" (abgerundet)
    feedback.update((std::size_t)hmr);
    //------------------------------------------------
    // Sort spaces in all_spaces[hmr] by their invariants |s| and |M1|.
    // To speed up sorting, turn the deque all_spaces[hmr] into a list r_spaces.
    struct miniSpace {
      // slightly higher memory usage than tinySpace: two more variables (s & M1)
      INT_P d;
      INT_P n;
      INT_P k1;
      INT_P k2;
      INT_R s;
      int_least8_t M1; // only values are +1, 0, -1
      int_least8_t M2; // only values are 1, 0
      bool operator<(const miniSpace& otherspace) const {
	if (abs(s) > abs(otherspace.s)) return false;
	if (abs(s) < abs(otherspace.s)) return true;
	if (abs(M1) > abs(otherspace.M1)) return false;
	if (abs(M1) < abs(otherspace.M1)) return true;
	if (sign(s)*sign(M1) > sign(otherspace.s)*sign(otherspace.M1)) return false;
	if (sign(s)*sign(M1) < sign(otherspace.s)*sign(otherspace.M1)) return true;
	if (M2 > otherspace.M2) return false;
	if (M2 < otherspace.M2) return true;
	return false;
      }
    };
    vector< struct miniSpace > r_spaces(all_spaces[hmr].size());
    INT_R mr = 2*hmr+1;    
    for(std::size_t i = 0; i < all_spaces[hmr].size(); ++i){
      const tinySpace& E = all_spaces[hmr][i];
      r_spaces[i].d =  E.d;
      r_spaces[i].n =  E.n;
      r_spaces[i].k1 = E.k1;
      r_spaces[i].k2 = E.k2;
      r_spaces[i].s =  signed_mod(E.unreduced_s(),mr);
      r_spaces[i].M1 = E.M1();
      r_spaces[i].M2 = E.M2();
    }
    std::sort(r_spaces.begin(),r_spaces.end());
    // Our list of paris (r_spaces) is now sorted.
    // Now find spaces where s-values match.
    for(std::size_t i1 = 0; i1 < r_spaces.size(); )// i1 is incremented indirectly via i2
      {
	std::size_t i2 = i1+1;
	while (i2 < r_spaces.size() 
	       && abs(r_spaces[i1].s) == abs(r_spaces[i2].s)
	       && abs(r_spaces[i1].M1) == abs(r_spaces[i2].M1)
	       && sign(r_spaces[i1].s)*sign(r_spaces[i1].M1) == sign(r_spaces[i2].s)*sign(r_spaces[i2].M1)
	       && r_spaces[i1].M2 == r_spaces[i2].M2
	       )
	  ++i2;
	++counter_distinct_rs_values;
	if(i2 > i1+1)
	  {
	    //r_spaces[i1], ..., r_spaces[i2] are spaces with the same invariants
            SpaceTuple new_tuple;
	    for(std::size_t j = i1; j < i2; ++j)
	      {
		miniSpace& e = r_spaces[j];
		Space E({e.k1, e.k2, -e.n-e.d},{e.k2-e.n, e.k1-e.d, 0});
		new_tuple.push_back(E);
	      }
	    this->push_back(new_tuple);
	  }
	else
	  ++singletons;
	i1 = i2;
      }
    all_spaces[hmr].clear();  // free up memory space!
  }
  feedback.finish();
  printf(">> %9lld distinct values of (r,s) in this range;\n", (long long)counter_distinct_rs_values);
  printf(">> %9lld values occur exactly once.\n\n", (long long)singletons);
}

