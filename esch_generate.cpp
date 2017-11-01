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
  Feedback feedback;
  feedback.start(100);

  std::size_t c_spaces = 0; // counter
  struct miniSpace {
    INT_p d;
    INT_p n;
    INT_p k1;
    INT_p k2;
    INT_R r() const { return -(k1*d + n*d + n*k2); };
    INT_R unreduced_s() const { return -(k1*k2*(n+d)); };
  };
  deque< deque< struct miniSpace > > all_spaces((INT_R)((R+1)/2)); 
  // deque of deque of spaces that we find: one deque for each value of |r| <= R
  // allocated dynamically, so size is only limited by OS/hardware 
  // (see https://stackoverflow.com/a/216731/3611932)  
  
  // Step (a)
  INT_p D = (INT_p)(sqrt(R-3/4) - 1/2) + 1; 
                         // + 1 because result of sqrt is rounded down
                         // sqrt(integer) returns double, which has 53 bit significand
  INT_p old_n = 0;	 // first pair in Farey sequence is 0/1
  INT_p old_d = 1;	 // -- only used to start the algorithm
  INT_p n = 1;		 // second pair in Farey sequence is 1/D
  INT_p d = D;		 //
  while (n <= D)
    {
      if(d == 101) feedback.update_percent((int)(n*100/101));
      // Step (b.1) 
      INT_p K1	= (INT_p)((R-n*n)/d - n) + 1 ; // + 1 to guard against rounding errors
      INT_p K1_ = min(d+n-1, K1);
      for(INT_p k1_ = d; k1_ <= K1_; ++k1_)
	{
	  if(gcd(k1_,n) != 1) continue; // i.e. try next k1_
	  
	  // Step (b.2)
	  for(INT_p k1 = k1_; k1 <= K1; k1+=n)
	    {
	      // Step (c.1)
	      INT_p K2 = min((INT_p)((R-k1*d)/n - d + 1), k1); // + 1 to guard against rounding errors 
	      INT_p K2_ = min(k1+n-1, K2);
	      for(INT_p k2_ = k1+n-d; k2_ <= K2_; ++k2_)
		{
		  if(gcd(k2_,d) != 1) continue; // i.e. try next k2_
		  
		  // Step (c.2)
		  for(INT_p k2 = k2_; k2 <= K2; k2+=d)
		    {
		      // Step (d)
		      INT_p l1 = k2 - n;
		      INT_p l2 = k1 - d;

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
		      miniSpace new_mini = {d, n, k1, k2};
		      //  0 < mr <= R   and   mr is odd
		      //  0 <= (mr-1)/2 <= (R-1)/2
		      all_spaces[(INT_R)((mr-1)/2)].push_back(new_mini);
		    }
		}
	    }
	}
      // Step (a) continued -- generate next Farey pair:
      INT_p k = (INT_p)((D + old_d) / d);
      INT_p new_n = k*n - old_n;
      INT_p new_d = k*d - old_d;
      old_n = n;  
      old_d = d;
      n = new_n;  
      d = new_d;
    }
  feedback.finish();
  printf("Found %ld spaces in this range.\n", c_spaces);
 
  //////////////////////////////////////////////////
  // List of spaces (all_spaces) is now complete.
  // Now look for tuples of spaces whose invariants r & s agree!
  printf("\nLooking for tuples of spaces whose invariants r & s agree ...\n");
 
  std::size_t counter_distinct_rs_values = 0;
  std::size_t counter_singletons = 0;
  feedback.start((std::size_t)(R+1)/2);

  for(INT_R hmr = 0; hmr < (R+1)/2; ++hmr){  //hmr = "half minus r" (abgerundet)
    feedback.update((std::size_t)hmr);
    //------------------------------------------------
    // Sort spaces in all_spaces[hmr] by their s-invariant.
    //
    // To speed up sorting, we make a list (vector) of pairs (index, s),
    // where s is the s-invariant of the space all_spaces[hmr][index].
    // Then we only sort this list of pairs.
    struct i_s_pair {
      std::size_t i; INT_R abs_s;  
      bool operator<(const i_s_pair& otherpair) const {return (abs_s < otherpair.abs_s);}
    };
    vector< struct i_s_pair > i_s_pairs(all_spaces[hmr].size());
    INT_R mr = 2*hmr+1;    
    for(std::size_t i = 0; i < all_spaces[hmr].size(); ++i){
      i_s_pairs[i].i = i;
      i_s_pairs[i].abs_s = abs(signed_mod(all_spaces[hmr][i].unreduced_s(), mr));
    }
    std::sort(i_s_pairs.begin(),i_s_pairs.end());
    // Our list of paris (i_s_pairs) is now sorted.
    // Now find spaces where s-values match.
    for(std::size_t i1 = 0; i1 < i_s_pairs.size(); )// i1 is incremented indirectly via i2
      {
	std::size_t i2 = i1+1;
	while (i2 < i_s_pairs.size() 
	       && i_s_pairs[i1].abs_s == i_s_pairs[i2].abs_s)
	  ++i2;
	++counter_distinct_rs_values;
	if(i2 > i1+1)
	  {
	    //all_spaces with indexes i1,...,i2 define spaces with the same invariants
            SpaceTuple new_tuple;
	    for(std::size_t j = i1; j < i2; ++j)
	      {
		miniSpace& e = all_spaces[hmr][i_s_pairs[j].i];
		//Space E;
		//E.setParameters({e.k1, e.k2, -e.n-e.d},{e.k2-e.n, e.k1-e.d, 0});
		Space E({e.k1, e.k2, -e.n-e.d},{e.k2-e.n, e.k1-e.d, 0});
		new_tuple.push_back(E);
	      }
	    //tuples_rs.push_back(new_tuple); //xxx
	    this->push_back(new_tuple);
	  }
	else
	  ++counter_singletons;
	i1 = i2;
      }
    all_spaces[hmr].clear();  // free up memory space!
  }
  feedback.finish();
  printf("Found %lld distinct values of (r,s) in this range;\n", (long long)counter_distinct_rs_values);
  printf("      %lld values occur exactly once.\n\n", (long long)counter_singletons);
}

