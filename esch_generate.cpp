#include "esch_space.h"
#include "esch_tuples.h"
using std::deque;
#include <vector>
using std::vector;
#include "aux_feedback.h"
#include "aux_math.h"
#include <algorithm>
using std::sort;


SpaceTupleList::SpaceTupleList(const long& R)
{
  printf("\nLooking for Eschenburg spaces with |r| <= %ld ... \n", R);
  Feedback feedback;
  feedback.start(100);

  double epsilon = 0.1; // used as "safety buffer" against rounding erros
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
  // deque of deque of spaces that we find: one deque for each value of |r| <= R
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
      if(d == 101) feedback.update_percent((int)(n*100/101));
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
  feedback.finish();
  printf("Found %ld spaces in this range.\n", c_spaces);
 
  //////////////////////////////////////////////////
  // List of spaces (all_spaces) is now complete.
  // Now look for tuples of spaces whose invariants r & s agree!
  printf("\nLooking for tuples of spaces whose invariants r & s agree ...\n");
 
  long c_rs_spaces = 0;
  feedback.start((size_t)(R+1)/2);

  for(long hmr = 0; hmr < (R+1)/2; ++hmr){  //hmr = "half minus r" (abgerundet)
    feedback.update((size_t)hmr);
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
            SpaceTuple new_tuple;
	    for(long j = i1; j < i2; ++j)
	      {
		miniSpace& e = all_spaces[hmr][i_s_pairs[j].i];
		//Space E;
		//E.setParameters({e.k1, e.k2, -e.n-e.d},{e.k2-e.n, e.k1-e.d, 0});
		Space E({e.k1, e.k2, -e.n-e.d},{e.k2-e.n, e.k1-e.d, 0});
		new_tuple.push_back(E);
		++c_rs_spaces;
	      }
	    //tuples_rs.push_back(new_tuple); //xxx
	    this->push_back(new_tuple);
	  }
	i1 = i2;
      }
    all_spaces[hmr].clear();  // free up memory space!
  }
  feedback.finish();
  printf("Found %ld spaces with non-unique values of (r,s) in this range.\n\n", c_rs_spaces);
}

