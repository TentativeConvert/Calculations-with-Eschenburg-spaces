/* compile with: 

   g++ -std=c++11 -O3 esch-gen-2.c -o esch-gen-2
   g++ -std=c++11 esch-gen-2.c -o esch-gen-2

   Flags:
   -std=c++11 specifies the "language version" (I think)
   -O3        enables optimisation for speed 
*/
#include <cstdio>
#include <vector>
#include <cmath>
#include <boost/math/common_factor.hpp>

//////////////////////////////////////////////////
// Auxiliary mathematics:
using boost::math::gcd;
using std::min;
using std::max;
using std::sqrt;
using std::vector;
long absolute (long a){ return a >= 0 ? a : -a; }
long signed_mod (long a, long base)
{
  long remainder = a % base;
  if(2*remainder > base)  // remainder > base/2, in "integer language"
    remainder -= base;
  else if(2*remainder <= -base) // remainder < -base/2, in "integer language"
    remainder += base;
  return remainder;
}
long absolute_mod (long a, long base)
{
  long remainder = a % base;
  if(remainder < 0)
    remainder += absolute(base);
  return remainder;
}

//////////////////////////////////////////////////
//  Spaces:
struct Space {
  long k[3];
  long l[3];
  long mr; // mr = -r = "minus r"
  long s;
  long p1;
  void print(void){
    printf("[%4ld,%4ld,%4ld |%4ld,%4ld,%4ld] -> r = %4ld, s = %4ld, p1 = %4ld\n", k[0],k[1],k[2],l[0],l[1],l[2],mr, s, p1);
  }
  void print_for_human(FILE* file){
    fprintf(file, "[%4ld,%4ld,%4ld |%4ld,%4ld,%4ld] -> r = %4ld, s = %4ld, p1 = %4ld", k[0],k[1],k[2],l[0],l[1],l[2],mr, s, p1);
  }
  void print_for_maple(FILE* file){
    fprintf(file, "[%ld,%ld,%ld,%ld,%ld,%ld]",k[0],k[1],k[2],l[0],l[1],l[2]);
  }
  void compute_s_and_p1(void){
    int sigma1_k = k[0] + k[1] + k[2];
    int sigma2_k = k[0]*k[1] + k[0]*k[2] + k[1]*k[2];
    int sigma3_k = k[0]*k[1]*k[2];
    int sigma1_l = l[0] + l[1] + l[2];
    int sigma2_l = l[0]*l[1] + l[0]*l[2] + l[1]*l[2];
    int sigma3_l = l[0]*l[1]*l[2];
    s = signed_mod(sigma3_k, mr);  // note: sigma3_l = 0
    p1 = absolute_mod(2*sigma1_k*sigma1_k - 6*sigma2_k, mr);
  }
};


//////////////////////////////////////////////////
// Main routine:

main(){
  double epsilon = 0.1; // used as "safety buffer" against rounding erros
  
  long R;
  printf("\n Maximum value of |r|: ");
  scanf("%ld",&R); 

  vector< struct Space> spaces [R+1]; // vector of lists of spaces that we find, one list for each value of |r| <= R

  // Step (a)
  long D = (long)(sqrt(R-3/4) - 1/2 + epsilon);
  long old_n = 0;       // first pair in Farey sequence is 0/1
  long old_d = 1;       // -- only used to start the algorithm
  long n = 1;           // second pair in Farey sequence is 1/D
  long d = D;           //
  while (n <= D)
    {
      if(d == 101) // progress feedback for user
      	printf("%ld%%  ", n);  // will display 1%  2%  3%  ...
      
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
		      struct Space new_space = {{k1, k2, l1 + l2 - k1 - k2}, {l1, l2, 0}, mr,0,0};
		      new_space.compute_s_and_p1();
		      spaces[new_space.mr].push_back(new_space);
		      // printf("%ld, %ld, %ld, %ld, %ld, %ld\n", gcd(k1-l1, k2-l2), gcd(k1-l1, k2), gcd(k1, k2-l1), gcd(k1-l2, k2-l1), gcd(k1-l2, k2), gcd(k1, k2-l2));
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
  printf("\n");


  //////////////////////////////////////////////////
  // List of spaces is now complete.
  // Now find pairs whose basic invariants agree 
  // and write them to a file:

  FILE *file_list_maple = fopen("list_basic_maple.txt", "w");  // r, s, p1 agree
  FILE *file_list_human = fopen("list_basic_human.txt", "w");  // r, s, p1 agree
  if ((file_list_maple == NULL) || (file_list_human == NULL))
    printf("Error opening file!\n");
    
  int c_basic = 0; // counts pairs whose basic invariants agree
  for(int mr = 0; mr <= R; mr++){ 
    for(int i1 = 0; i1 < spaces[mr].size(); ++i1){
      for(int i2 = i1+1; i2 < spaces[mr].size(); ++i2){
	struct Space E1 = spaces[mr][i1];
	struct Space E2 = spaces[mr][i2];
	if ( (E1.s == E2.s || E1.s == -E2.s) && E1.p1 == E2.p1 ) 
	  {
	    ++c_basic;
	    printf("Found pair %d.\n", c_basic);
	    // print for maple:
	    fprintf(file_list_maple,"[");
	    E1.print_for_maple(file_list_maple);
	    fprintf(file_list_maple,",");
	    E2.print_for_maple(file_list_maple);
	    fprintf(file_list_maple,"]\n");
	    // print for human:
	    fprintf(file_list_human,"\n Pair %d:\n ", c_basic);
	    E1.print_for_human(file_list_human);
	    fprintf(file_list_human,"\n ");
	    E2.print_for_human(file_list_human);
	    fprintf(file_list_human,"\n");
	  }
      }
    }
  }

  fclose(file_list_maple);
  fclose(file_list_human);
}
