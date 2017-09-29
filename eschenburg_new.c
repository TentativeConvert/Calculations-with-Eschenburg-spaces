/* compile with: 

   g++ -std=c++11 -O3 eschenburg_new.c -o e_new
   g++ -std=c++11 eschenburg_new.c -o e_new

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
  /*
    "basic invariants":  r, s, p_1
    r, s, s_22      coincide  =>  spaces homotopy equivalent
    r, s, s_22, p_1 coincide  =>  spaces tangentially homotopy equivalent
    If values of s_2 differ, then spaces are not homeomorphic.
  */

  // by [CEZ06, Lemma 1.4] we have "standard parametrization" by 
  //    k1 >= k2 > l1 >= l2 >= 0;
  // by [CEZ06, proof of Prop. 1.7],  
  //    r(k1,k2,l1,l2) = R
  // implies that moreover R >= k1.
  //
  // in the code, k1 is called k0 or k[0]
  //              k2 is called k1 or k[1]
  //              k3 is called k2 or k[2]
  // etc.

  long R;
  printf("\n Maximum value of |r|: ");
  scanf("%ld",&R); 

  vector< struct Space> spaces [R+1];

  // (1) Find coprime pairs (n, d), such that  R >= d >= n > 0 using Farey algorithm.
  //     We will take 
  //        d = k1 - l2
  //        n = k2 - l1
  //
  long old_n = 0;       // first pair in Farey sequence is 0/1
  long old_d = 1;       //   discard this pair -- it's only needed to start the algorithm
  long n = 1;          // second pair in Farey sequence is 1/R
  //**long R_ = R;
  long R_ = (long)std::sqrt(R); // experimental 
  long d = R_;
  while (n <= R_)
    {
       if(d == 101) // pick a prime number here
      	printf("%3ld/%ld\n", n, d);
      
      // (2) For the current Farey pair (n,d), 
      //     find values of k2 ...
      if( n*(d+n) <= R ) // n*(d+n) is a lower bound for |r|
	{
	  for(long k2_ = 1; k2_ <= d; ++k2_)
	    {
	      if(gcd(k2_,d) == 1) 
		{
		  for(long k2 = k2_; k2 <= R; k2+=d)
		    {
		      // (3) For current pair (n,d) and current k2, 
		      //     find values of k1
		      long min_k1 = max(k2, d);
		      long max_k1 = min(R, d - n + k2);
		      long max_k1_ = min(min_k1 + n - 1, max_k1);
		      for(long k1_ = min_k1; k1_ <= max_k1_; ++k1_)
			{
			if(gcd(k1_,n) == 1)
			  {
			  for(long k1 = k1_; k1 <= max_k1; k1+=n)
			    {
			      long l1 = k2 - n;
			      long l2 = k1 - d;
			      
			      // Check that |r| <= R:
			      long mr = k1*(k1-l1)+(k2-l2)*(k1+k2-l1); // see [CEZ06, proof of Prop. 1.7]
			      if (mr > R) break; // note that r is monotonous in k1
			      
			      // Check remaining conditions for (k1,k2,k3 | l1,l2,l3)
			      // to define an Eschenburg space 
			      // (conditions (1.1) of [CEZ06]):
			      if (gcd(k1, k2 - l2)      != 1) continue;
			      if (gcd(k1 - l1, k2 - l2) != 1) continue;
			      if (gcd(k1 - l1, k2)      != 1) continue;
			      
			      //if (d > (long)std::sqrt(R))
			      //	printf(" %ld/%ld ---!--- \n",n,d);
				      
			      // If we get this far, we've found an actual 
			      // positively curved Eschenburg space with |r| < R.
			      // Add it to our list:
			      struct Space new_space = {{k1, k2, l1 + l2 - k1 - k2}, {l1, l2, 0}, mr,0,0};
			      new_space.compute_s_and_p1();
			      /* fprintf(file_list_human,"%ld/%ld: ",n,d);
				 new_space.print_for_human(file_list_human);
				 fprintf(file_list_human,"\n");*/
			      spaces[new_space.mr].push_back(new_space);
			      // printf("%ld, %ld, %ld, %ld, %ld, %ld\n", gcd(k1-l1, k2-l2), gcd(k1-l1, k2), gcd(k1, k2-l1), gcd(k1-l2, k2-l1), gcd(k1-l2, k2), gcd(k1, k2-l2));
			    }
			  }
			}
		    }
		}
	    }
	}
      // Generate next Farey pair:
      long k = (long)((R_ + old_d) / d);
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
