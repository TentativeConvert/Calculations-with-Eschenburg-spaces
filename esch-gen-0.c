/* compile with: 

   g++ -std=c++11 -O3 eschenburg.c

   Flags:
   -std=c++11 specifies the "language version" (I think)
   -O3        enables optimisation for speed 
*/
#include <cstdio>
#include <vector>
#include <cmath>
#include <boost/rational.hpp>
#include <boost/math/common_factor.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions/round.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>
//////////////////////////////////////////////////
// Auxiliary mathematics:
// should use templates or libary ...
long long square(long long a){ return a*a; }
int absolute (int a){ return a >= 0 ? a : -a; }
long long absolute (long long a){ return a >= 0 ? a : -a; }
int signed_mod (int a, int base)
// In general, sign of % is machine dependent, 
// see https://stackoverflow.com/a/4003287/3611932.
// For the computation of s mod r,
// I want modulus numbers to be "symmetric around zero".
// Input:   a  in ZZ, base in ZZ
// Output:  a' in (-base/2, base/2] in ZZ such that 
//          a' = a mod base
{
  int remainder = a % base;
  if(2*remainder > base)  // remainder > base/2, in "integer language"
    remainder -= base;
  else if(2*remainder <= -base) // remainder < -base/2, in "integer language"
    remainder += base;
  return remainder;
}
long long signed_mod (long long a, long long base)
{
  long long remainder = a % base;
  if(2*remainder > base)  // remainder > base/2, in "integer language"
    remainder -= base;
  else if(2*remainder <= -base) // remainder < -base/2, in "integer language"
    remainder += base;
  return remainder;
}
int absolute_mod (int a, int base)
{
  int remainder = a % base;
  if(remainder < 0)
    remainder += absolute(base);
  return remainder;
}
long long absolute_mod (long long a, long long base)
{
  long long remainder = a % base;
  if(remainder < 0)
    remainder += absolute(base);
  return remainder;
}
/*int gcd (int a, int b)
// Euclid's algorithm
// see https://codereview.stackexchange.com/a/39110
{
  int x;
  while (b)
    {
      x = a % b;
      a = b;
      b = x;
    }
  return absolute(a);
  }*/
boost::rational<long long> reduce_mod_ZZ(boost::rational<long long> q)
{
  // Input:   q  in QQ 
  // Output:  q' in (-1/2, 1/2] in QQ such that 
  //          q' = q in QQ/ZZ.
  int q_n = q.numerator();
  int q_d = q.denominator();
  if(q_d < 0) // probaby unnecessary, but depends on implementation of boost::rational
    {
      q_n = -q_n;
      q_d = -q_d;
    }
  q_n = signed_mod(q_n,q_d);  
  q.assign(q_n,q_d); // this should simplify the fraction
  return q;
}
//////////////////////////////////////////////////
// Lens space invariants

boost::rational<long long> lens_s2(int p, int param[4])
{
  boost::multiprecision::cpp_dec_float_100 a = 0;
  for(int k = 1; k < absolute(p); ++k)
    {
      boost::multiprecision::cpp_dec_float_100 s = 
	(cos(2*M_PI*k/absolute(p)) - 1)
	/sin(k*M_PI*param[0]/p)
	/sin(k*M_PI*param[1]/p)
	/sin(k*M_PI*param[2]/p)
	/sin(k*M_PI*param[3]/p);
      a += s;
    }
  ///printf("     16*p*45*s2 = %.6f  (This should be a integer)\n",a*45);
  std::cout << "    ";
  std::cout << std::setprecision(std::numeric_limits<boost::multiprecision::cpp_dec_float_100>::max_digits10)
	    << a*45 
	    << std::endl; // print a
  long long rounded_45_a = (long long)boost::math::round(a*45);
  printf("     round(...) = %lld\n",rounded_45_a);
  boost::rational<long long> s2(rounded_45_a,45*16*p);
  //// printf("     s2 = %lld/%lld\n",s2.numerator(),s2.denominator());
  return s2;
}


//////////////////////////////////////////////////
//  Spaces:
struct Space {
  int k[3];
  int l[3];
  int mr; // mr = -r = "minus r"
  int s;
  int p1;
  boost::rational<long long> s22;
  boost::rational<long long> s2;
  int good_row; // for "condition C"
  int good_col; // for "condition C"

  void print(FILE *file){
    printf(         "[%4d,%4d,%4d |%4d,%4d,%4d] --> r = %4d, s = %4d, p1 = %4d, s22 = %lld/%lld, s2 = %lld/%lld\n",
		    k[0],k[1],k[2],l[0],l[1],l[2],-mr,s,p1,s22.numerator(),s22.denominator(),s2.numerator(),s2.denominator());
    fprintf(  file, "[%4d,%4d,%4d |%4d,%4d,%4d] --> r = %4d, s = %4d, p1 = %4d, s22 = %lld/%lld, s2 = %lld/%lld\n",
	      k[0],k[1],k[2],l[0],l[1],l[2],-mr,s,p1,s22.numerator(),s22.denominator(),s2.numerator(),s2.denominator());
  }
  void compute_r_s_p1(void){
    int sigma1_k = k[0] + k[1] + k[2];
    int sigma2_k = k[0]*k[1] + k[0]*k[2] + k[1]*k[2];
    int sigma3_k = k[0]*k[1]*k[2];
    int sigma1_l = l[0] + l[1] + l[2];
    int sigma2_l = l[0]*l[1] + l[0]*l[2] + l[1]*l[2];
    int sigma3_l = l[0]*l[1]*l[2];
    mr = sigma2_l - sigma2_k;
    s = signed_mod(sigma3_k, mr);  // note: sigma3_l = 0
    p1 = absolute_mod(2*sigma1_k*sigma1_k - 6*sigma2_k, mr);
  }
  bool is_positively_curved_space(void) {
    // Is it an Eschenburg space? 
    // -- Test conditions of [CEZ06] (1.1):
    if(boost::math::gcd(k[0] - l[0], k[1] - l[1]) > 1) return false;
    if(boost::math::gcd(k[0] - l[0], k[1] - l[2]) > 1) return false;
    if(boost::math::gcd(k[0] - l[2], k[1] - l[0]) > 1) return false;
    if(boost::math::gcd(k[0] - l[1], k[1] - l[0]) > 1) return false;
    if(boost::math::gcd(k[0] - l[1], k[1] - l[2]) > 1) return false;
    if(boost::math::gcd(k[0] - l[2], k[1] - l[1]) > 1) return false;
    // Is the space positively curved? 
    // -- Test conditions of [CEZ06] (1.2):
    int min_l = std::min(l[0],std::min(l[1],l[2]));
    if (k[0] == min_l) return false;
    if (k[1] == min_l) return false;
    if (k[2] == min_l) return false;
    int max_l = std::max(l[0],std::max(l[1],l[2]));
    if (k[0] == max_l) return false;
    if (k[1] == max_l) return false;
    if (k[2] == max_l) return false;
    return true;
  }
  int gcdA(int i,int j,int ii,int jj){ 
    // auxiliary function for find_good_row_or_col
    // should be declared private, so I probably ought to make Space a proper class, not just a struct
     return boost::math::gcd(k[i]-l[j],k[ii]-l[jj]); 
  } 
  void find_good_col_or_row(void)
  {
    // see "condition C" in [CEZ06]
    good_row = -1;
    good_col = -1;
    for(int c = 0; c < 3; ++c){
      if(gcdA(0,c,1,c) == 1 && gcdA(0,c,2,c) == 1 && gcdA(1,c,2,c) == 1)
	{
	  good_col = c;
	  return;
	}
    }
    for(int r = 0; r < 3; ++r){
      if(gcdA(r,0,r,1) == 1 && gcdA(r,0,r,2) == 1 && gcdA(r,1,r,2) == 1)
	{
	  good_row = r;
	  return;
	}
    }
  }
  void compute_s2_col(int j)
  {
    // [CEZ06] (2.1)
    // Again, this should really be a private function.
    // The "public" function is compute_s22 below.
    //
    //     s2  = (q-2)/d + SUM, 
    //     s22 = 2|r|s2
    //
    int jp1 = absolute_mod(j+1,2);
    printf("I'm using column %d: ",j+1);
    printf("(%d,%d,%d)\n",k[0]-l[j],k[1]-l[j],k[2]-l[j]);
    //// printf("jp1 = %d",jp1+1);
    long long q = 
      square(k[0]-l[j])   + square(k[1]-l[j])   + square(k[2]-l[j]) +
      square(k[0]-l[jp1]) + square(k[1]-l[jp1]) + square(k[2]-l[jp1])
      - square(l[j]-l[jp1]);
    long long d = 16*3*(long long)(-mr)*(k[0]-l[j])*(k[1]-l[j])*(k[2]-l[j]);
    s2.assign(q-2,d);
    printf("q = %lld, r = %lld, d = %lld\n", q, (long long)(-mr), d);
    printf("=> (q-2)/d = %lld/%lld\n",s2.numerator(),s2.denominator());
    for(int i = 0; i < 3; ++i)
      {
	printf("Lens space invariant s_2 for i=%d: \n",i+1);
	int ip1 = absolute_mod(i+1,3);
	int ip2 = absolute_mod(i+2,3);
	int params[4] = {k[ip1]-l[j], k[ip2]-l[j], k[ip1]-l[jp1], k[ip2]-l[jp1]};
        printf("    parameters: %d; %d, %d, %d, %d\n", k[i]-l[j], params[0], params[1], params[2], params[3]);
	s2 += -lens_s2(k[i]-l[j], params);
      }
    s2 = reduce_mod_ZZ(s2);
    printf("=> s2(E)  = %lld/%lld\n",s2.numerator(),s2.denominator());
  }
  void compute_s2_row(int j)
  {
    int jp1 = absolute_mod(j+1,2);
    printf("I'm using row %d",j+1);
    printf(" (%d,%d,%d)\n",k[j]-l[0],k[j]-l[1],k[j]-l[2]);
    //// printf("jp1 = %d",jp1+1);
    long long q = 
      square(k[j]-l[0])   + square(k[j]-l[1])   + square(k[j]-l[2]) +
      square(k[jp1]-l[0]) + square(k[jp1]-l[1]) + square(k[jp1]-l[2])
      - square(k[j]-k[jp1]);
    long long d = 16*3*(long long)(-mr)*(k[j]-l[0])*(k[j]-l[1])*(k[j]-l[2]);
    s2.assign(q-2,d);
    //// printf("q = %lld, d = %lld\n", q, d);
    //// printf("=> (q-2)/d = %lld/%lld\n",s2.numerator(),s2.denominator());
    for(int i = 0; i < 3; ++i)
      {
	//// printf("Lens space invariant s_2 for i=%d: \n",i+1);
	int ip1 = absolute_mod(i+1,3);
	int ip2 = absolute_mod(i+2,3);
	int params[4] = {k[j]-l[ip1], k[j]-l[ip2], k[jp1]-l[ip1], k[jp1]-l[ip2]};
        //// printf("     parameters: %d; %d, %d, %d, %d\n", k[i]-l[j], params[0], params[1], params[2], params[3]);
	s2 += lens_s2(k[j]-l[i], params);
      }
    s2 = reduce_mod_ZZ(s2);
    //// printf("=> s2(E)  = %lld/%lld\n",s2.numerator(),s2.denominator());
  }
  void compute_s2(void)
  {
    printf("\n===== E(%d, %d, %d | %d, %d, %d) ===== \n", k[0], k[1], k[2], l[0], l[1], l[2]);
    find_good_col_or_row();
    if (good_col >= 0)
	compute_s2_col(good_col);
    else if(good_row >= 0)
      compute_s2_row(good_row);
    else 
      s2.assign(666,1); // need better error handling here
  }
  void compute_s22(void) // only possible AFTER computing s2
  {
    //// printf("=> s2(E)  = %lld/%lld\n",s2.numerator(),s2.denominator());
    s22 = 2*absolute(mr)*s2;
    //// printf("=> s22(E) = %lld/%lld (in QQ)\n",s22.numerator(),s22.denominator());
    s22 = reduce_mod_ZZ(s22);
    //// printf("=> s22(E) = %lld/%lld (in QQ/ZZ)\n",s22.numerator(),s22.denominator());
  }
};


//////////////////////////////////////////////////
// Main routine:

main(){
  /* 
     printf("gcd(3,6) = %d\n",gcd(3,6));
     printf("gcd(3,3) = %d\n",gcd(3,3));
     printf("gcd(3,5) = %d\n",gcd(3,5));
     printf("gcd(-3,6) = %d\n",gcd(-3,6));
     printf("gcd(-3,3) = %d\n",gcd(-3,3));
     printf("gcd(-3,5) = %d\n",gcd(-3,5));
     printf("gcd(5,-5) = %d\n",gcd(5,-5));
     
     float a = M_PI;
     printf("%.6f\n",a);

     boost::rational<long long> myr(10,6);
     printf("%lld/%lld\n",myr.numerator(),myr.denominator());
     boost::rational<long long> myr2(0,6);
     printf("%lld/%lld\n",myr2.numerator(),myr2.denominator());

     boost::rational<long long> q(-49,3);
     printf("q = %lld/%lld\n",q.numerator(),q.denominator());
     q = reduce_mod_ZZ(q);
     printf("&q' = %lld/%lld\n",q.numerator(),q.denominator());
     printf("signed_mod(-49,3) = %d\n",signed_mod(-49,3));
  */

  /*
    "basic invariants":  r, s, p_1
    r, s, s_22      coincide  =>  spaces homotopy equivalent
    r, s, s_22, p_1 coincide  =>  spaces tangentially homotopy equivalent
    If values of s_2 differ, then spaces are not homeomorphic.
  */

  FILE *file_list_basic = fopen("list_basic.txt", "w");  // r, s, p1 agree
  FILE *file_list_the   = fopen("list_the.txt", "w");    // r, s, p1, s22 agree
  if (file_list_basic == NULL || file_list_the == NULL)
    {
      printf("Error opening file!\n");
      //exit(1);
    }

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

  int max_mr;
  printf("\n Maximum value of |r|: ");
  scanf("%d",&max_mr); 
  std::vector< struct Space> spaces [max_mr + 1];

  for (int l1 = 0; l1 < max_mr; ++l1){
    printf(" %d . ", l1);
    for (int l0 = l1; l0 < max_mr; ++l0){
      for (int k1 = l0+1; k1 <= max_mr; k1++){
	for (int k0 = k1; k0 <= max_mr; k0++){
	  struct Space new_space = {{k0, k1, l0 + l1 - k0 - k1}, {l0, l1, 0},};
	  new_space.compute_r_s_p1();
	  if (new_space.mr > max_mr) break; // note:  r is monotonous in k0 and k1
	  if (new_space.is_positively_curved_space())
	    spaces[new_space.mr].push_back(new_space);
	}}}}
  printf("\n");


  
  /*
  //////////////////////////////////////////////////
  // find homotopy equivalent pairs:
  
  //std::pair<struct Space, struct Space> spair = std::make_pair(E1,E2);
  //spair.first.print();

  // find pairs with matching r & |s|;
  //std::vector< std::pair<struct Space, struct Space> > basic_pairs; 
  for(int mr = 0; mr <= max_mr; mr++){ 
    for(int i1 = 0; i1 < spaces[mr].size(); ++i1){
      for(int i2 = i1+1; i2 < spaces[mr].size(); ++i2){
	struct Space E1 = spaces[mr][i1];
	struct Space E2 = spaces[mr][i2];
	if (E1.s == E2.s || E1.s == -E2.s) {
	  E1.compute_s2();
	  E1.compute_s22();
	  E2.compute_s2();
	  E2.compute_s22();
	  printf("\n");
	  boost::rational<long long> a_half(1,2);
	  if ((E1.s == E2.s && E1.s22 == E2.s22) // oriented case
	      || (E1.s == -E2.s && (E1.s22 == -E2.s22 || E1.s22 == a_half && E2.s22 == a_half))) // non-oriented case
	    {
	      printf("----------------------\n");
	      fprintf(file_list_basic,"----------------------\n");
	      E1.print(file_list_basic);
	      E2.print(file_list_basic);
	    }
	  ////std::pair<struct Space, struct Space> new_pair = std::make_pair(E1,E2);
	  ////spaces_pairs.push_back(new_pair);
	}
      }
    }
  }
  //////////////////////////////////////////////////
  */

  //////////////////////////////////////////////////
  // find pairs whose basic invariants agree
  int c_basic = 0; // counts pairs whose basic invariants agree
  int c_the   = 0; // counts tangentially homotopy equivalent pairs
  for(int mr = 0; mr <= max_mr; mr++){ 
    for(int i1 = 0; i1 < spaces[mr].size(); ++i1){
      for(int i2 = i1+1; i2 < spaces[mr].size(); ++i2){
	struct Space E1 = spaces[mr][i1];
	struct Space E2 = spaces[mr][i2];
	if ( (E1.s == E2.s || E1.s == -E2.s) && E1.p1 == E2.p1 ) 
	  {
	    ++c_basic;
	    printf("Found pair %d.\n", c_basic);
      
	    E1.compute_s2();
	    E1.compute_s22();
	    E2.compute_s2();
	    E2.compute_s22();

	    fprintf(file_list_basic,"%d: ----------------------\n", c_basic);
	    E1.print(file_list_basic);
	    E2.print(file_list_basic);

	    boost::rational<long long> a_half(1,2);
	    if ((E1.s == E2.s && E1.s22 == E2.s22) // oriented case
		|| (E1.s == -E2.s && (E1.s22 == -E2.s22 || E1.s22 == a_half && E2.s22 == a_half))) // non-oriented case
	      {
		++c_the;
		printf("--- This pair is tangentially homotopy equivalent. (%d)\n", c_the);
		fprintf(file_list_the,"%d: ----------------------\n", c_the);
		E1.print(file_list_the);
		E2.print(file_list_the);
	      }
	    ////std::pair<struct Space, struct Space> new_pair = std::make_pair(E1,E2);
	  ////spaces_pairs.push_back(new_pair);
	}
      }
    }
  }
  fclose(file_list_basic);
  fclose(file_list_the);
}