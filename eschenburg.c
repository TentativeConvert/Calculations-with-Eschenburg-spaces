/* compile with: 
   g++ -std=c++11 eschenburg.c
*/
#include<cstdio>
#include<vector>
#include<cmath>
#include<boost/rational.hpp>

//////////////////////////////////////////////////
// Auxiliary mathematics:

int signed_mod (int a, int base)
// In general, sign of % is machine dependent, 
// see https://stackoverflow.com/a/4003287/3611932.
// For the computation of s mod r,
// I want modulus numbers to be "symmetric around zero",
// e.g. s in {2,1,0,-1,-2} for r = 5.
{
  int remainder = a % base;
  if(remainder > base/2)
    remainder -= base;
  else if(remainder < -base/2)
    remainder += base;
  return remainder;
}
int absolute (int a){ return a >= 0 ? a : -a; }
int absolute_mod (int a, int base)
{
  int remainder = a % base;
// if (a < 0)
    // printf("\n%d mod %d = %d, |base| = %d \n", a, base, remainder, absolute(base));
  if(remainder < 0)
    remainder += absolute(base);
// if (a < 0) 
    // printf("%d mod %d = %d\n", a, base, remainder);
  return remainder;
}
int maximum (int a, int b){ return a > b ? a : b; }
int minimum (int a, int b){ return a < b ? a : b; }
int gcd (int a, int b)
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
}

//////////////////////////////////////////////////
// Lens space invariants

boost::rational<int> lens_s2(int p, int param[4])
{
  float a = 0;
  for(int k = 1; k < absolute(p); ++k)
    {
      float s = 
	(cos(2*M_PI*k/absolute(p)) - 1)
	/sin(k*M_PI*param[0]/p)
	/sin(k*M_PI*param[1]/p)
	/sin(k*M_PI*param[2]/p)
	/sin(k*M_PI*param[3]/p);
      a += s;
    }
  //// printf("16*p*45*s2 = %.6f -- ",a*45);// this is supposed to be an integer
  //// printf("round(a*45) = %d -- ",(int)round(a*45));
  //// printf("45*16*p = %d --", 45*16*p);
  boost::rational<int> s2((int)round(a*45),45*16*p);
  //// printf("s2 = %d/%d\n",s2.numerator(),s2.denominator());
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
  boost::rational<int> s22;
  int good_row; // for "condition C"
  int good_col; // for "condition C"
  void print(FILE *file){
    printf(      "[%4d,%4d,%4d |%4d,%4d,%4d] --> r = %4d, s = %4d, p1 = %4d, s22 = %d/%d\n",k[0],k[1],k[2],l[0],l[1],l[2],-mr,s,p1,s22.numerator(),s22.denominator());
    fprintf(file,"[%4d,%4d,%4d |%4d,%4d,%4d] --> r = %4d, s = %4d, p1 = %4d, s22 = %d/%d\n",k[0],k[1],k[2],l[0],l[1],l[2],-mr,s,p1,s22.numerator(),s22.denominator());
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
    if(gcd(k[0] - l[0], k[1] - l[1]) > 1) return false;
    if(gcd(k[0] - l[0], k[1] - l[2]) > 1) return false;
    if(gcd(k[0] - l[2], k[1] - l[0]) > 1) return false;
    if(gcd(k[0] - l[1], k[1] - l[0]) > 1) return false;
    if(gcd(k[0] - l[1], k[1] - l[2]) > 1) return false;
    if(gcd(k[0] - l[2], k[1] - l[1]) > 1) return false;
    // Is the space positively curved? 
    // -- Test conditions of [CEZ06] (1.2):
    int min_l = minimum(l[0],minimum(l[1],l[2]));
    if (k[0] == min_l) return false;
    if (k[1] == min_l) return false;
    if (k[2] == min_l) return false;
    int max_l = maximum(l[0],maximum(l[1],l[2]));
    if (k[0] == max_l) return false;
    if (k[1] == max_l) return false;
    if (k[2] == max_l) return false;
    return true;
  }
  int gcdA(int i,int j,int ii,int jj){ 
    // auxiliary function for find_good_row_or_col
    // should be declared private, so I probably ought to make Space a proper class, not just a struct
    /*int ret = gcd(k[i]-l[j],k[ii]-l[jj]); 
    if (ret==1) 
      {
	printf("\n k[%d]-l[%d] = %d; ", i,j,k[i]-l[j]);
	printf("k[%d]-l[%d] = %d; ", ii,jj,k[ii]-l[jj]);
	}*/
    return gcd(k[i]-l[j],k[ii]-l[jj]); 
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
	  //// printf("(%d,%d,%d | %d,%d,%d) -- ",k[0],k[1],k[2],l[0],l[1],l[2]);
	  //// printf("using column %d: (%d,%d,%d)\n",c+1,k[0]-l[c],k[1]-l[c],k[2]-l[c]);
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
  void compute_s22_col(int j)
  {
    // [CEZ06] (2.1)
    // Again, this should really be a private function.
    // The "public" function is compute_s22 below.
    //
    //     s2  = (q-2)/d + SUM, 
    //     s22 = 2|r|s2
    //
    int jp1 = absolute_mod(j+1,2);
    int q = 
      (k[0]-l[j])^2   + (k[1]-l[j])^2   + (k[2]-l[j])^2 +
      (k[0]-l[jp1])^2 + (k[1]-l[jp1])^2 + (k[2]-l[jp1])^2 
      - (l[j]-l[jp1])^2; 
    int d = 16*3*(-mr)*(k[0]-l[j])*(k[1]-l[j])*(k[2]-l[j]);
    boost::rational<int> s2(q-2,d);
    //// printf("1: %d/%d\n",s2.numerator(),s2.denominator());
    for(int i = 0; i < 3; ++i)
      {
	int ip1 = absolute_mod(i+1,3);
	int ip2 = absolute_mod(i+2,3);
	int params[4] = {k[ip1]-l[j], k[ip2]-l[j], k[ip1]-l[jp1], k[ip2]-l[jp1]};
	s2 += -lens_s2(k[i]-l[j], params);
      }
    //// printf("2: %d/%d\n",s2.numerator(),s2.denominator());
    s22 = 2*absolute(mr)*s2;
    //// printf("3: %d/%d\n",s22.numerator(),s22.denominator());
    // s2 is now computed as an element of QQ.
    // The final value should be in QQ/ZZ, so we forget integral part:
    int s22_n = s22.numerator();
    int s22_d = s22.denominator();
    if(s22_d < 0) // probaby unnecessary, but depends on implementation of boost::rational
      {
      s22_n = -s22_n;
      s22_d = -s22_d;
      }
    //// printf("4: %d/%d\n",s22_n,s22_d);
    s22_n = absolute_mod(s22_n,s22_d);  
    //// printf("5: %d/%d\n",s22_n,s22_d);
    s22.assign(s22_n,s22_d); // this should simplify the fraction
    //// printf("6: %d/%d\n",s22.numerator(),s22.denominator());
  }
  void compute_s22_row(int j)
  {
    // still need to write this
    s22.assign(33,1);
  }
  void compute_s22(void)
  {
    find_good_col_or_row();
    if (good_col >= 0)
	compute_s22_col(good_col);
    else if(good_row >= 0)
      compute_s22_row(good_row);
    else 
      s22.assign(666,1); // need better error handling here
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

     boost::rational<int> myr(10,6);
     printf("%d/%d\n",myr.numerator(),myr.denominator());
     boost::rational<int> myr2(0,6);
     printf("%d/%d\n",myr2.numerator(),myr2.denominator());
  */

  FILE *output_file = fopen("output.txt", "w");
  if (output_file == NULL)
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

  /*struct Space E1 {.k1 = 3, .k2 = 3, .l1 = 2, .l2 = 2};
  E1 = compute_r_s_p1(E1);
  E1.print();
  struct Space E2 {.k1 = 5, .k2 = 4, .l1 = 3, .l2 = 2};
  E2 = compute_r_s_p1(E2);
  E2.print();
  std::pair<struct Space, struct Space> spair = std::make_pair(E1,E2);
  spair.first.print();*/

  // find pairs with matching r & |s|;
  ////  std::vector< std::pair<struct Space, struct Space> > space_pairs; 
  for(int mr = 0; mr <= max_mr; mr++){ 
    for(int i1 = 0; i1 < spaces[mr].size(); ++i1){
      for(int i2 = i1+1; i2 < spaces[mr].size(); ++i2){
	struct Space E1 = spaces[mr][i1];
	struct Space E2 = spaces[mr][i2];
	if (E1.s == E2.s || E1.s == -E2.s) {
	  printf("----------------------\n");
	  fprintf(output_file,"----------------------\n");
	  E1.compute_s22();
	  E2.compute_s22();
	  E1.print(output_file);
	  E2.print(output_file);
	  ////std::pair<struct Space, struct Space> new_pair = std::make_pair(E1,E2);
	  ////spaces_pairs.push_back(new_pair);
	}
      }
    }
  }
  fclose(output_file);
}
