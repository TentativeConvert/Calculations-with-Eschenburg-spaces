/* compile with: 
   g++ -std=c++11 eschenburg.c
*/
#include<cstdio>
#include<vector>

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
int absolute_mod (int a, int base)
{
  int remainder = a % base;
  if(remainder < 0)
    remainder += base;
  return remainder;
}
int maximum (int a, int b)
{
  a > b ? a : b;
}
int minimum (int a, int b)
{
  a < b ? a : b;
}
int absolute (int a)
{
  a >= 0 ? a : -a;
}
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
//  Spaces:
struct Space {
  int k1;
  int k2;
  int l1;
  int l2;
  int mr; // mr = -r = "minus r"
  int s;
  int p1;
  void print(FILE *file){
    printf(      "[%4d,%4d,%4d |%4d,%4d, 0] --> r = %4d, s = %4d, p1 = %4d\n",k1,k2,l1+l2-k1-k2,l1,l2,-mr,s,p1);
    fprintf(file,"[%4d,%4d,%4d |%4d,%4d, 0] --> r = %4d, s = %4d, p1 = %4d\n",k1,k2,l1+l2-k1-k2,l1,l2,-mr,s,p1);
  }
  void compute_r_s_p1(void){
     int k3 = l1 + l2 - k1 - k2;
     int sigma1_k = l1 + l2;
     int sigma2_k = k1*k2 + (k1 + k2)*k3;
     int sigma3_k = k1*k2*k3;
     int sigma2_l = l1*l2;
     mr = sigma2_l - sigma2_k;
     s = signed_mod(sigma3_k, mr);  // note: sigma3_l = 0
     p1 = absolute_mod(2*sigma1_k*sigma1_k - 6*sigma2_k, mr);
  }
  bool is_positively_curved_space(void) {
    // Is it an Eschenburg space? 
    // -- Test conditions of [CEZ06] (1.1):
    int k3 = l1 + l2 - k1 - k2;
    int l3 = 0;
    // printf("\n[%4d,%4d,%4d |%4d,%4d,%4d]",k1,k2,k3,l1,l2,l3);
    if(gcd(k1 - l1, k2 - l2) > 1) return false;
    if(gcd(k1 - l1, k2 - l3) > 1) return false;
    if(gcd(k1 - l3, k2 - l1) > 1) return false;
    if(gcd(k1 - l2, k2 - l1) > 1) return false;
    if(gcd(k1 - l2, k2 - l3) > 1) return false;
    if(gcd(k1 - l3, k2 - l2) > 1) return false;
    // Is the space positively curved? 
    // -- Test conditions of [CEZ06] (1.2):
    int min = minimum(l1,l2);
    // printf("min: %d  ", min);
    if (k1 == min) return false;
    if (k2 == min) return false;
    if (k3 == min) return false;
    int max = maximum(l1,l2);
    // printf("max: %d  ", max);
    if (k1 == max) return false;
    if (k2 == max) return false;
    if (k3 == max) return false;
    printf("positively curved");
    return true;
  }
};


//////////////////////////////////////////////////
// Main routine:

main(){
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
 
  int max_mr;
  printf("\n Maximum value of |r|: ");
  scanf("%d",&max_mr); 
  std::vector< struct Space> spaces [max_mr + 1];

  int k1, k2, l1, l2;
  for (l2 = 0; l2 < max_mr; l2++){
    printf(" %d . ", l2);
    for (l1 = l2; l1 < max_mr; l1++){
      for (k2 = l1+1; k2 <= max_mr; k2++){
	for (k1 = k2; k1 <= max_mr; k1++){
	  struct Space new_space{.k1 = k1, .k2 = k2, .l1 = l1, .l2 = l2};
	  new_space.compute_r_s_p1();
	  //new_space = compute_r_s_p1(new_space);
	  if (new_space.mr > max_mr) break; // note:  r is monotonous in k1 and k2
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
    for(int s1 = 0; s1 < spaces[mr].size(); s1++){
      for(int s2 = s1+1; s2 < spaces[mr].size(); s2++){
	struct Space E1 = spaces[mr][s1];
	struct Space E2 = spaces[mr][s2];
	if (E1.s == E2.s || E1.s == -E2.s) {
	  printf("----------------------\n");
	  fprintf(output_file,"----------------------\n");
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
