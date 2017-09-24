/* compile with: 
   g++ -std=c++11 eschenburlg.c
*/
#include<cstdio>
#include<vector>

int mod (int a, int b){ // see https://stackoverflow.com/a/4003287/3611932
   int ret = a % b;
   if(ret < 0)
     ret+=b;
   return ret;
}

struct Space {
  int k1;
  int k2;
  int l1;
  int l2;
  int mr; // mr = -r = "minus r"
  int s;
  int p1;
  void print(void){
    printf("[%4d,%4d,%4d |%4d,%4d, 0] --> r = %4d, s = %4d, p1 = %4d\n",k1,k2,l1+l2-k1-k2,l1,l2,-mr,s,p1);
  }
};

struct Space compute_r_s_p1(struct Space E){
  int k3 = E.l1 + E.l2 - E.k1 - E.k2;
  int sigma1_k = E.l1 + E.l2;
  int sigma2_k = E.k1*E.k2 + (E.k1 + E.k2)*k3;
  int sigma3_k = E.k1*E.k2*k3;
  int sigma2_l = E.l1*E.l2;
  E.mr = sigma2_l - sigma2_k;
  E.s = mod(sigma3_k, E.mr);  // note: sigma3_l = 0
  E.p1 = mod(2*sigma1_k*sigma1_k - 6*sigma2_k, E.mr);
  return E;
}

main(){
  // by Lemma 1.4 we have "standard parametrization" by 
  //    k1 >= k2 > l1 >= l2 >= 0;
  // by proof of Prop. 1.7,  
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
	  new_space = compute_r_s_p1(new_space);
	  if (new_space.mr > max_mr) break; //note:  r is monotonous in k1 and k2
	  spaces[new_space.mr].push_back(new_space);
	}}}}
  printf("\n");

  //print all spaces, sorted according to value of r:
  /*
    for(int mr = 0; mr <= max_mr; mr++){ 
      for(int s = 0; s < spaces[mr].size(); s++){
	spaces[mr][s].print();
      }
      printf("----------------------\n");
    }
  */

  /*struct Space E1 {.k1 = 3, .k2 = 3, .l1 = 2, .l2 = 2};
  E1 = compute_r_s_p1(E1);
  E1.print();
  struct Space E2 {.k1 = 5, .k2 = 4, .l1 = 3, .l2 = 2};
  E2 = compute_r_s_p1(E2);
  E2.print();
  std::pair<struct Space, struct Space> spair = std::make_pair(E1,E2);
  spair.first.print();*/

  // find pairs with matching r & s;
  ////  std::vector< std::pair<struct Space, struct Space> > space_pairs; 
  for(int mr = 0; mr <= max_mr; mr++){ 
    for(int s1 = 0; s1 < spaces[mr].size(); s1++){
      for(int s2 = s1+1; s2 < spaces[mr].size(); s2++){
	struct Space E1 = spaces[mr][s1];
	struct Space E2 = spaces[mr][s2];
	if (E1.s == E2.s) {
	  printf("----------------------\n");
	  E1.print();
	  E2.print();
	  ////std::pair<struct Space, struct Space> new_pair = std::make_pair(E1,E2);
	  ////spaces_pairs.push_back(new_pair);
	}
      }
    }
  }

}
