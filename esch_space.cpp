#include "esch_space.h"
using std::array;
#include "aux_math.h"

#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions/round.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>

//////////////////////////////////////////////////
// Lens space invariants:
rational<long long> lens_s2(long p, array<long,4> param)
{
  boost::multiprecision::cpp_dec_float_100 a = 0;
  for(long k = 1; k < abs(p); ++k)
    {
      boost::multiprecision::cpp_dec_float_100 s = 
	(cos(2*M_PI*k/abs(p)) - 1)
	/sin(k*M_PI*param[0]/p)
	/sin(k*M_PI*param[1]/p)
	/sin(k*M_PI*param[2]/p)
	/sin(k*M_PI*param[3]/p);
      a += s;
    }
  /*std::cout << "    ";
  std::cout << std::setprecision(std::numeric_limits<boost::multiprecision::cpp_dec_float_100>::max_digits10)
  	    << a*45 
  	    << std::endl; // print a*/
  long long rounded_45_a = (long long)boost::math::round(a*45);
  //printf("     round(...) = %lld\n",rounded_45_a);
  rational<long long> s2(rounded_45_a,45*16*p);
  //printf("     s2 = %lld/%lld\n",s2.numerator(),s2.denominator());
  return s2;
}

//////////////////////////////////////////////////
void Space::print_for_human(FILE* file) const 
{
  fprintf(file, "[%4ld,%4ld,%4ld |%4ld,%4ld,%4ld] -> r = %4ld, s = %4ld, p1 = %4ld, s22 = %lld/%lld, s2 = %lld/%lld", k_[0],k_[1],k_[2],l_[0],l_[1],l_[2],r_, s_, p1_, s22_.numerator(), s22_.denominator(), s2_.numerator(), s2_.denominator());
}

void Space::print_for_maple(FILE* file) const 
{
  fprintf(file, "[%ld,%ld,%ld,%ld,%ld,%ld]",k_[0],k_[1],k_[2],l_[0],l_[1],l_[2]);
}

void Space::setParameters(array<long,3> kkk, array<long,3> lll){ //compute_s_and_p1(void){
  k_ = kkk;
  l_ = lll;
  int sigma1_k = k_[0] + k_[1] + k_[2];
  int sigma2_k = k_[0]*k_[1] + k_[0]*k_[2] + k_[1]*k_[2];
  int sigma3_k = k_[0]*k_[1]*k_[2];
  int sigma1_l = l_[0] + l_[1] + l_[2];
  int sigma2_l = l_[0]*l_[1] + l_[0]*l_[2] + l_[1]*l_[2];
  int sigma3_l = l_[0]*l_[1]*l_[2];
  r_ = sigma2_k - sigma2_l;
  s_ = signed_mod(sigma3_k, abs(r_));  // note: sigma3_l = 0
  p1_ = absolute_mod(2*sigma1_k*sigma1_k - 6*sigma2_k, abs(r_));
  good_col_or_row = -1; // -1   means "not yet computed"
  s2_.assign(-1,1);     // -1/1 means "not yet computed"
  s22_.assign(-1,1);    // -1/1 means "not yet computed"
}

bool Space::is_positively_curved_space(void) {
    // Is it an Eschenburg space? 
    // -- Test conditions of [CEZ06] (1.1):
    if(boost::math::gcd(k_[0] - l_[0], k_[1] - l_[1]) > 1) return false;
    if(boost::math::gcd(k_[0] - l_[0], k_[1] - l_[2]) > 1) return false;
    if(boost::math::gcd(k_[0] - l_[2], k_[1] - l_[0]) > 1) return false;
    if(boost::math::gcd(k_[0] - l_[1], k_[1] - l_[0]) > 1) return false;
    if(boost::math::gcd(k_[0] - l_[1], k_[1] - l_[2]) > 1) return false;
    if(boost::math::gcd(k_[0] - l_[2], k_[1] - l_[1]) > 1) return false;
    // Is the space positively curved? 
    // -- Test conditions of [CEZ06] (1.2):
    int min_l = std::min(l_[0],std::min(l_[1],l_[2]));
    if (k_[0] == min_l) return false;
    if (k_[1] == min_l) return false;
    if (k_[2] == min_l) return false;
    int max_l = std::max(l_[0],std::max(l_[1],l_[2]));
    if (k_[0] == max_l) return false;
    if (k_[1] == max_l) return false;
    if (k_[2] == max_l) return false;
    return true;
  }

int Space::gcdA(int i, int j, int ii, int jj)
{ 
  // auxiliary function for find_good_row_or_col
  /*std::cout << k_[i] << " " << l_[j] << " " 
	    << k_[ii] << " " << l_[jj] << "--->" 
	    << boost::math::gcd(k_[i]-l_[j],k_[ii]-l_[jj]) << std::endl;*/
  return boost::math::gcd(k_[i]-l_[j],k_[ii]-l_[jj]); 
} 

bool Space::conditionC(void)
{
  // see "condition C" in [CEZ06]
  if (good_col_or_row == -1) // condition has not yet been checked
    {
      good_col_or_row = 6;  // if values stays 6, condition C is not satisfied
      for(int i = 0; i < 3; ++i)
	{
	  if(gcdA(i,0,i,1) == 1 && gcdA(i,0,i,2) == 1 && gcdA(i,1,i,2) == 1)
	    {
	      good_col_or_row = i + 3;// 3/4/5 = row 1/2/3 satisfies condition C
	      break;
	    }
	  else if(gcdA(0,i,1,i) == 1 && gcdA(0,i,2,i) == 1 && gcdA(1,i,2,i) == 1)
	    {
	      good_col_or_row = i;// 0/1/2 = col 1/2/3 satisfies condition C
	      break;
	    }
	}
    }
  return (good_col_or_row < 6);
}

void Space::compute_s2_col(int j)
{
  // [CEZ06] (2.1)
  // Again, this should really be a private function.
  // The "public" function is compute_s22 below.
  //
  //     s2  = (q-2)/d + SUM, 
  //     s22 = 2|r|s2
  //
  int jp1 = absolute_mod(j+1,2);
  //printf("I'm using column %d: ",j+1);
  //printf("(%ld,%ld,%ld)\n",k_[0]-l_[j],k_[1]-l_[j],k_[2]-l_[j]);
  //printf("jp1 = %d",jp1+1);
  long long q = 
    square(k_[0]-l_[j])   + square(k_[1]-l_[j])   + square(k_[2]-l_[j]) +
    square(k_[0]-l_[jp1]) + square(k_[1]-l_[jp1]) + square(k_[2]-l_[jp1])
    - square(l_[j]-l_[jp1]);
  long long d = 16*3*(long long)(r_)*(k_[0]-l_[j])*(k_[1]-l_[j])*(k_[2]-l_[j]);
  s2_.assign(q-2,d);
  //printf("q = %lld, r = %ld, d = %lld\n", q, (long)(r_), d);
  //printf("=> (q-2)/d = %lld/%lld\n",s2_.numerator(),s2_.denominator());
  for(int i = 0; i < 3; ++i)
    {
      // printf("Lens space invariant s_2 for i=%d: \n",i+1);
      int ip1 = absolute_mod(i+1,3);
      int ip2 = absolute_mod(i+2,3);
      array<long,4> params = {k_[ip1]-l_[j], k_[ip2]-l_[j], k_[ip1]-l_[jp1], k_[ip2]-l_[jp1]};
      // printf("    parameters: %ld; %ld, %ld, %ld, %ld\n", k_[i]-l_[j], params[0], params[1], params[2], params[3]);
      s2_ += -lens_s2(k_[i]-l_[j], params);
    }
  s2_ = reduce_mod_ZZ(s2_);
  //printf("=> s2(E)  = %lld/%lld\n",s2_.numerator(),s2_.denominator());
}

void Space::compute_s2_row(int j)
  {
    int jp1 = absolute_mod(j+1,2);
    //printf("I'm using row %d",j+1);
    //printf(" (%ld,%ld,%ld)\n",k_[j]-l_[0],k_[j]-l_[1],k_[j]-l_[2]);
    //printf("jp1 = %d",jp1+1);
    long long q = 
      square(k_[j]-l_[0])   + square(k_[j]-l_[1])   + square(k_[j]-l_[2]) +
      square(k_[jp1]-l_[0]) + square(k_[jp1]-l_[1]) + square(k_[jp1]-l_[2])
      - square(k_[j]-k_[jp1]);
    long long d = 16*3*(long long)(r_)*(k_[j]-l_[0])*(k_[j]-l_[1])*(k_[j]-l_[2]);
    s2_.assign(q-2,d);
    //printf("q = %lld, d = %lld\n", q, d);
    // printf("=> (q-2)/d = %lld/%lld\n",s2_.numerator(),s2_.denominator());
    for(int i = 0; i < 3; ++i)
      {
	//   printf("Lens space invariant s_2 for i=%d: \n",i+1);
	int ip1 = absolute_mod(i+1,3);
	int ip2 = absolute_mod(i+2,3);
	array<long,4> params = {k_[j]-l_[ip1], k_[j]-l_[ip2], k_[jp1]-l_[ip1], k_[jp1]-l_[ip2]};
        //   printf("     parameters: %ld; %ld, %ld, %ld, %ld\n", k_[j]-l_[i], params[0], params[1], params[2], params[3]);
	s2_ += lens_s2(k_[j]-l_[i], params);
      }
    s2_ = reduce_mod_ZZ(s2_);
    // printf("=> s2(E)  = %lld/%lld\n",s2_.numerator(),s2_.denominator());
  }

const rational<long long>& Space::s2(void)
{
  if (s2_ == rational<long long>(-1)) // s2 not yet computed
    {
      //printf("\n===== E(%ld, %ld, %ld | %ld, %ld, %ld) ===== \n", k_[0], k_[1], k_[2], l_[0], l_[1], l_[2]);
      if (conditionC() == true)
	{
	  if (good_col_or_row < 3)
	    compute_s2_col(good_col_or_row);
	  else
	    compute_s2_row(good_col_or_row-3);
	}
      else 
	s2_.assign(1,1); // need better error handling here
    }
  return s2_;
}
const rational<long long>& Space::s22(void)
{
  if (s22_ == rational<long long>(-1)) // s22 not yet computed
    {
      //// printf("=> s2(E)  = %ld/%ld\n",s2.numerator(),s2.denominator());
      s22_ = 2*abs(r_)*s2();  // use s2() instead of s2_ instead s2_ has not yet been computed
      //// printf("=> s22(E) = %ld/%ld (in QQ)\n",s22_.numerator(),s22_.denominator());
      s22_ = reduce_mod_ZZ(s22_);
      //// printf("=> s22(E) = %ld/%ld (in QQ/ZZ)\n",s22_.numerator(),s22_.denominator());
    }
  return s22_;
}

//////////////////////////////////////////////////
bool compare_spaces(class Space E1, class Space E2)
{
  return (E1.p1() < E2.p1())  
    || (E1.p1() == E2.p1() && abs(E1.s()) < abs(E2.s()));
}


