#include "esch_space.h"
using std::array;
const boost::rational<long long> Space::KS_UNKNOWN = boost::rational<long long>(-1,1);
const boost::rational<long long> Space::KS_UNCOMPUTABLE = boost::rational<long long>(1,1); 

#include "aux_math.h"
#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions/round.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>

//////////////////////////////////////////////////

void Space::print(FILE* file) const 
{
  fprintf(file, " [%4ld,%4ld,%4ld, %4ld,%4ld,%4ld] -> r = %5ld, s = %5ld, p1 = %5ld",k_[0],k_[1],k_[2],l_[0],l_[1],l_[2],r_, s_, p1_);
  if (s2_ == KS_UNKNOWN)
    fprintf(file, "\n");
  else if (s2_ == KS_UNCOMPUTABLE)
    fprintf(file, "  |!| WARNING: Condition C not satisfied |!|\n");
  else
    fprintf(file, ", s22 = %lld/%lld, s2 = %lld/%lld\n", s22_.numerator(), s22_.denominator(), s2_.numerator(), s2_.denominator());
}

void Space::print(void) const 
{
  if (is_positively_curved_space()) 
    {
      printf("\nInvariants of the positively curved Eschenburg space with parameters [%ld,%ld,%ld, %ld,%ld,%ld]:\n", k_[0],k_[1],k_[2],l_[0],l_[1],l_[2]);
      printf("  r = %ld, s = %ld, p1 = %ld\n",r_, s_, p1_);      
      if (s2_ == KS_UNCOMPUTABLE)
	printf("! Condition C is not satisfied !\n\n");
      else
	printf("  s22 = %lld/%lld, s2 = %lld/%lld\n\n", s22_.numerator(), s22_.denominator(), s2_.numerator(), s2_.denominator());
    }
  else
    printf("\nThe supplied parameters do not describe a positively curved Eschenburg space.\n\n");
}	

//////////////////////////////////////////////////

Space::Space(array<long,3> kkk, array<long,3> lll){
  k_ = kkk;
  l_ = lll;
  int sigma1_k = k_[0] + k_[1] + k_[2];
  int sigma2_k = k_[0]*k_[1] + k_[0]*k_[2] + k_[1]*k_[2];
  int sigma3_k = k_[0]*k_[1]*k_[2];
  int sigma1_l = l_[0] + l_[1] + l_[2];
  int sigma2_l = l_[0]*l_[1] + l_[0]*l_[2] + l_[1]*l_[2];
  int sigma3_l = l_[0]*l_[1]*l_[2];
  r_ = sigma2_k - sigma2_l;
  s_ = signed_mod(sigma3_k-sigma3_l, abs(r_))*sign(r_);
  p1_ = absolute_mod(2*sigma1_k*sigma1_k - 6*sigma2_k, abs(r_));
  good_col_or_row = GOOD_CoR_UNKNOWN; 
  s2_ =  KS_UNKNOWN;
  s22_ = KS_UNKNOWN;
}

bool Space::is_positively_curved_space(void) const {
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

//////////////////////////////////////////////////
// KS-invariants

bool Space::compute_KS_invariants() // compute s2 & s22
{
  if (test_condition_C())
    {
      if (good_col_or_row < 3)
	compute_s2_from_col(good_col_or_row);
      else
	compute_s2_from_row(good_col_or_row-3);
      s22_ = 2*abs(r_)*s2_;  
      s22_ = reduce_mod_ZZ(s22_);
      return true;
    }
  else 
    {
      s2_ = KS_UNCOMPUTABLE; 
      s22_ = KS_UNCOMPUTABLE;
      return false;
    }
} 

bool Space::test_condition_C(void) // check condition C & find good column or row
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

int Space::gcdA(int i, int j, int ii, int jj) // helper function for testing condition C
{ 
  // auxiliary function for find_good_row_or_col
  /*std::cout << k_[i] << " " << l_[j] << " " 
	    << k_[ii] << " " << l_[jj] << "--->" 
	    << boost::math::gcd(k_[i]-l_[j],k_[ii]-l_[jj]) << std::endl;*/
  return boost::math::gcd(k_[i]-l_[j],k_[ii]-l_[jj]); 
} 

void Space::compute_s2_from_col(int j)
{
  // see [CEZ06] (2.1):
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
      int ip1 = absolute_mod(i+1,3);
      int ip2 = absolute_mod(i+2,3);
      array<long,4> params = {k_[ip1]-l_[j], k_[ip2]-l_[j], k_[ip1]-l_[jp1], k_[ip2]-l_[jp1]};
      s2_ += -lens_s2(k_[i]-l_[j], params);
    }
  s2_ = reduce_mod_ZZ(s2_);
  //printf("=> s2(E)  = %lld/%lld\n",s2_.numerator(),s2_.denominator());
}

void Space::compute_s2_from_row(int j)
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


//////////////////////////////////////////////////
// Lens space invariants:
rational<long long> Space::lens_s2(long p, array<long,4> param)
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

Space::comp Space::compareBasicType(const Space& E1, const Space& E2)
{
  if (abs(E1.p1_) > abs(E2.p1_)) return comp::GREATER;
  if (abs(E1.p1_) < abs(E2.p1_)) return comp::SMALLER;
  return comp::EQUAL;
}

Space::comp Space::compareHomotopyType(const Space& E1, const Space& E2)
{ 
  if (abs(E1.s22_) > abs(E2.s22_)) return comp::GREATER;
  if (abs(E1.s22_) < abs(E2.s22_)) return comp::SMALLER;
  if (sign(E1.s22_)*sign(E1.s_) > sign(E2.s22_)*sign(E2.s_)) return comp::GREATER;
  if (sign(E1.s22_)*sign(E1.s_) < sign(E2.s22_)*sign(E2.s_)) return comp::SMALLER;
  if (abs(E1.s_) > abs(E2.s_)) return comp::GREATER;
  if (abs(E1.s_) < abs(E2.s_)) return comp::SMALLER;
  if (abs(E1.r_) > abs(E2.r_)) return comp::GREATER;
  if (abs(E1.r_) < abs(E2.r_)) return comp::SMALLER;
  return comp::EQUAL;
}

Space::comp Space::compareTangentialHomotopyType(const Space& E1, const Space& E2)
{
  if (E1.p1_ > E2.p1_) return comp::GREATER;
  if (E1.p1_ < E2.p1_) return comp::SMALLER;
  return compareHomotopyType(E1,E2);
}

Space::comp Space::compareHomeomorphismType(const Space& E1, const Space& E2)
{
  if (E1.s2_ > E2.s2_) return comp::GREATER;
  if (E1.s2_ < E2.s2_) return comp::SMALLER;
  if (sign(E1.s2_)*sign(E1.s_) > sign(E2.s2_)*sign(E2.s_)) return comp::GREATER;
  if (sign(E1.s2_)*sign(E1.s_) < sign(E2.s2_)*sign(E2.s_)) return comp::SMALLER;
  return compareHomotopyType(E1,E2);
}
