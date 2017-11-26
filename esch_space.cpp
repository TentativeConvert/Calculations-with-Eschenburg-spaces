#include "esch_space.h"
using std::array;
const boost::rational<INT_KS> Space::KS_UNKNOWN = boost::rational<INT_KS>(-1,1);
const boost::rational<INT_KS> Space::KS_UNCOMPUTABLE = boost::rational<INT_KS>(1,1); 

#include "aux_math.h"
#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions/round.hpp>
#include <boost/math/special_functions/sin_pi.hpp>
#include <boost/math/special_functions/cos_pi.hpp>

//////////////////////////////////////////////////

void Space::print(FILE* file) const 
{
  fprintf(file, " [%4ld,%4ld,%4ld, %4ld,%4ld,%4ld] -> |r| = %5ld,  s = %7ld,  M1 = %2d,  M2 = %d,  p1 = %5ld",
	  (long)k_[0],(long)k_[1],(long)k_[2],(long)l_[0],(long)l_[1],(long)l_[2],(long)(abs(r_)), (long)s_, (int)(M1_), (int)(M2_), (long)p1_);
  if (s2_ == KS_UNKNOWN)
    fprintf(file, "\n");
  else if (s2_ == KS_UNCOMPUTABLE)
    fprintf(file, "  |!| WARNING: condition C not satisfied |!|\n");
  else
    {
      boost::rational<INT_KS> s22 = this->s22();
      fprintf(file, ",  s22 = %3lld/%lld,  s2 = %7lld/%lld\n", (long long)s22.numerator(), (long long)s22.denominator(), (long long)s2_.numerator(), (long long)s2_.denominator());
    }
}

void Space::print(void) const 
{
  if (is_space()) 
    {
      printf("\nInvariants of the Eschenburg space with parameters [%ld,%ld,%ld, %ld,%ld,%ld]:\n", 
	     (long)k_[0],(long)k_[1],(long)k_[2],(long)l_[0],(long)l_[1],(long)l_[2]);
      printf("  |r| = %ld,  s = %ld,  M1 = %d,  M2 = %d\n  p1  = %ld\n",(long)(abs(r_)), (long)s_, (int)(M1_), (int)(M2_), (long)p1_);
      if (s2_ == KS_UNCOMPUTABLE)
	printf("! Condition C is not satisfied !\n");
      else
	{
	  boost::rational<INT_KS> s22 = this->s22();
	  printf("  s2  = %lld/%lld,  s22 = %lld/%lld\n", (long long)s2_.numerator(), (long long)s2_.denominator(), (long long)s22.numerator(), (long long)s22.denominator());
	}
      if (is_positively_curved())
	printf("  The space is positively curved.\n\n");
      else
	printf("  The space is NOT positively curved.\n\n");
    }
  else
    printf("\nThe supplied parameters do not describe an Eschenburg space.\n\n");
}	

//////////////////////////////////////////////////

Space::Space(array<INT_P,3> kkk, array<INT_P,3> lll){
  k_ = kkk;
  l_ = lll;
  INT_R sigma1_k = (INT_R)k_[0] + (INT_R)k_[1] + (INT_R)k_[2];
  INT_R sigma2_k = (INT_R)k_[0]*(INT_R)k_[1] + (INT_R)k_[0]*(INT_R)k_[2] + (INT_R)k_[1]*(INT_R)k_[2];
  INT_R sigma3_k = (INT_R)k_[0]*(INT_R)k_[1]*(INT_R)k_[2];
  INT_R sigma1_l = (INT_R)l_[0] + (INT_R)l_[1] + (INT_R)l_[2];
  INT_R sigma2_l = (INT_R)l_[0]*(INT_R)l_[1] + (INT_R)l_[0]*(INT_R)l_[2] + (INT_R)l_[1]*(INT_R)l_[2];
  INT_R sigma3_l = (INT_R)l_[0]*(INT_R)l_[1]*(INT_R)l_[2];
  r_ = sigma2_k - sigma2_l;
  s_ = signed_mod(sigma3_k-sigma3_l, abs(r_))*sign(r_);
  M1_ = signed_mod(-sigma1_l,3);
  M2_ = absolute_mod(sigma1_l + sigma2_l,2);
  p1_ = absolute_mod(2*sigma1_k*sigma1_k - 6*sigma2_k, abs(r_));
  good_col_or_row = GOOD_CoR_UNKNOWN; 
  s2_ =  KS_UNKNOWN;
}

bool Space::is_space(void) const {
  // Is it an Eschenburg space? 
  // -- Test conditions of [CEZ07] (1.1):
  using boost::math::gcd;
  if(k_[0]+k_[1]+k_[2] != l_[0]+l_[1]+l_[2]) return false;
  if(gcd(k_[0] - l_[0], k_[1] - l_[1]) > 1) return false;
  if(gcd(k_[0] - l_[0], k_[1] - l_[2]) > 1) return false;
  if(gcd(k_[0] - l_[2], k_[1] - l_[0]) > 1) return false;
  if(gcd(k_[0] - l_[1], k_[1] - l_[0]) > 1) return false;
  if(gcd(k_[0] - l_[1], k_[1] - l_[2]) > 1) return false;
  if(gcd(k_[0] - l_[2], k_[1] - l_[1]) > 1) return false;
  return true;
}

bool Space::is_positively_curved(void) const {
  // Is the space positively curved? 
  // -- Test conditions of [CEZ07] (1.2):
  INT_P min_l = std::min(l_[0],std::min(l_[1],l_[2]));
  INT_P max_l = std::max(l_[0],std::max(l_[1],l_[2]));
  if (min_l <= k_[0] && k_[0] <= max_l) return false;
  if (min_l <= k_[1] && k_[1] <= max_l) return false;
  if (min_l <= k_[2] && k_[2] <= max_l) return false;
  return true;
}

//////////////////////////////////////////////////
// KS-invariants

boost::rational<INT_KS> Space::s22() const
{
  if (s2_ == KS_UNKNOWN || s2_ == KS_UNCOMPUTABLE)
    return s2_;
  else
    return reduce_mod_ZZ( 2*abs(r_)*s2_ );
}

bool Space::compute_KS_invariants() // compute s2 & s22
{
  if (test_condition_C())
    {
      if (good_col_or_row < 3)
	compute_s2_from_col(good_col_or_row);
      else
	compute_s2_from_row(good_col_or_row-3);
      return true;
    }
  else 
    {
      s2_ = KS_UNCOMPUTABLE; 
      return false;
    }
} 

bool Space::test_condition_C(void) // check condition C & find good column or row
{
  // see "condition C" in [CEZ07]
  if (good_col_or_row == GOOD_CoR_UNKNOWN) // condition has not yet been checked
    {
      good_col_or_row = GOOD_CoR_NONEXISTENT; // if value stays, condition C is not satisfied
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
  return (good_col_or_row != GOOD_CoR_NONEXISTENT);
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
  // see [CEZ07] (2.1):
  //
  //     s2  = (q-2)/d + SUM, 
  //     s22 = 2|r|s2
  //
  int jp1 = absolute_mod(j+1,2);
  //printf("I'm using column %d: ",j+1);
  //printf("(%ld,%ld,%ld)\n",k_[0]-l_[j],k_[1]-l_[j],k_[2]-l_[j]);
  //printf("jp1 = %d",jp1+1);
  INT_KS q = 
    square(k_[0]-l_[j])   + square(k_[1]-l_[j])   + square(k_[2]-l_[j]) +
    square(k_[0]-l_[jp1]) + square(k_[1]-l_[jp1]) + square(k_[2]-l_[jp1])
    - square(l_[j]-l_[jp1]);
  INT_KS d = 16*3*(INT_KS)(r_)*(k_[0]-l_[j])*(k_[1]-l_[j])*(k_[2]-l_[j]);
  s2_.assign(q-2,d);
  //printf("q = %lld, r = %lld, d = %lld\n", q, (long long)(r_), d);
  //printf("=> (q-2)/d = %lld/%lld\n",s2_.numerator(),s2_.denominator());
  for(int i = 0; i < 3; ++i)
    {
      int ip1 = absolute_mod(i+1,3);
      int ip2 = absolute_mod(i+2,3);
      array<INT_P,4> params = {k_[ip1]-l_[j], k_[ip2]-l_[j], k_[ip1]-l_[jp1], k_[ip2]-l_[jp1]};
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
  INT_KS q = 
    square(k_[j]-l_[0])   + square(k_[j]-l_[1])   + square(k_[j]-l_[2]) +
    square(k_[jp1]-l_[0]) + square(k_[jp1]-l_[1]) + square(k_[jp1]-l_[2])
    - square(k_[j]-k_[jp1]);
  INT_KS d = 16*3*(INT_KS)(r_)*(k_[j]-l_[0])*(k_[j]-l_[1])*(k_[j]-l_[2]);
  s2_.assign(q-2,d);
  //printf("q = %lld, d = %lld\n", q, d);
  // printf("=> (q-2)/d = %lld/%lld\n",s2_.numerator(),s2_.denominator());
  for(int i = 0; i < 3; ++i)
    {
      //   printf("Lens space invariant s_2 for i=%d: \n",i+1);
      int ip1 = absolute_mod(i+1,3);
      int ip2 = absolute_mod(i+2,3);
      array<INT_P,4> params = {k_[j]-l_[ip1], k_[j]-l_[ip2], k_[jp1]-l_[ip1], k_[jp1]-l_[ip2]};
      //   printf("     parameters: %ld; %ld, %ld, %ld, %ld\n", k_[j]-l_[i], params[0], params[1], params[2], params[3]);
      s2_ += lens_s2(k_[j]-l_[i], params);
    }
  s2_ = reduce_mod_ZZ(s2_);
  // printf("=> s2(E)  = %lld/%lld\n",s2_.numerator(),s2_.denominator());
}


//////////////////////////////////////////////////
// Lens space invariants:
rational<INT_KS> Space::lens_s2(INT_P p, array<INT_P,4> param)
{
  using boost::math::cos_pi;
  using boost::math::sin_pi;
  FLOAT_KS a = 0;
  for(INT_P k = 1; k < abs(p); ++k)
    {
      FLOAT_KS s = 
	(cos_pi((FLOAT_KS)2*k/abs(p)) - 1)
	/sin_pi((FLOAT_KS)k*param[0]/p)
	/sin_pi((FLOAT_KS)k*param[1]/p)
	/sin_pi((FLOAT_KS)k*param[2]/p)
	/sin_pi((FLOAT_KS)k*param[3]/p);
      a += s;
    }
  /*std::cout << "    ";
  std::cout << std::setprecision(std::numeric_limits<FLOAT_KS>::max_digits10)
  	    << a*45 
  	    << std::endl; // print a*/
  INT_KS rounded_45_a = (INT_KS)boost::math::round(a*45);
  //printf("     round(...) = %lld\n",rounded_45_a);
  rational<INT_KS> s2(rounded_45_a,45*16*p);
  //printf("     s2 = %lld/%lld\n",s2.numerator(),s2.denominator());
  return s2;
}

//////////////////////////////////////////////////

Space::comp Space::compareHomotopyType(const Space& E1, const Space& E2)
{ 
  if (abs(E1.r_ ) > abs(E2.r_ ))  return comp::GREATER;
  if (abs(E1.r_ ) < abs(E2.r_ ))  return comp::SMALLER;
  if (abs(E1.s_ ) > abs(E2.s_ ))  return comp::GREATER;
  if (abs(E1.s_ ) < abs(E2.s_ )) return comp::SMALLER;
  if (abs(E1.M1_) > abs(E2.M1_)) return comp::GREATER;
  if (abs(E1.M1_) < abs(E2.M1_)) return comp::SMALLER;
  if (    E1.M2_  >     E2.M2_ ) return comp::GREATER;
  if (	  E1.M2_  <     E2.M2_ ) return comp::SMALLER;
  if (sign(E1.M1_)*sign(E1.s_) > sign(E2.M1_)*sign(E2.s_)) return comp::GREATER;
  if (sign(E1.M1_)*sign(E1.s_) < sign(E2.M1_)*sign(E2.s_)) return comp::SMALLER;
  return comp::EQUAL;
}

Space::comp Space::compareHomotopyType_using_KS(const Space& E1, const Space& E2)
{ 
  if (abs(E1.r_) > abs(E2.r_)) return comp::GREATER;
  if (abs(E1.r_) < abs(E2.r_)) return comp::SMALLER;
  if (abs(E1.s_) > abs(E2.s_)) return comp::GREATER;
  if (abs(E1.s_) < abs(E2.s_)) return comp::SMALLER;
  // If both KS-invariants are unkown, return MAYBE_EQUAL;
  // otherwise, the unknown KS-invariant is   MAYBE_GREATER 
  // than the known one.			
  if (E1.s22() == KS_UNKNOWN || E1.s22() == KS_UNCOMPUTABLE)
    if (E2.s22() == KS_UNKNOWN || E2.s22() == KS_UNCOMPUTABLE)
      return comp::MAYBE_EQUAL;
    else return comp::MAYBE_GREATER;
  if (E2.s22() == KS_UNKNOWN || E2.s22() == KS_UNCOMPUTABLE)
    return comp::MAYBE_SMALLER;
  if (abs(E1.s22()) > abs(E2.s22())) return comp::GREATER;
  if (abs(E1.s22()) < abs(E2.s22())) return comp::SMALLER;
  if (sign(E1.s22())*sign(E1.s_) > sign(E2.s22())*sign(E2.s_)) return comp::GREATER;
  if (sign(E1.s22())*sign(E1.s_) < sign(E2.s22())*sign(E2.s_)) return comp::SMALLER;
  return comp::EQUAL;
}

Space::comp Space::compareTangentialHomotopyType(const Space& E1, const Space& E2)
{
  comp homotopy = compareHomotopyType(E1,E2);
  if (homotopy != comp::EQUAL) return homotopy;
  if (abs(E1.p1_) > abs(E2.p1_)) return comp::GREATER;
  if (abs(E1.p1_) < abs(E2.p1_)) return comp::SMALLER;
  return comp::EQUAL;
}

Space::comp Space::compareHomeomorphismType(const Space& E1, const Space& E2)
{
  comp tangential_homotopy = compareTangentialHomotopyType(E1,E2);
  if (tangential_homotopy != comp::EQUAL) return tangential_homotopy;
  // If KS-invariants of both spaces are unkown, return MAYBE_EQUAL;
  // otherwise, the unknown KS-invariant is MAYBE_GREATER 
  // than the known one.			
  if (E1.s2_ == KS_UNKNOWN || E1.s2_ == KS_UNCOMPUTABLE)
    if (E2.s2_ == KS_UNKNOWN || E2.s2_ == KS_UNCOMPUTABLE)
      return comp::MAYBE_EQUAL;
    else return comp::MAYBE_GREATER;
  if (E2.s2_ == KS_UNKNOWN || E2.s2_ == KS_UNCOMPUTABLE)
    return comp::MAYBE_SMALLER;
  if (abs(E1.s2_) > abs(E2.s2_)) return comp::GREATER;
  if (abs(E1.s2_) < abs(E2.s2_)) return comp::SMALLER;
  if (sign(E1.s2_)*sign(E1.s_) > sign(E2.s2_)*sign(E2.s_)) return comp::GREATER;
  if (sign(E1.s2_)*sign(E1.s_) < sign(E2.s2_)*sign(E2.s_)) return comp::SMALLER;
  return comp::EQUAL;
}
