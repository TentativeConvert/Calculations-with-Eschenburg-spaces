#include "aux_math.h"

//#include <boost/math/special_functions/round.hpp>

INT_KS square(INT_KS a){ return a*a; }
INT_R signed_mod (INT_R a, INT_R base)
{
  long remainder = a % base;
  if(2*remainder > base)  // remainder > base/2, in "integer language"
    remainder -= base;
  else if(2*remainder <= -base) // remainder < -base/2, in "integer language"
    remainder += base;
  return remainder;
}
INT_KS signed_mod (INT_KS a, INT_KS base) // exact copy of function above, for bigger integers
{
  INT_KS remainder = a % base;
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
    remainder += abs(base);
  return remainder;
}

rational<INT_KS> reduce_mod_ZZ(const rational<INT_KS>& q)
{
  // Input:   q  in QQ 
  // Output:  q' in (-1/2, 1/2] in QQ such that 
  //          q' = q in QQ/ZZ.
  INT_KS q_n = q.numerator();
  INT_KS q_d = q.denominator();
  if(q_d < 0) // probaby unnecessary, but depends on implementation of boost::rational
    {
      q_n = -q_n;
      q_d = -q_d;
    }
  q_n = signed_mod(q_n,q_d);  
  return rational<INT_KS>(q_n,q_d);
}

int sign(const INT_P &i){
  if (i > 0) return 1;
  if (i < 0) return -1;
  return 0;
}
int sign(const INT_R &i){
  if (i > 0) return 1;
  if (i < 0) return -1;
  return 0;
}
int sign(const INT_KS &i){
  if (i > 0) return 1;
  if (i < 0) return -1;
  return 0;
}

int sign(const rational<INT_KS> &r){
  // implementation is really for rationals mod ZZ
  // (We need this only for the invariants s22 and s2.)
  //
  //  r =  1/4          ---> sign(r) = +1
  //  r = -1/4          ---> sign(r) = -1
  //  r = 0, 1/2, -1/2  ---> sign(r) =  0
  rational<INT_KS> r_ = reduce_mod_ZZ(r);
  if (0 < r && r < rational<INT_KS>(1,2)) return 1;
  if (0 > r && r > -rational<INT_KS>(1,2)) return -1;
  return 0;
}
