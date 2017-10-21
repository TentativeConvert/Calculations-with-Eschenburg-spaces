#include "aux_math.h"
//#include <cmath>
//#include <algorithm>
using std::abs;
using std::sqrt;
using std::min;
using std::max;
using boost::rational;

#include <boost/math/special_functions/round.hpp>
using boost::math::gcd;

long square(long long a){ return a*a; }
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
    remainder += abs(base);
  return remainder;
}
rational<long long> reduce_mod_ZZ(const rational<long long>& q)
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
  return rational<long long>(q_n,q_d);
}
