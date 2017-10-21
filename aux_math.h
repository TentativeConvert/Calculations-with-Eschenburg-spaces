#pragma once

#include <cmath>     // std::abs,sqrt
#include <algorithm> // std::min, max

#include <boost/math/common_factor.hpp>
// boost::math::gcd;

#include <boost/rational.hpp>
// boost::rational;

long square(long long a);
long signed_mod (long a, long base);
long absolute_mod (long a, long base);
boost::rational<long long> reduce_mod_ZZ(const boost::rational<long long>& q);
