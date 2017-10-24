#pragma once

#include <cmath>     
using std::abs; //*
using std::sqrt;//*

#include <algorithm>
using std::max; //* 
using std::min; //*

#include <boost/math/common_factor.hpp>
using boost::math::gcd; //*

long square(long long a);
long signed_mod (long a, long base);
long absolute_mod (long a, long base);
int sign(const long &i);

// rationals:
#include <boost/rational.hpp>
using boost::rational; //*
//template<class Typ> rational<Typ> reduce_mod_ZZ(const rational<Typ> &q);
rational<long long> reduce_mod_ZZ(const rational<long long> &q);
int sign(const rational<long long> &r);



/* 

* Using "using" directives and declarations in header files is generally discouraged, for good reason.  However, here I am simply including std::abs, std::sqrt, ...  in place of custom implementations of these functions.  The effect on a file with an 

  #include "aux_math.h"

directive is exactly the same as if I had implemented by own functions abs(...), sqrt(...), in aux_math.h.

*/
