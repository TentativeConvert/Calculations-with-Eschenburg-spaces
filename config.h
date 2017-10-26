#pragma once

#define MAX_TUPLES_PER_TUPLESIZE_PER_FILE 100

// DATA TYPES
// The following types should be large enough to work for |r| < 10^6.

typedef long INT_R;  // should hold all possible values of r
typedef long INT_p;  // should hold all possible values of sqrt(r) 
typedef long long INT_KS; // should hold integers necessary to compute KS-invariants

#include <boost/multiprecision/cpp_dec_float.hpp>
typedef boost::multiprecision::cpp_dec_float_50 FLOAT_KS;
//#include <boost/multiprecision/float128.hpp> 
//typedef boost::multiprecision::float128 FLOAT_KS;


