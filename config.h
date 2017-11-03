#pragma once

#define MAX_TUPLES_PER_TUPLESIZE_PER_FILE 100

//////////////////////////////////////////////////
// DATA TYPES
// The following types should be large enough to work for |r| < 10^6.

#include <boost/multiprecision/cpp_int.hpp>
#include <boost/cstdint.hpp>
// http://www.boost.org/doc/libs/1_65_1/libs/multiprecision/doc/html/boost_multiprecision/tut/ints/cpp_int.html

typedef boost::int_t<32>::exact INT_P;  // should hold all possible values of sqrt(r) 
typedef boost::int_t<64>::exact INT_R; // should hold all possible values of r
typedef boost::multiprecision::int128_t INT_KS; 
                    // should hold integers necessary to compute KS-invariants

//#
//typedef int_fast32_t INT_P; 
//typedef int_fast32_t INT_R; 
//typedef int_fast64_t INT_KS;

typedef long double FLOAT_KS;
//#include <boost/multiprecision/cpp_dec_float.hpp>
//#typedef boost::multiprecision::cpp_dec_float_50 FLOAT_KS;

//////////////////////////////////////////////////
// LIBRARIES
// boost chrono library available?
#define CHRONO_AVAILABLE

