#pragma once
#define DEFAULT_MAX_TUPLES_PER_TUPLE_SIZE_PER_FILE 1010

//////////////////////////////////////////////////
// DATA TYPES
// see docs/limits.pdf

#include <cstdint>
#include <boost/multiprecision/cpp_int.hpp>

typedef std::int_least32_t INT_P;                // should hold all possible values of sqrt(r) 
typedef std::int_least64_t INT_R;                // should hold all possible values of r
typedef boost::multiprecision::int128_t INT_KS;  // should hold integers necessary to compute KS-invariants
typedef long double FLOAT_KS;

//#include <boost/multiprecision/cpp_dec_float.hpp>
//#typedef boost::multiprecision::cpp_dec_float_50 FLOAT_KS;
/* Integer types supplied by boost::multiprecision:

   http://www.boost.org/doc/libs/1_65_1/libs/multiprecision/doc/html/boost_multiprecision/tut/ints/cpp_int.html

*/


//////////////////////////////////////////////////
// LIBRARIES
// boost chrono library available?
#define CHRONO_AVAILABLE

