#include <boost/rational.hpp>
using boost::rational;

long square(long long a);
long signed_mod (long a, long base);
long absolute_mod (long a, long base);
rational<long long> reduce_mod_ZZ(const rational<long long>& q);
