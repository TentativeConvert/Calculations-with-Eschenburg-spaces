#pragma once
#include <array>               // std::array;
#include <boost/rational.hpp>  // boost::rational;

class Space {
 private:
  std::array<long,3> k_;
  std::array<long,3> l_;
  long r_;  // = r(k,l); it can be positive or negative.
            // (in [CEZ07], "r := |r(k,l)|")
  long s_;  // = s(k,l) modulo |r(k,l)|
  long p1_; // = p_1(k,l) modulo |r(k,l)|

  boost::rational<long long> s2_;
  boost::rational<long long> s22_;
  // KS-invariants s2 and s22 take values in interval (-1/2,1/2]
  static const boost::rational<long long> KS_UNKNOWN;      // (set to -1/1 in .cpp)
  static const boost::rational<long long> KS_UNCOMPUTABLE; // (set to 1/1 in .cpp)
  // KS-invariants s2, s22 are uncomputable if condition C is not satisfied

  int good_col_or_row;
  static const int GOOD_CoR_UNKNOWN = -1;
  static const int GOOD_CoR_NONEXISTENT = 6;
  // 0,1,2  = col 1,2,3 satisfies condition C
  // 3,4,5  = row 1,2,3 satisfies condition C
  void find_good_col_or_row(void);
  int gcdA(int i, int j, int ii, int jj);
  void compute_s2_from_row(int j);
  void compute_s2_from_col(int j);
  boost::rational<long long> lens_s2(long p, std::array<long,4> param);

public:
  // Constructor:
  Space(std::array<long,3> kkk, std::array<long,3>);

  // Methods that change class members:
  bool test_condition_C(); 
  bool compute_KS_invariants();  // same return value as conditionC, but also computes s22 & s2

  // Getters that don't change class members:
  bool is_positively_curved_space(void) const;
  const std::array<long,3>& k() const { return k_; }
  const std::array<long,3>& l() const { return k_; }
  const long& r() const { return r_; }
  const long& s() const { return s_; }
  const long& p1() const { return p1_; }
  const boost::rational<long long>& s22() const {return s22_; }
  const boost::rational<long long>& s2() const {return s2_; }

  // Other methods that don't change class members:
  // (note that wrong values of s2 & s22 will be printed
  //  if they haven't yet been computed)
  void print(FILE* file) const;
  void print(void) const;

  // static methods:
  enum class comp {EQUAL, SMALLER, GREATER};
  static comp compareBasicType(const Space& E1, const Space& E2);
  static comp compareHomotopyType(const Space& E1, const Space& E2);
  static comp compareTangentialHomotopyType(const Space& E1, const Space& E2);
  static comp compareHomeomorphismType(const Space& E1, const Space& E2);
};


