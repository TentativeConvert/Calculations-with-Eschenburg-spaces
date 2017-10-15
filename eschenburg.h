#include <array>
using std::array;

#include <boost/rational.hpp>
using boost::rational;

class Space {
 public:
  // Setters:
  void setParameters(array<long,3> kkk, array<long,3>);  //xx

  // Getters that change class members:
  bool conditionC(); 
  const rational<long long>& s22(); 
  const rational<long long>& s2(); 

  // Getters that don't change class members:
  bool is_positively_curved_space(void);
  const array<long,3>& k() const { return k_; }
  const array<long,3>& l() const { return k_; }
  const long& s() const { return s_; }
  const long& p1() const { return p1_; }

  // Other methods that don't change class members:
  void print_for_human(FILE* file) const;
  void print_for_maple(FILE* file) const;
  
 private:
  array<long,3> k_;
  array<long,3> l_;
  long r_;  // = r(k,l); it can be positive or negative.
            // (in [CEZ07], "r := |r(k,l)|")
  long s_;  // = s(k,l) modulo |r(k,l)|
  long p1_; // = p_1(k,l) modulo |r(k,l)|

  rational<long long> s2_;
  // values in interval (-1/2,1/2]
  // -1 = not yet computed
  //  1 = condition C not satisfied
  rational<long long> s22_;
  // values in interval (-1/2,1/2]
  // -1 = not yet computed

  int good_col_or_row;
  // -1     = not checked
  // 0,1,2  = col 1,2,3 satisfies condition C
  // 3,4,5  = row 1,2,3 satisfies condition C
  // 6      =  condition C not satisfied
  int gcdA(int i, int j, int ii, int jj);
  void find_good_col_or_row(void);
  void compute_s2_row(int j);
  void compute_s2_col(int j);
};
bool compare_spaces(struct Space E1, struct Space E2);

