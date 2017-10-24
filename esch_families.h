#pragma once

#include "esch_space.h"
#include <deque>
#include <vector>
#include <functional>

class Space_family : public std::deque< class Space >
{
 public:
  size_t compute_KS_invariants(); // (return value = number of spaces for which condition C fails)

  void print_for_human(FILE* file) const;
  bool operator<(const Space_family& otherfam) const {return (this->size() > otherfam.size());}
};

class Deque_of_Space_families : public std::deque< Space_family >
{
 private:
  std::vector< long > counters_;  // counter_[i-1] = number of families with i members
  
 public:
  size_t compute_KS_invariants(); // (return value = number of spaces for which condition C fails)
  void sort_and_count_families();

  void print(const char* filename) const;
  
  static void filter(Deque_of_Space_families& original_deque, Deque_of_Space_families& new_deque, std::function<Space::comp(const Space& E1, const Space& E2)> compareFunction);
};
