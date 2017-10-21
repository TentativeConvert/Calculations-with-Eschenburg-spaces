#pragma once

#include "esch_space.h"
#include <deque>
#include <vector>

class Space_family : public std::deque< class Space >
{
 public:
  void print_for_human(FILE* file) const;
  bool operator<(const Space_family& otherfam) const {return (this->size() > otherfam.size());}
};

class Deque_of_Space_families : public std::deque< Space_family >
{
 private:
  std::vector< long > counters_;
  
 public:
  void sort_and_count_families();
  void print(const char* filename) const;
};
