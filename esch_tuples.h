#pragma once

#include "esch_space.h"
#include <deque>
#include <functional>
#include <cstddef> // for std::size_t

class SpaceTuple : public std::deque< class Space >
{
 public:
  std::size_t compute_KS_invariants(); // (return value = number of spaces for which condition C fails)

  void print(FILE* file) const;
  bool operator<(const SpaceTuple& otherfam) const {return (this->size() > otherfam.size());}
};

class SpaceTupleList : public std::deque< SpaceTuple >
{
 public:
  // CONSTRUCTOR 1 (implemented in esch_generate.cpp):
  // fills list with pairs of spaces whose invariants r & s coincide
  SpaceTupleList(const INT_R& max_R);

  // CONSTRUCTOR 2 (implement in esch_tuples.cpp):
  // filters tuples from original_list according to the supplied function;
  // original list may be resorted, so cannot be passed as const
  SpaceTupleList(SpaceTupleList& original_list, 
		 std::function<Space::comp(const Space& E1, const Space& E2)> compareFunction, 
		 const char* description);

  // Other methods:
  std::size_t compute_KS_invariants(); // (return value = number of spaces for which condition C fails)
  void print(const char* filename, const char* description);
};
