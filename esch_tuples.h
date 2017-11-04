#pragma once

#include "esch_space.h"
#include <string>
#include <deque>
#include <functional>
#include <cstddef> // for std::size_t

class SpaceTuple : public std::deque< class Space >
{
 public:
  std::size_t test_condition_C(); // return value = number of spaces for which condition C fails
  std::size_t compute_KS_invariants(); // return value as for test_condition_C()
  void print(FILE* file) const;

//bool operator<(const SpaceTuple& otherfam) const {return (this->size() > otherfam.size());}
  friend bool operator<(const SpaceTuple& T1, const SpaceTuple& T2) {return (T1.size() > T2.size());}
  friend bool operator>(const SpaceTuple& T1, const SpaceTuple& T2) {return T2 < T1;}
  friend bool sort(const SpaceTuple& T1, const SpaceTuple& E2,
	  std::function<Space::comp(const Space& E1, const Space& E2)> compareFunction);
};




class SpaceTupleList : public std::deque< SpaceTuple >
{
 public:
  std::size_t singletons;  // number of 1-tuples (1-tuples themselves are not stored)
  std::string description;
  
  // CONSTRUCTOR 1 (implemented in esch_generate.cpp):
  // fills list with pairs of spaces whose invariants r & s coincide
  SpaceTupleList(const INT_R& max_R);

  // CONSTRUCTOR 2 (implement in esch_tuples.cpp):
  // filters tuples from original_list according to the supplied function;
  // original list may be resorted, so cannot be passed as const
  SpaceTupleList(SpaceTupleList& original_list, 
		 std::function<Space::comp(const Space& E1, const Space& E2)> compareFunction, 
		 std::string description);

  // Other methods:
  std::size_t test_condition_C();// return value = number of spaces for which condition C fails
  std::size_t compute_KS_invariants();// return value of as for test_condition_C()
  void print(const char* filename);
};
