#include "config.h"      // defines MAX_TUPLES_PER_TUPLESIZE_PER_FILE
#include "esch_tuples.h"
using std::deque;
#include <vector>
using std::vector;
#include "aux_feedback.h"
#include <algorithm>

//////////////////////////////////////////////////

void SpaceTuple::print(FILE* file) const
{
  for(const Space& E : *this)
    {
      E.print(file);
    }
  return;
}

void SpaceTupleList::print(const char* filename)
{
  FILE *file = fopen(filename, "w");
  if (file == NULL) printf("Error opening file!\n");
  
  // tuples should already by sorted by size, in descending order, but let's make sure:
  printf("Sorting and counting ...\r");
  std::stable_sort(this->begin(),this->end());

  // count total number of pairs, triples, tuples of length 3, ...:
  // counter_[0] = number of singletons (already known)
  // counter_[1] = number of pairs
  // conuter_[2] = ...
  size_t counter_total = singletons;
  vector< std::size_t > counter(1);  
  counter[0] = singletons;
  for(SpaceTuple F : *this)
    {
      std::size_t s = F.size();
      if (s > counter.size()) 
	counter.resize(s,0);
      ++counter[s-1];
      ++counter_total;
    }
  
  // print "statistics" to screen and to file (number of tuples of each length):
  printf(       (">> %9lld different " + description + " in this range\n").c_str(),(long long)counter_total);
  fprintf(file, (">> %9lld different " + description + " in this range\n").c_str(),(long long)counter_total);
  
  for(std::size_t c = 1; c <= counter.size(); ++c) 
    {
      std::string text;
      if (c == 1)
	text = ">> %9lld " + description + " defined by just %ld parameter vector (k1,k2,k3,l1,l2,l3)\n";
      else if (c >= 2)
	text = ">> %9lld " + description + " defined by %ld different parameter vectors\n";
      printf(text.c_str(), counter[c-1], (long long)c);
      fprintf(file,text.c_str(), counter[c-1], (long long)c);
    }
  
  // Print tuples to file, grouped by tuple length
  printf("Writing to file %s ...",filename);
  for(std::size_t c = counter.size(); c >= 2; --c)
    {
      fprintf(file,"\n\n\n###################### Tuples of length %ld ######################\n",(long)c);
      // As tuples are ordered by size in DESCENDING order, 
      // start position = number of tuples of larger sizes:
      std::size_t start = 0;
      for(std::size_t i = c; i < counter.size(); ++i)
	start += counter[i];
      for(std::size_t i = 0; i < counter[c-1]; ++i)
	{
	  fprintf(file," \nTuple %ld: \n", (long)i+1);
	  this->at(start + i).print(file);
	  if  (i+11 >= MAX_TUPLES_PER_TUPLESIZE_PER_FILE && i+11 < counter[c-1])
	    {
	      fprintf(file,"\n .\n .\n . skipping some tuples of length %ld \n .\n .\n",(long)c);
	      i = counter[c-1]-11;
	    }
	}
      fprintf(file,"\n\n\n");
    }
  fclose(file);
  printf(" ... done.\n\n");
}

std::size_t SpaceTuple::compute_KS_invariants() 
// (return value = number of spaces for which condition C fails)
{
  std::size_t failures = 0;
  for(Space& E : *this)
    if (! E.compute_KS_invariants())
      ++failures;
  return failures;
}

std::size_t SpaceTuple::test_condition_C()
// (return value = number of spaces for which condition C fails)
{
  std::size_t failures = 0;
  for(Space& E : *this)
    if (! E.test_condition_C())
      ++failures;
  return failures;
}

std::size_t SpaceTupleList::compute_KS_invariants()
// (return value = number of spaces for which condition C fails)
{
  printf("Computing Kreck-Stolz-invariants s2 and s22 for non-singletons ...\n");
  std::size_t failures = 0;
  Feedback feedback;
  feedback.start(this->size());
  for(std::size_t i = 0; i < this->size(); ++i)
    {
      feedback.update(i);
      failures += this->at(i).compute_KS_invariants();
    }
  feedback.finish();
  printf("Invariants computed.  ");
  if(failures == 0) 
    printf("All spaces satisfy condition C.\n\n");
  else
    printf("\nWARNING: Condition C fails for %ld non-singletons.\n\n",(long)failures);
  return failures;
}

std::size_t SpaceTupleList::test_condition_C()
// (return value = number of spaces for which condition C fails)
{
  printf("Checking conditin C for non-singletons ...\n");
  std::size_t failures = 0;
  Feedback feedback;
  feedback.start(this->size());
  for(std::size_t i = 0; i < this->size(); ++i)
    {
      feedback.update(i);
      failures += this->at(i).test_condition_C();
    }
  feedback.finish();
  if(failures == 0) 
    printf("All non-singletons satisfy condition C.\n\n");
  else
    printf("\nWARNING: Condition C fails for %ld non-singletons.\n\n",(long)failures);
  return failures;
}

bool sort(const SpaceTuple& T1, const SpaceTuple& T2,
	  std::function<Space::comp(const Space& E1, const Space& E2)> compareFunction)
{
  if (T1 > T2) return true; // comparison by size
  if (T1 < T2) return false;
  if (!(T1.empty() || T2.empty()))
    {
      return (compareFunction(T1[0],T2[0]) == Space::comp::GREATER);
    }
  return false;
}


SpaceTupleList::SpaceTupleList(SpaceTupleList& original_list, 
			       std::function<Space::comp(const Space& E1, const Space& E2)> compareFunction,
			       std::string description)
{
  this->description = description;
  printf("%s\n",("Looking for " + description + " ...").c_str());
  singletons = original_list.singletons;
  Feedback feedback;
  feedback.start(original_list.size());
  for(std::size_t i = 0; i < original_list.size(); ++i)
    {
      feedback.update(i);
      SpaceTuple& T = original_list[i];

      std::stable_sort(T.begin(),T.end(),[&compareFunction](const Space& E1, const Space& E2) -> bool
	   {
	     return (compareFunction(E1,E2) == Space::comp::GREATER);
	   } );

      for(std::size_t i1 = 0; i1 < T.size(); )// i1 is incremented indirectly via i2
      {
	std::size_t i2 = i1+1;
	while (i2 < T.size() && compareFunction(T[i1],T[i2]) == Space::comp::EQUAL)
	  ++i2;
	if(i2 > i1+1)
	  {
	    //all_spaces with indexes i1,...,i2 define spaces with the same invariants
            SpaceTuple new_tuple;
	    for(std::size_t j = i1; j < i2; ++j)
		new_tuple.push_back(T[j]);
	    this->push_back(new_tuple);
	  }
	else if (i2 = i1+1)
	  ++singletons;
	i1 = i2;
      }
    }
  feedback.finish();
  printf("Sorting ...\r");
  // sort tuples by size, in DESCENDING ORDER, and by compareFunction
  std::stable_sort(this->begin(),this->end(),[&compareFunction](const SpaceTuple& T1, const SpaceTuple& T2) -> bool
	    {
	      return sort(T1,T2,compareFunction);
	    } );
  printf("Done.      \r");
}
