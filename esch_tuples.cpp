#include "config.h"      // defines MAX_TUPLES_PER_TUPLESIZE_PER_FILE
#include "esch_tuples.h"
using std::deque;
#include <vector>
using std::vector;
#include "aux_feedback.h"
#include <algorithm>
using std::stable_sort;

//////////////////////////////////////////////////

void SpaceTuple::print(FILE* file) const
{
  for(const Space& E : *this)
    {
      E.print(file);
    }
  return;
}

void SpaceTupleList::print(const char* filename, const char* description)
{
  FILE *file = fopen(filename, "w");
  if (file == NULL) printf("Error opening file!\n");
  
  // sort tuples by size, in DESCENDING ORDER:
  stable_sort(this->begin(),this->end());

  // count total number of pairs, triples, tuples of length 3, ...:
  // counter_[i-1] = number of tuples of length i
  vector< long > counters_;  
  for(SpaceTuple F : *this)
    {
      size_t s = F.size();
      if (s > counters_.size()) 
	counters_.resize(s,0);
      ++counters_[s-1];
    }
  
  // print "statistics" to screen and to file (number of tuples of each length):
  printf("%s\n", description);
  fprintf(file, "%s\n", description);
  for(size_t c = counters_.size(); c >= 2; --c) 
    {
      // note that size_t is unsigned, so can't start at counters_.size()-1!
      if (c == 3)
	{
	  printf(      "  %5ld triples.\n", counters_[c-1]);
	  fprintf(file,"  %5ld triples.\n", counters_[c-1]);
	}
      else if (c == 2)
	{
	  printf(      "  %5ld pairs.\n", counters_[c-1]);
	  fprintf(file,"  %5ld pairs.\n", counters_[c-1]);
	}
      else
	{
	  printf(       "  %5ld tuples of length %ld.\n", counters_[c-1], (long)c);      
	  fprintf(file, "  %5ld tuples of length %ld.\n", counters_[c-1], (long)c);
	}
    }
  
  // Print tuples to file, grouped by tuple length
  printf("Writing to file ...");
  for(size_t c = counters_.size(); c >= 2; --c)
    {
      fprintf(file,"\n\n\n####### Tuples of length %ld #######\n",(long)c);
      // As tuples are ordered by size in DESCENDING order, 
      // start position = number of tuples of larger sizes:
      size_t start = 0;
      for(size_t i = c; i < counters_.size(); ++i)
	start += counters_[i];
      for(size_t i = 0; i < counters_[c-1] && i < MAX_TUPLES_PER_TUPLESIZE_PER_FILE; ++i)
	{
	  fprintf(file," \nTuple %ld: \n", (long)i+1);
	  this->at(start + i).print(file);
	}
      if (counters_[c-1] > MAX_TUPLES_PER_TUPLESIZE_PER_FILE)
	fprintf(file,"\n [ ... skipping %ld further tuples of length %ld ... ] \n\n", (long)(counters_[c-1]-MAX_TUPLES_PER_TUPLESIZE_PER_FILE),(long)c);
    }
  fclose(file);
  printf(" ... done.\n\n");
}

size_t SpaceTuple::compute_KS_invariants() 
// (return value = number of spaces for which condition C fails)
{
  size_t failures = 0;
  for(Space& E : *this)
    if (! E.compute_KS_invariants())
      ++failures;
}

size_t SpaceTupleList::compute_KS_invariants()
// (return value = number of spaces for which condition C fails)
{
  printf("Computing Kreck-Stolz-invariants s2 and s22 ...\n");
  size_t failures = 0;
  Feedback feedback;
  feedback.start(this->size());
  for(size_t i = 0; i < this->size(); ++i)
    {
      feedback.update(i);
      failures += this->at(i).compute_KS_invariants();
    }
  feedback.finish();
  printf("Invariants computed.  ");
  if(failures == 0) 
    printf("All spaces satisfy condition C.\n\n");
  else
    printf("\nWARNING: Condition C fails for %ld spaces\n\n",(long)failures);
}

SpaceTupleList::SpaceTupleList(SpaceTupleList& original_list, 
			       std::function<Space::comp(const Space& E1, const Space& E2)> compareFunction,
			       const char* description)
{
  printf("%s\n",description);
  Feedback feedback;
  feedback.start(original_list.size());
  for(std::size_t i = 0; i < original_list.size(); ++i)
    {
      feedback.update(i);
      SpaceTuple& F = original_list[i];

      sort(F.begin(),F.end(),[&compareFunction](const Space& E1, const Space& E2) -> bool
	   {
	     return (compareFunction(E1,E2) == Space::comp::GREATER);
	   } );

      for(std::size_t i1 = 0; i1 < F.size(); )// i1 is incremented indirectly via i2
      {
	std::size_t i2 = i1+1;
	while (i2 < F.size() && compareFunction(F[i1],F[i2]) == Space::comp::EQUAL)
	  ++i2;
	if(i2 > i1+1)
	  {
	    //all_spaces with indexes i1,...,i2 define spaces with the same invariants
            SpaceTuple new_tuple;
	    for(std::size_t j = i1; j < i2; ++j)
	      {
		new_tuple.push_back(F[j]);
	      }
	    this->push_back(new_tuple);
	  }
	i1 = i2;
      }
    }
  feedback.finish();
}
