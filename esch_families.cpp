#include "esch_families.h"
using std::deque;
using std::vector;

#include "aux_feedback.h"

#include <algorithm>
using std::stable_sort;

void Space_family::print_for_human(FILE* file) const
{
  for(const Space& E : *this)
    {
      E.print_for_human(file);
      fprintf(file,"\n ");
    }
  // printf("%d",(int)(this->size()));
  return;
}

void Deque_of_Space_families::print(const char* filename) const
{
  printf("Writing to file ...");
  FILE *file = fopen(filename, "w");  // r, s agree
  if (file == NULL)
    printf("Error opening file!\n");
  for(size_t c = counters_.size(); c >= 2; --c) 
    // note that size_t is unsigned, so can't start at counters_size()-1!
    {
      fprintf(file,"Found %4ld families with %ld members.\n", counters_[c-1], (long)c);
    }
  for(size_t i = 0; i < this->size(); ++i)
    {
      fprintf(file,"\nFamily %ld:\n ", (long)i+1);
      this->at(i).print_for_human(file);
    }
  fclose(file);
  printf(" ... done.\n");
}

size_t Space_family::compute_KS_invariants() 
// (return value = number of spaces for which condition C fails)
{
  size_t failures = 0;
  for(Space& E : *this)
    if (! E.compute_KS_invariants())
      ++failures;
}

size_t Deque_of_Space_families::compute_KS_invariants()
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

void Deque_of_Space_families::sort_and_count_families()
{ 
  // sort families by size:
  stable_sort(this->begin(),this->end());
  // count number of families with 2 members, with 3 members, etc.:
  counters_.clear();
  for(Space_family F : *this)
    {
      size_t s = F.size();
      if (s > counters_.size()) 
	counters_.resize(s,0);
      ++counters_[s-1];
    }
}

void Deque_of_Space_families::filter(Deque_of_Space_families& original_families, Deque_of_Space_families& new_families, std::function<Space::comp(const Space& E1, const Space& E2)> compareFunction)
{
  Feedback feedback;
  feedback.start(original_families.size());
  for(std::size_t i = 0; i < original_families.size(); ++i)
    {
      feedback.update(i);
      Space_family& F = original_families[i];

      /*struct MyComp {
	//bool operator() (const Space& E1, const Space& E2){ return (compareFunction(E1,E2) > 0); }
	bool operator() (const Space& E1, const Space& E2){ return (E1.p1() > E2.p1()); }
	} myComp;*/
      //sort(F.begin(),F.end(),Filter_rs_to_homotopy_equivalent::sort);
      //sort(F.begin(),F.end(),myComp);
      sort(F.begin(),F.end(),[&compareFunction](const Space& E1, const Space& E2) -> bool
	   {
	     return (compareFunction(E1,E2) == Space::comp::GREATER);
	   } );

      for(std::size_t i1 = 0; i1 < F.size(); )// i1 is incremented indirectly via i2
      {
	std::size_t i2 = i1+1;
	//while (i2 < F.size() && Filter_rs_to_homotopy_equivalent::equal(F[i1],F[i2]))
	while (i2 < F.size() && compareFunction(F[i1],F[i2]) == Space::comp::EQUAL)
	  ++i2;
	if(i2 > i1+1)
	  {
	    //all_spaces with indexes i1,...,i2 define spaces with the same invariants
            Space_family new_family;
	    for(std::size_t j = i1; j < i2; ++j)
	      {
		new_family.push_back(F[j]);
	      }
	    new_families.push_back(new_family);
	  }
	i1 = i2;
      }
    }
  feedback.finish();
}
