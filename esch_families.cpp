#include "esch_families.h"
using std::deque;
using std::vector;

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
