#pragma once
#include "config.h"
#include <cstddef> // for std::size_t

#ifdef CHRONO_AVAILABLE
#include <chrono>
#endif

class Feedback
{
 private:
#ifdef CHRONO_AVAILABLE
  std::chrono::system_clock::time_point t_start_;
#endif
  std::size_t step_size_;
  
 public:
  void start(std::size_t max);
  void update(std::size_t progress) const;    // passes progress/max to update_percent
  void update_percent(int percent) const; 
  void finish(void) const;
};
