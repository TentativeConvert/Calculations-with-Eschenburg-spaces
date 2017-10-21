#pragma once

#include <chrono>

class Feedback
{
 private:
  std::chrono::system_clock::time_point t_start_;
  size_t step_size_;
 
 public:
  void start(size_t max);
  void update(size_t progress) const;    // passes progress/max to update_percent
  void update_percent(int percent) const; 
  void finish(void) const;
};
