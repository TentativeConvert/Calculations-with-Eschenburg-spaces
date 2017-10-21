#include "aux_feedback.h"
using namespace std::chrono;
#include <algorithm> // for std::max

void Feedback::update_percent(int percent) const
{
  system_clock::time_point t = system_clock::now();
  float duration = (float)duration_cast<milliseconds>(t - t_start_).count()/1000;
  printf(" %3d%%  (%.2f s)\r", percent, duration);
  fflush(stdout);
}

void Feedback::start(size_t max)
{ 
  t_start_ = system_clock::now();
  step_size_ = std::max((long)max/100,(long)1);
  update_percent(0);
};

void Feedback::update(size_t progress) const
{
  if(progress % step_size_ == 0) 
    update_percent(std::min((int)(progress/step_size_),(int)99));
}; 

void Feedback::finish() const
{
  update_percent(100);
  printf("\n");
}
