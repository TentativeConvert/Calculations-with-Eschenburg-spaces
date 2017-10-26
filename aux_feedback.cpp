#include "aux_feedback.h"
#include <algorithm> // for std::max

void Feedback::update_percent(int percent) const
{
#ifdef CHRONO_AVAILABLE
  using namespace std::chrono; 
  system_clock::time_point t = system_clock::now();
  float duration = (float)duration_cast<milliseconds>(t - t_start_).count()/1000;
  printf(" %3d%%  (%.2f s)\r", percent, duration);
#else 
  printf(" %3d%%\r", percent);
#endif
  fflush(stdout);
}

void Feedback::start(std::size_t max)
{ 
#ifdef CHRONO_AVAILABLE
  t_start_ = std::chrono::system_clock::now();
#endif
  step_size_ = std::max((long)max/100,(long)1);
  update_percent(0);
};

void Feedback::update(std::size_t progress) const
{
  if(progress % step_size_ == 0) 
    update_percent(std::min((int)(progress/step_size_),(int)99));
}; 

void Feedback::finish() const
{
  update_percent(100);
  printf("\n");
}
