#pragma once
#include <chrono>
#include "esch_space.h"
#include "esch_families.h"

void generate_rs_families(Deque_of_Space_families& families_rs, const long& R, const std::chrono::system_clock::time_point& t1);
