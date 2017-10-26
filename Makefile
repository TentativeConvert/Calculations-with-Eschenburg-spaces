# Makefile taken from 
# https://stackoverflow.com/a/2481326/3611932
#
CXX=g++
RM=rm -f
CPPFLAGS=-std=c++11 -I=/usr/local/include/boost_1_65_1 -O3
LDFLAGS=-std=c++11 -I=/usr/local/include/boost_1_65_1 -O3

SRCS=main.cpp esch_space.cpp esch_tuples.cpp esch_generate.cpp aux_math.cpp aux_feedback.cpp
OBJS=$(subst .cpp,.o,$(SRCS))

all: main

main: $(OBJS)
	$(CXX) $(LDFLAGS) -o main $(OBJS)

depend: .depend

.depend: $(SRCS)
	$(RM) ./.depend
	$(CXX) $(CPPFLAGS) -MM $^>>./.depend;

clean:
	$(RM) $(OBJS)

distclean: clean
	$(RM) *~ .depend

include .depend
