CXX?=g++
CFLAGS=-Wall -Werror -Wextra -Wno-ignored-attributes -Wunused-parameter -Wpointer-arith -DLARGEFILE_SOURCE -D_LARGEFILE64_SOURCE -D_FILE_OFFSET_BITS=64 -funroll-loops -O3
CPPFLAGS=-I ../src
LD=${CXX}
LDLIBS=-lm -lz -lpthread -lstdc++
CXXFLAGS=-std=c++17

GTTL?=../../gttl
GTTLFLAGS=-I ${GTTL}/src

%.o:%.cpp
	$(CXX) -o $@ -c $< ${CXXFLAGS} ${CFLAGS} ${GTTLFLAGS} -MT $@ -MMD -MP -MF $(@:.o=.d)

%.x:%.o
	$(CXX) ${LDFLAGS} $< -o $@ ${LDLIBS}

.PHONY:clean
clean:
	${RM} -r *.[odx] *.x.dSYM/ __pycache__