CXX?=g++
CFLAGS=-Wall -Werror -Wextra -Wno-ignored-attributes -Wunused-parameter -Wpointer-arith -DLARGEFILE_SOURCE -D_LARGEFILE64_SOURCE -D_FILE_OFFSET_BITS=64 -funroll-loops -O3
COMPILATION_RESOURCES=-ftemplate-depth=2000 -fconstexpr-depth=2000
CPPFLAGS=-I ../src
LD=${CXX}
LDLIBS=-lm -lz -lpthread -lstdc++
CXXFLAGS=-std=c++17

GTTL?=../../gttl
GTTLFLAGS=-I ${GTTL}/src

TMPFILE=temp.txt

QMER=sorted_q_mer
SPACED_SEEDS=spaced_seeds

all:${QMER}.x ${SPACED_SEEDS}.x

test: test_sorted_qmer_nuc_3 test_sorted_qmer_prot_2 test_spaced_seeds

test_sorted_qmer_nuc_3:${QMER}.x
	@./${QMER}.x -q 3 -s >${TMPFILE}
	@diff ./ref/sorted_qmer_nuc_3.txt ${TMPFILE}
	@${RM} ${TMPFILE}
	@echo "$@ passed."

test_sorted_qmer_prot_2:${QMER}.x
	@./${QMER}.x -p -q 2 -s >${TMPFILE}
	@diff ./ref/sorted_qmer_prot_2.txt ${TMPFILE}
	@${RM} ${TMPFILE}
	@echo "$@ passed."

test_spaced_seeds:${SPACED_SEEDS}.x
	@./${SPACED_SEEDS}.x -p -s ../testdata/uniprot_sprot.fasta >${TMPFILE}
	@diff ./ref/test_uniprot_sprot.txt ${TMPFILE}
	@${RM} ${TMPFILE}
	@echo "$@ passed."

%.o:%.cpp
	$(CXX) -o $@ -c $< ${CXXFLAGS} ${COMPILE_CONST} ${COMPILATION_RESOURCES} ${CPPFLAGS} ${GTTLFLAGS} -MT $@ -MMD -MP -MF $(@:.o=.d)

%.x:%.o
	$(CXX) ${LDFLAGS} $< -o $@ ${LDLIBS}

.PHONY:clean
clean:
	${RM} -r *.[odx] *.x.dSYM/ __pycache__ temp.txt
