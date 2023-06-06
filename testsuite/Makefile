CXX?=g++
CFLAGS=-O3 -Wall -Wextra -Wno-ignored-attributes -Wunused-parameter -Wpointer-arith -DLARGEFILE_SOURCE -D_LARGEFILE64_SOURCE -D_FILE_OFFSET_BITS=64 -funroll-loops -O3
COMPILATION_RESOURCES=-ftemplate-depth=10000 -fconstexpr-depth=10000
CPPFLAGS=-I ../src
LD=${CXX}
LDLIBS=-lm -lz -lpthread -lstdc++
CXXFLAGS=-std=c++17

GTTL?=../../gttl
GTTLFLAGS=-I ${GTTL}/src

TMPFILE=temp.txt
REFDIR?=../testdata/ref

SQGRAM=sorted_qgram
UQGRAM=unsorted_qgram
SPACED_SEEDS=spaced_seeds
SMATRIX=score_matrix

all:${SQGRAM}.x ${SPACED_SEEDS}.x ${SMATRIX}.x ${UQGRAM}.x

test: test_sorted_qgram test_unsorted_qgram test_score_matrix test_spaced_seeds 

test_sorted_qgram:${SQGRAM}.x
	@./${SQGRAM}.x -q 3 -s >${TMPFILE}
	@diff ${REFDIR}/sorted_qgram_nuc_3.txt ${TMPFILE}
	@./${SQGRAM}.x -p -q 2 -s >${TMPFILE}
	@diff ${REFDIR}/sorted_qgram_prot_2.txt ${TMPFILE}
	@${RM} ${TMPFILE}
	@echo "$@ passed."

test_unsorted_qgram:${UQGRAM}.x
	@./${UQGRAM}.x -q 3 -s >${TMPFILE}
	@diff ${REFDIR}/unsorted_qgram_nuc_3.txt ${TMPFILE}
	@./${UQGRAM}.x -p -q 2 -s >${TMPFILE}
	@diff ${REFDIR}/unsorted_qgram_prot_2.txt ${TMPFILE}
	@${RM} ${TMPFILE}
	@echo "$@ passed."

test_score_matrix:${SMATRIX}.x
	@./${SMATRIX}.x -m 0 -q 2 >${TMPFILE}
	@diff ${REFDIR}/smatrix.tsv ${TMPFILE}
	@${RM} ${TMPFILE}
	@echo "$@ passed."

test_spaced_seeds:${SPACED_SEEDS}.x
	@./check_err.py ./${SPACED_SEEDS}.x
	@./${SPACED_SEEDS}.x -p -d 0 -s ${GTTL}/testdata/sw175.fna >${TMPFILE}
	@diff ${REFDIR}/test_ss_sw175.txt ${TMPFILE}
	@./${SPACED_SEEDS}.x -p -d 1 -s ../testdata/QUERY.fasta >${TMPFILE}
	@diff ${REFDIR}/test_ss_QUERY.txt ${TMPFILE}
	@${RM} ${TMPFILE}
	@echo "$@ passed."

%.o:%.cpp
	$(CXX) -o $@ -c $< ${CXXFLAGS} ${COMPILE_CONST} ${COMPILATION_RESOURCES} ${CFLAGS} ${CPPFLAGS} ${GTTLFLAGS} -MT $@ -MMD -MP -MF $(@:.o=.d)

%.x:%.o
	$(CXX) ${LDFLAGS} $< -o $@ ${LDLIBS}

.PHONY:clean
clean:
	${RM} -r *.[odx] *.x.dSYM/ __pycache__ temp.txt
