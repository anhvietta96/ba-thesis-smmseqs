#!/bin/sh

checkerror() {
if test $? -ne 0
then
  echo "failure: ${cmd}"
  exit 1
fi
}

VALGRIND=valgrind

if test ! -f ${VALGRIND}
then
  VALGRIND=valgrind
fi

valgrind --help > /dev/null
if test $? -eq 127
then
  echo "valgrind not available"
  exit 0
fi

# add --gen-suppressions=all to generate suppression format and
# add it to Admin/sk.supp

cmd="${VALGRIND} --quiet --tool=memcheck --memcheck:leak-check=full --memcheck:leak-resolution=high --error-exitcode=1 --log-fd=1 --error-limit=yes --gen-suppressions=all --dsymutil=yes $*"
${cmd}
checkerror

# for filename in `ls Valgrind-error.*`
# do
  # filesize=`cat ${filename} | wc -c`
  # if test ${filesize} -eq 0
  # then
    # echo "remove ${filename}"
    # rm -f ${filename}
  # else
    # echo "there seems to be an error"
  # fi
# done
