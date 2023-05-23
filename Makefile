all:
	@make -C testsuite test

debug:
	@make -C testsuite debug=yes test

.PHONY:clean	
clean:
	@make -C testsuite clean

.PHONY:tags
tags:
	@${GTTL}/scripts/ctags.sh

.PHONY:code_check
code_check:
	@scripts/code_check.py -wt `find . -name '*.[ch]pp'`
	@${GTTL}/scripts/check_ifndef.py `find . -name '*.hpp'`
