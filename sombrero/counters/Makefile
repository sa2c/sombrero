
TOPDIR=../../

SRCS=flop_count_template.py.c memory_count_template.py.c transfer_count_template.py.c
TARGETS=(patsubst %_template.py.c,%.o,$(SRCS))

$(TARGETS) :

%.o: %.c $(TOPDIR)/Include/global.h
	gcc -o $@ -I $(TOPDIR)/Include -c $<
