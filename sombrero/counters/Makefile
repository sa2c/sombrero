
TOPDIR=../..

SRCS=flop_count_template.py.c memory_count_template.py.c communication_count_template.py.c
TARGETS=$(patsubst %_template.py.c,lib%.o,$(SRCS))

all : $(TARGETS)

#$(TARGETS) :

lib%.o: %_template.py.c $(TOPDIR)/Include/global.h
	echo Compiling $@ && $(CC) $(CFLAGS) -I $(TOPDIR)/Include -c $< -o $@

clean:
	rm -f $(TARGETS)
