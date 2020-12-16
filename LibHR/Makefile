TOPDIR = ..
MKDIR = $(TOPDIR)/Make

SUBDIRS = $(shell find * -maxdepth 0 -type d)
LIBDIRS = $(SUBDIRS)

clean::
	rm -fr su2 su2adj su3 su3sym sp4 sp4adj

include $(MKDIR)/MkRules

