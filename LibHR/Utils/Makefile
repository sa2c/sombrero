TOPDIR = ../..
MKDIR = $(TOPDIR)/Make

SRCS = $(wildcard *.c)
LIBNAME = ../libhr.a

#Exclude this source files
EXCLUDE += $(if $(findstring 2,$(NG)),,HYP_smearing.c)

include $(MKDIR)/MkRules

clean::

