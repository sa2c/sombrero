TOPDIR = .
SUBDIRS := LibHR sombrero

MKDIR = $(TOPDIR)/Make

# Forces compiling one at a time even with -j
default:
	$(MAKE) sombrero/sombrero1
	$(MAKE) sombrero/sombrero2
	$(MAKE) sombrero/sombrero3
	$(MAKE) sombrero/sombrero4
	$(MAKE) sombrero/sombrero5
	$(MAKE) sombrero/sombrero6


include $(MKDIR)/MkRules

