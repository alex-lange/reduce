#
# Created by gmakemake (Ubuntu Mar 10 2012) on Thu Apr  4 23:05:48 2013
#

#
# Definitions
#

.SUFFIXES:
.SUFFIXES:	.a .o .c .C .cpp .s .S
.c.o:
		$(COMPILE.c) $<
.C.o:
		$(COMPILE.cc) $<
.cpp.o:
		$(COMPILE.cc) $<
.S.s:
		$(CPP) -o $*.s $<
.s.o:
		$(COMPILE.s) -o $@ $<
.c.a:
		$(COMPILE.c) -o $% $<
		$(AR) $(ARFLAGS) $@ $%
		$(RM) $%
.C.a:
		$(COMPILE.cc) -o $% $<
		$(AR) $(ARFLAGS) $@ $%
		$(RM) $%
.cpp.a:
		$(COMPILE.cc) -o $% $<
		$(AR) $(ARFLAGS) $@ $%
		$(RM) $%

AS =		as
CC =		gcc
CXX =		g++

RM = rm -f
AR = ar
LINK.c = $(CC) $(CFLAGS) $(CPPFLAGS) $(LDFLAGS)
LINK.cc = $(CXX) $(CXXFLAGS) $(CPPFLAGS) $(LDFLAGS)
COMPILE.s = $(AS) $(ASFLAGS)
COMPILE.c = $(CC) $(CFLAGS) $(CPPFLAGS) -c
COMPILE.cc = $(CXX) $(CXXFLAGS) $(CPPFLAGS) -c
CPP = $(CPP) $(CPPFLAGS)
########## Default flags (redefine these with a header.mak file if desired)
CXXFLAGS =	-ggdb -std=c++11 -I$(FLENSFILES)
CFLAGS =	-ggdb
CLIBFLAGS =	-lm
CCLIBFLAGS =	
########## End of default flags


g_path = /home/alex/research/archer

CPP_FILES =	reduce.cpp reduce_ramsey.cpp test.cpp lll.cpp gen_gr.cpp
C_FILES =	
PS_FILES =	
S_FILES =	
H_FILES =	lll.h
SOURCEFILES =	$(H_FILES) $(CPP_FILES) $(C_FILES) $(S_FILES)
.PRECIOUS:	$(SOURCEFILES)
OBJFILES =	lll.o
GOBJFILES = 	$(g_path)/g.o $(g_path)/set_operations.o $(g_path)/number_algorithms.o $(g_path)/g_core.o

FLENSFILES = /home/alex/software/flens/FLENS/

#
# Main targets
#

all:	reduce reduce_ramsey test


reduce:	reduce.o $(OBJFILES) $(GOBJFILES)
	$(CXX) $(CXXFLAGS) -o reduce reduce.o $(OBJFILES) $(GOBJFILES) $(CCLIBFLAGS)

reduce_ramsey:	reduce_ramsey.o $(OBJFILES) $(GOBJFILES)
	$(CXX) $(CXXFLAGS) -o reduce_ramsey reduce_ramsey.o $(OBJFILES) $(GOBJFILES) $(CCLIBFLAGS)

test:	test.o $(OBJFILES) $(GOBJFILES)
	$(CXX) $(CXXFLAGS) -o test test.o $(OBJFILES) $(GOBJFILES) $(CCLIBFLAGS)

gen_gr:	gen_gr.o $(OBJFILES) $(GOBJFILES)
	$(CXX) $(CXXFLAGS) -o gen_gr gen_gr.o $(OBJFILES) $(GOBJFILES) $(CCLIBFLAGS)


#
# Dependencies
#

reduce.o:

#
# Housekeeping
#

Archive:	archive.tgz

archive.tgz:	$(SOURCEFILES) Makefile
	tar cf - $(SOURCEFILES) Makefile | gzip > archive.tgz

clean:
	-/bin/rm -f $(OBJFILES) reduce.o reduce_ramsey.o test.o core

realclean:        clean
	-/bin/rm -f reduce 
