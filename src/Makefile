#####################################################
# Makefile for xyz2pov                              #
#####################################################
# You should not change anything here.              #
#####################################################

TOP=..

include $(TOP)/Make.conf

.SUFFIXES : .cpp .c .h .hpp $(O)

# First do compile-time flags:

DECLARES= -DTYPES_PATH="$(ATOMTYPES_DIR)"

# include paths starting with "." are all relative to the SUBDIRS!!!

INCLUDES= -I.. 
CXXFLAGS=$(CXXFLAGS_C) $(INCLUDES) $(DECLARES)
CFLAGS=$(CFLAGS_C) $(INCLUDES) $(DECLARES)
DEPFLAGS=$(CFLAGS_C) -I.

XYZ2POV=xyz2pov$(EXE)

POVTOOLS_LIB = libpovtools.a

XYZ2POV_LDFLAGS= -L.
XYZ2POV_LDADD= -lpovtools -lm

PREFS=-DXYZ2POV_HOME="$(XYZ2POV_HOME)"

# Then do the compile rules:

SRCPATH=$(PWD)

obj/%.o: %.cpp
	cd obj; $(CXX) $(CXXFLAGS) -c $(SRCPATH)/$<

obj/%.o: %.c
	cd obj; $(CC) $(CFLAGS) -c $(SRCPATH)/$<

# Then list what we have to compile

LIBC_FILES = \
	$(DIR)/atom_parser.c \
	$(DIR)/frameCount.c \
	$(DIR)/pov_writer.c

LIBCXX_FILES = 

XYZC_FILES = \
	$(DIR)/xyz2pov.c

XYZCXX_FILES = 

DIR:=.
LIBC_SRCS   :=   $(LIBC_FILES)
LIBCXX_SRCS :=	 $(LIBCXX_FILES)
LIB_SRCS    :=   $(LIBC_SRCS) $(LIBCXX_SRCS)

DIR:=.
XYZC_SRCS   :=  $(XYZC_FILES)
XYZCXX_SRCS :=	$(XYZCXX_FILES)
XYZ_SRCS    :=  $(XYZC_SRCS) $(XYZCXX_SRCS)

SRCS := $(LIB_SRCS) $(XYZ_SRCS)

DIR:=obj
LIBC_OBJS   :=     $(LIBC_FILES:.c=.o)
LIBCXX_OBJS :=	   $(LIBCXX_FILES:.cpp=.o)
LIB_OBJS    :=     $(LIBC_OBJS) $(LIBCXX_OBJS)

DIR:=obj
XYZC_OBJS   :=  $(XYZC_FILES:.c=.o)
XYZCXX_OBJS :=	$(XYZCXX_FILES:.cpp=.o)
XYZ_OBJS    :=  $(XYZC_OBJS) $(XYZCXX_OBJS)

LIBSUBDIRS = obj

all : Make.dep $(POVTOOLS_LIB) $(XYZ2POV)

depend: Make.dep

libpovtools.a: $(LIB_OBJS)
	$(RM) $@
	$(AR) cr $@ $(LIB_OBJS)
	$(RANLIB) $@

xyz2pov$(EXE) : $(XYZ_OBJS) $(POVTOOLS_LIB)
	$(CC) $(CFLAGS) $(XYZ_OBJS) -o $@ $(XYZ2POV_LDFLAGS) $(XYZ2POV_LDADD)

install : $(XYZ2POV)
	$(MKINSTALLDIRS) $(XYZ2POV_HOME)/bin
	$(INSTALL_PROGRAM) -s $(XYZ2POV) $(XYZ2POV_HOME)/bin/$(XYZ2POV)

tests : dummy

$(XYZ2POV_HOME)/bin/$(XYZ2POV) : $(XYZ2POV)
	$(MAKE) install

links : $(XYZ2POV_HOME)/bin/$(XYZ2POV)
	$(MKINSTALLDIRS) $(PREFIX)/bin
	cd $(PREFIX)/bin; $(RM) $(XYZ2POV); $(LN_S) $(XYZ2POV_HOME)/bin/$(XYZ2POV) $(XYZ2POV)

# cleaning

.PHONY : clean

clean : dummy
	for i in $(LIBSUBDIRS); do \
	  (cd $$i; $(RM) *.o) || exit 1; \
	done
	$(RM) $(POVTOOLS_LIB) $(XYZ2POV) *.o *~ Make.temp

distclean : dummy
	for i in $(LIBSUBDIRS); do \
	  (cd $$i; $(RM) *.o) || exit 1; \
	done
	$(RM) $(POVTOOLS_LIB) $(XYZ2POV) *.o *~ Make.temp Make.dep

devclean : distclean

Make.dep : $(SRCS)
	echo "# DO NOT DELETE THIS LINE - used by make depend" > Make.dep
	$(CC) $(DEPFLAGS) -MM $(SRCS) \
		| sed 's/\.o:/\$$(O)\ :/g' > Make.temp
	cat Make.temp | sed 's/^[0-9a-zA-Z]/obj\/&/' >> Make.dep
	$(RM) Make.temp

tags : $(SRCS)
	ctags $(SRCS)

dummy :

# DO NOT DELETE
include Make.dep
