#####################################################
# Makefile for xyz2pov AtomType files               #
#####################################################
# You should not change anything here.              #
#####################################################

TOP=..

include $(TOP)/Make.conf

TYPE_FILES= \
	AtomTypes

all: dummy

install : 
	$(MKINSTALLDIRS) $(ATOMTYPES_DIR)
	for i in $(TYPE_FILES); do \
	  echo "Installing $$i in $(ATOMTYPES_DIR)"; \
	  $(INSTALL_DATA) $$i $(ATOMTYPES_DIR)/$$i; \
	done

tests : dummy

clean : dummy
	$(RM) *~

distclean : dummy
	$(RM) *~

devclean : distclean

depend : dummy

links: dummy

tags : dummy

dummy :
