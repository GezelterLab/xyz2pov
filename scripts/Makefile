#####################################################
# Makefile for xyz2pov Scripts                      #
#####################################################
# You should not change anything here.              #
#####################################################

TOP=..

include $(TOP)/Make.conf

SCRIPT_DIR=$(XYZ2POV_HOME)/bin

SCRIPT_FILES= \
	renderAll\
	renderRange 

all: dummy

install : 
	$(MKINSTALLDIRS) $(SCRIPT_DIR)
	for i in $(SCRIPT_FILES); do \
	  echo "Installing $$i in $(SCRIPT_DIR)\n"; \
	  $(INSTALL_PROGRAM) $$i $(SCRIPT_DIR)/$$i; \
	done

tests : dummy

clean : dummy
	$(RM) *~

distclean : dummy
	$(RM) *~

devclean : distclean

depend : dummy

links:
	$(MKINSTALLDIRS) $(PREFIX)/bin
	for i in $(SCRIPT_FILES); do \
	  cd $(PREFIX)/bin; $(LN_S)  -f $(SCRIPT_DIR)/$$i $$i; \
	done

tags : dummy

dummy :
