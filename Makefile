#####################################################
# Top-level Makefile for Grace                      #
#####################################################
# You should not change anything here.              #
#####################################################

include Make.conf

subdirs : configure Make.conf
	@set -e; for i in $(SUBDIRS); do (cd $$i; $(MAKE)) || exit 1; done

all : subdirs

install : subdirs
	@set -e; for i in $(SUBDIRS); do (cd $$i; $(MAKE) install) || exit 1; done
#	$(MKINSTALLDIRS) $(XYZ2POV_HOME)
#	$(INSTALL_DATA) README $(XYZ2POV_HOME)
#	$(INSTALL_DATA) LICENSE $(XYZ2POV_HOME)
#	$(INSTALL_DATA) NEWS $(XYZ2POV_HOME)
#	$(INSTALL_DATA) AUTHORS $(XYZ2POV_HOME)
#	$(INSTALL_DATA) ChangeLog $(XYZ2POV_HOME)

tests : subdirs
	@set -e; for i in $(SUBDIRS); do (cd $$i; $(MAKE) tests) || exit 1; done

check : tests

links : subdirs
	@set -e; for i in $(SUBDIRS); do (cd $$i; $(MAKE) links) || exit 1; done

clean :
	@set -e; for i in $(SUBDIRS); do (cd $$i; $(MAKE) clean) || exit 1; done

distclean :
	@set -e; for i in $(SUBDIRS); do (cd $$i; $(MAKE) distclean) || exit 1; done
	$(RM) config.log config.status config.cache Make.conf

devclean :
	@set -e; for i in $(SUBDIRS); do (cd $$i; $(MAKE) devclean) || exit 1; done
	$(RM) config.log config.status config.cache Make.conf configure

texts : CHANGES ChangeLog

CHANGES : doc/CHANGES.html
	@lynx -dump $? > CHANGES

ChangeLog : 
	./scripts/cvs2cl.pl

Make.conf : ac-tools/Make.conf.in configure
	@echo
	@echo 'Please re-run ./configure'
	@echo
	@exit 1

configure : ac-tools/configure.in ac-tools/aclocal.m4
	autoconf ac-tools/configure.in > configure
	chmod +x configure

dummy :

