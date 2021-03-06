# /**
# *
# * *** Copyright Notice ***
# *
# * PIOK - Parallel I/O Kernels, Copyright (c) 2015-2016, 
# * The Regents of the University of California, through Lawrence
# * Berkeley National Laboratory (subject to receipt of any required approvals
# * from the U.S. Dept. of Energy).  All rights reserved.
# *  
# * If you have questions about your rights to use or distribute this software,
# * please contact Berkeley Lab's Innovation & Partnerships Office at
# * IPO@lbl.gov.
# *  
# * NOTICE.  This Software was developed under funding from the U.S. Department
# * of Energy and the U.S. Government consequently retains certain rights. As
# * such, the U.S. Government has been granted for itself and others acting on
# * its behalf a paid-up, nonexclusive, irrevocable, worldwide license in the
# * Software to reproduce, distribute copies to the public, prepare derivative
# * works, and perform publicly and display publicly, and to permit other to do
# * so.
# *
# */
# *
# * Email questions to SByna@LBL.GOV
# * Scientific Data Management Research Group
# * Lawrence Berkeley National Laboratory
# *
# * Fri Mar 25 09:59:39 PDT 2016
# *
###################################################################
# Change the following directory definition as appropriate
###################################################################
SHELL	= /bin/sh
HOST	= `uname -n`

###################################################################
# C++ compilation definition
###################################################################

TOP	=	/Users/asim/Desktop/flash/src/piok
SRC	= $(TOP)/benchmarks/write/vpicio
#BUILDSUBDIRS = vpicio_hdf5 vpicio_nonuni vpicio_uni vpicio_uni_md
BUILDSUBDIRS = vpicio_hdf5 vpicio_nonuni vpicio_uni
DEPENDSUBDIRS = $(BUILDSUBDIRS)
INSTALLSUBDIRS = $(BUILDSUBDIRS)
###################################################################
# For ALL build for C
###################################################################
all::	$(BUILDSUBDIRS)
	for i in $(BUILDSUBDIRS) ; \
	do \
	(cd $$i ; echo "making" all "in $(TOP)..." ; \
	$(MAKE) all ) ; \
	done

###################################################################
# clean : removes intermediate files
# cleanall : removes intermediate files and executables
###################################################################
clean:: $(BUILDSUBDIRS)
	for i in $(BUILDSUBDIRS) ; \
	do \
	(cd $$i ; echo "making" clean "in $(TOP)/$$i..." ; \
	$(MAKE) clean ) ; \
	done

cleanall:: $(BUILDSUBDIRS)
	for i in $(BUILDSUBDIRS) ; \
	do \
	(cd $$i ; echo "making" cleanall "in $(TOP)/$$i..." ; \
	$(MAKE) cleanall ) ; \
	done

###################################################################
# For install
###################################################################
install::	$(INSTALLSUBDIRS)
	for i in $(INSTALLSUBDIRS) ; \
	do \
	(cd $$i ; echo "making" install "in $(TOP)/$$i..." ; \
	${MAKE} install) ; \
	done

######################
depend:	$(DEPENDSUBDIRS)
	for i in $(DEPENDSUBDIRS) ; \
	do \
	(cd $$i ; echo building Makefile "in $$i..." ; \
	$(MAKE) depend $(BUILDFLAG) ) ; \
	done

###################################################################
# to make a Makefile, make -f thisfile Makefile FRC=rebuild
# for dependency generation, make sure -xM is for CC only
#                       and use -MM or -M for g++
###################################################################
Makefile	: $(FRC)
	@ \rm -f $@
	@ \cp Makefile.in $@
	@ chmod 440 $@
	@ for i in $(BUILDSUBDIRS) ; \
	do \
	(cd $$i ; echo  building Makefile "in $$i..." ; \
	$(MAKE) -f Makefile.in Makefile FRC=rebuild) ; \
	done

###################################################################
# for forceful rebuild
# automatically generates dependencies of each source files
###################################################################
rebuild	:
