###############################################################################
## Makefile for ULySS
## ULySS is developed and maintained by Mina Koleva and collaborators
##            tel. (+33) 04 78 86 83 83 fax. (+33) 04 72 72 87 87
##            mina.koleva@obs.univ-lyon1.fr
## Copyright (C) 2006-2011; 
###############################################################################
## Installation procedure:
##   In short, type:   make install
##   Play the file: pgm/uly_startup.pro in your GDL/IDL startup program
###############################################################################

.DEFAULT help:  
	@ echo "Usage:"; 
	@ echo "   make install"; 
	@ echo "      Test the environment (availability of dependences)"; 
	@ echo "      and generate the file: pgm/uly_startup.pro"; 
	@ echo "   make test"; 
	@ echo "      Execute the test program uly_test.pro"; 
	@ echo "   make update"; 
	@ echo "      Update the package with respect to the cvs repository"; 
	@ echo "      (a login/password is required)"; 
	@ echo "   make tar"; 
	@ echo "      Prepare ulyss.tar to release the package"; 
	@ echo "------"; 

install: pgm/uly_startup.pro
pgm/uly_startup.pro: pgm/uly_startup.pro.in configure
	./configure

configure: configure.in; autoconf  # create/update 'configure'

test:
	@ cd pgm && echo "uly_test & exit" | idl

#### The following targets are for maintenance and development ################

up update:
	@ cvs up

CTG1="; CATEGORY:\s*ULY\s*$$"
CTG2="; CATEGORY:\s*ULY_SSP\s*$$"
CTG3="; CATEGORY:\s*ULY_TGM\s*$$"

docu doc: pgm/*.pro
# Build the cross reference index
	@rm -f xref.txt; 

	@a="[" ;\
	for file in `grep -l ${CTG1} pgm/*.pro`; do \
	  a=$$a\'$$file\',; \
	done; \
	a=`echo $$a | sed s/.$$/]/`; \
	echo "uly_HTML_HELP, $$a, 'doc/uly.html', /CROSS, XREF='xref.txt'& exit" | idl > /dev/null && echo doc ctg1 OK

	@a="[" ;\
	for file in `grep -l ${CTG2} pgm/*.pro`; do \
	  a=$$a\'$$file\',; \
	done; \
	a=`echo $$a | sed s/.$$/]/`; \
	echo "uly_HTML_HELP, $$a, 'doc/ssp.html', /CROSS, XREF='xref.txt' & exit" | idl > /dev/null && echo doc ctg2 OK

	@a="[" ;\
	for file in `grep -l ${CTG3} pgm/*.pro`; do \
	  a=$$a\'$$file\',; \
	done; \
	a=`echo $$a | sed s/.$$/]/`; \
	echo "uly_HTML_HELP, $$a, 'doc/tgm.html', /CROSS, XREF='xref.txt'& exit" | idl > /dev/null && echo doc ctg3 OK

	@a="[" ;\
	for file in `grep -L -e ${CTG1} -e ${CTG2} -e ${CTG3} pgm/*.pro`; do \
	  a=$$a\'$$file\',; \
	done; \
	a=`echo $$a | sed s/.$$/]/`; \
	echo "uly_HTML_HELP, $$a, 'doc/util.html', /CROSS, XREF='xref.txt'& exit" | idl > /dev/null && echo doc util OK

	@echo "Cross reference index OK"

# Build the html files
	@a="[" ;\
	for file in `grep -l ${CTG1} pgm/*.pro`; do \
	  a=$$a\'$$file\',; \
	done; \
	a=`echo $$a | sed s/.$$/]/`; \
	echo "uly_HTML_HELP, $$a, 'doc/uly.html', TIT='ULySS main documentation', STY='ulyss.css', HEAD='<!--#include virtual=\"header.html\" -->', XREF='xref.txt' & exit" | idl > /dev/null && echo doc ctg1 OK

	@a="[" ;\
	for file in `grep -l ${CTG2} pgm/*.pro`; do \
	  a=$$a\'$$file\',; \
	done; \
	a=`echo $$a | sed s/.$$/]/`; \
	echo "uly_HTML_HELP, $$a, 'doc/ssp.html', TIT='ULySS SSP documentation', STY='ulyss.css', HEAD='<!--#include virtual=\"header.html\" -->', XREF='xref.txt' & exit" | idl > /dev/null && echo doc ctg2 OK

	@a="[" ;\
	for file in `grep -l ${CTG3} pgm/*.pro`; do \
	  a=$$a\'$$file\',; \
	done; \
	a=`echo $$a | sed s/.$$/]/`; \
	echo "uly_HTML_HELP, $$a, 'doc/tgm.html', TIT='ULySS TGM documentation', STY='ulyss.css', HEAD='<!--#include virtual=\"header.html\" -->', XREF='xref.txt' & exit" | idl > /dev/null && echo doc ctg3 OK

	@a="[" ;\
	for file in `grep -L -e ${CTG1} -e ${CTG2} -e ${CTG3} pgm/*.pro`; do \
	  a=$$a\'$$file\',; \
	done; \
	a=`echo $$a | sed s/.$$/]/`; \
	echo "uly_HTML_HELP, $$a, 'doc/util.html', TIT='ULySS utilities documentation', STY='ulyss.css', HEAD='<!--#include virtual=\"header.html\" -->', XREF='xref.txt'& exit" | idl > /dev/null && echo doc util OK
	@ touch doc

tut: 
	@for file in `ls tuto/*.pro`; do \
	  echo "pro2html, SIG='', /INCL, XREF='xref.txt', PROM='GDL>', '$$file', 'doc'" | idl -arg quiet> /dev/null; \
	done; \
	echo "pro2html, SIG='', /INCL, XREF='xref.txt', 'tuto/tuto_power.pro', 'doc'" | idl -arg quiet > /dev/null

	@ echo "@tuto/tuto" | idl  -arg quiet> /dev/null


tar: optional clean
	@ test -f models/PHR_Elodie31.fits -a -f models/elodie32_flux_tgm.fits -a -f models/sun.fits -a -f models/Vaz_Miles.fits || { echo 'miss some model fils'; exit 1; }
	@ cd .. && tar --exclude CVS  --exclude local --exclude TODO --exclude \*.tar.gz --exclude \*~ --exclude \*.res --exclude output.\* --exclude \*.cache --exclude xref.txt --exclude \*.log  --exclude pgm/uly_startup.pro  --exclude .\#\* --exclude doc --exclude store -czf ulyss.tar.gz ulyss/

clean:	
	rm -f *~ pgm/*~ doc/*~ dep/*~ tuto/*~ doc/optional/*~ doc/optional/aux_*~

# Build the index of optional model files and web page with the auxillary routines
optional:doc/optional/index.html aux
doc/optional/index.html: doc/optional/*fits* doc/optional/makefile doc/optional/aux_*.pro
	cd doc/optional; make index.html auxhtml

### End of ULySS makefile #####################################################
