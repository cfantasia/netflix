# Makefile for the ROOT test programs.
# This Makefile shows nicely how to compile and link applications
# using the ROOT libraries on all supported platforms.
#
# Copyright (c) 2000 Rene Brun and Fons Rademakers
#
# Author: Fons Rademakers, 29/2/2000

include Makefile.arch

#------------------------------------------------------------------------------
STORAGEO       = Storage_dict.$(ObjSuf)
STORAGES       = Storage_dict.$(SrcSuf) 
STORAGE        = Storage_dict$(ExeSuf)

NETFLIXO	= netflix.$(ObjSuf)
NETFLIXS       = netflix.$(SrcSuf)
NETFLIX        = netflix$(ExeSuf)

OBJS	      = $(STORAGEO) $(NETFLIXO) 
#OBJS          = $(EVENTO) $(MAINEVENTO) $(HWORLDO) $(HSIMPLEO) $(MINEXAMO) \

PROGRAMS      = $(NETFLIX) 
#PROGRAMS      = $(EVENT) $(HWORLD) $(HSIMPLE) $(MINEXAM) $(TSTRING) \

#OBJS         += $(GUITESTO) $(GUIVIEWERO) $(TETRISO)
#PROGRAMS     += $(GUITEST) $(GUIVIEWER) $(TETRISSO)

#------------------------------------------------------------------------------

.SUFFIXES: .$(SrcSuf) .$(ObjSuf) .$(DllSuf)
.PHONY:    Aclock Hello Tetris $(STORAGE)

all:		$(PROGRAMS)
		@echo "$@ done"

$(NETFLIX):	$(STORAGEO) $(NETFLIXO)
		$(LD) $(LDFLAGS) $^ $(LIBS) $(OutPutOpt)$@  
		@echo "$@ done"

$(STORAGE):     
		rootcint -f Storage_dict.cxx -c $(CORY_INC) Storage.h Storage.cxx Storage_LinkDef.h
		@echo "$@ done"

clean:
		@rm -f $(OBJS) core

distclean:      clean
		@rm -f $(PROGRAMS) $(EVENTSO) $(EVENTLIB) *Dict.* *.def *.exp \
		   *.root *.ps *.so *.lib *.dll *.d .def so_locations
		@rm -rf cxx_repository
		-@cd RootShower && $(MAKE) distclean

.SUFFIXES: .$(SrcSuf)


###

CORY_INC = -I/home/fantasia/work/netflix/include

.$(SrcSuf).$(ObjSuf):
	$(CXX) $(CXXFLAGS) $(CORY_INC) -c $<
