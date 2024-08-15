# $Id: Makefile.rdkit_align3d_test.app 673324 2023-09-28 13:21:53Z thiessen $

REQUIRES = RDKIT MT

APP = rdkit_align3d_test

SRC = rdkit_align3d_test

CPPFLAGS = $(ORIG_CPPFLAGS) $(RDKIT_INCLUDE)

LIB = align3d

LIBS = $(ORIG_LIBS) $(RDKIT_LIBS) 
