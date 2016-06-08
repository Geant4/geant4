#/*+---------------------- Copyright notice -------------------------------+ */
#/*| Copyright (C) 1995, Guy Barrand, LAL Orsay. (barrand@lalcls.in2p3.fr) | */
#/*|   Permission to use, copy, modify, and distribute this software       | */
#/*|   and its documentation for any purpose and without fee is hereby     | */
#/*|   granted, provided that the above copyright notice appear in all     | */
#/*|   copies and that both that copyright notice and this permission      | */
#/*|   notice appear in supporting documentation.  This software is        | */
#/*|   provided "as is" without express or implied warranty.               | */
#/*+---------------------- Copyright notice -------------------------------+ */
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# Examples :
#  UNIX> make OSF1
#  UNIX> make OSF1 target=all
#  UNIX> make OSF1 target=clean
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SHELL = /bin/sh

default :
	@echo "  Type :"
	@echo "     UNIX> make <config>"
	@echo " where, for the distribution, <config> could be :"
	@echo "     OSF1-cxx"
	@echo "     HP-UX-aCC"
	@echo "     Linux-gxx"
	@echo "     SunOS-CC"
	@echo "  For example : "
	@echo "     UNIX> make OSF1"
	@echo "  If none of the above 'config' matches you environment,"
	@echo " edit the Co/v3/mgr/Config.mk file."

OSF1-cxx : This_needed
	@CONFIG=$@;export CONFIG;MAKEFLAGS='';export MAKEFLAGS;$(MAKE) -f This.mk $(t) $(target) \
"a_or_so       = so" \
"so            = so" \
"CC            = cc" \
"CFLAGS        = -std1 -O" \
"CCLD          = cc" \
"CCLDFLAGS     = -std1" \
"CCLDEND       = -lots -lSM -lICE -ldnet_stub" \
"CXX           = cxx" \
"CXXFLAGS      = -std strict_ansi -ptr ${OREPOSITORY} -DG4_HAVE_BOOL -g" \
"CXXLD         = cxx" \
"CXXLDFLAGS    = -ptr ${OREPOSITORY} -non_shared" \
"CXXLDEND      = -lots -lSM -lICE -ldnet_stub" \
"APP_CXXFLAGS  = " \
"F77           = f77" \
"F77FLAGS      = -w" \
"F77LD         = f77" \
"X11_CPPFLAGS  = " \
"Xext_CPPFLAGS = -I/usr/include/X11/extensions" \
"Xmu_CPPFLAGS  = -I/usr/include/X11" \
"Xt_CPPFLAGS   = " \
"Xaw_CPPFLAGS  = " \
"Xm_CPPFLAGS   = " \
"libc          = -lc" \
"libm          = -lm" \
"libf77        = -lfor -lFutil -lUfor" \
"libX11        = -L/usr/shlib -lX11" \
"libXext       = -L/usr/shlib -lXext" \
"libXt         = -L/usr/shlib -lXt" \
"libXmu        = -L/usr/shlib -lXmu" \
"libXaw        = -L/usr/shlib -lXaw" \
"libXm         = -L/usr/shlib -lXm" \
$(This_macros)

HP-UX-aCC : This_needed
	@CONFIG=$@;export CONFIG;MAKEFLAGS='';export MAKEFLAGS;$(MAKE) -f This.mk $(t) $(target) \
"a_or_so       = sl" \
"so            = sl" \
"CC            = c89" \
"CFLAGS        = +Z -O" \
"CCLD          = c89" \
"CCLDFLAGS     = " \
"CCLDEND       = -ldld -lisamstub" \
"CXX           = aCC" \
"CXXFLAGS      = +Z -O -DRW_NO_STL -DG4_HAVE_BOOL -I/opt/aCC/include" \
"CXXLD         = aCC" \
"CXXLDFLAGS    = +Z -O" \
"CXXLDEND      = -lisamstub" \
"APP_CXXFLAGS  = " \
"F77           = f77" \
"F77FLAGS      = -w +ppu" \
"F77LD         = fort77" \
"X11_CPPFLAGS  = " \
"Xext_CPPFLAGS = -I/usr/include/X11/extensions" \
"Xmu_CPPFLAGS  = -I/usr/contrib/X11R6/include/X11" \
"Xt_CPPFLAGS   = " \
"Xaw_CPPFLAGS  = -I/usr/contrib/X11R6/include -I/usr/contrib/X11R6/include/X11" \
"Xm_CPPFLAGS   = " \
"libc          = -lc" \
"libm          = -lm" \
"libf77        = -lf" \
"libX11        = -L/usr/lib -lX11" \
"libXext       = -L/usr/lib -lXext" \
"libXt         = -L/usr/lib -lXt" \
"libXmu        = -L/usr/contrib/X11R6/lib -lXmu" \
"libXaw        = -L/usr/contrib/X11R6/lib -lXaw" \
"libXm         = -L/usr/lib -lXm" \
$(This_macros)

Linux-gxx : This_needed
	@CONFIG=$@;export CONFIG;MAKEFLAGS='';export MAKEFLAGS;$(MAKE) -f This.mk $(t) $(target)  \
"a_or_so       = so" \
"so            = so" \
"CC            = g++" \
"CFLAGS        = -x c++ -Wall -ansi -O" \
"CCCC          = -x c" \
"CCLD          = g++" \
"CCLDFLAGS     = " \
"CCLDEND       = " \
"CXX           = g++" \
"CXXFLAGS      = -pipe -fno-for-scope -DGNU_GCC -DG4_HAVE_BOOL -O" \
"CXXLD         = g++" \
"CXXLDFLAGS    = " \
"CXXLDEND      = " \
"F77           = f77" \
"F77FLAGS      = " \
"F77LD         = f77" \
"X11_CPPFLAGS  = " \
"Xext_CPPFLAGS = -I/usr/include/X11/extensions" \
"Xmu_CPPFLAGS  = -I/usr/include/X11" \
"Xt_CPPFLAGS   = " \
"Xaw_CPPFLAGS  = " \
"Xm_CPPFLAGS   = -I/usr/X11R6/include" \
"libc          = -lc" \
"libm          = -lm" \
"libf77        = -lf" \
"libX11        = -L/usr/X11R6/lib -lX11 -lSM -lICE" \
"libXext       = -L/usr/X11R6/lib -lXext" \
"libXt         = -L/usr/X11R6/lib -lXt" \
"libXmu        = -L/usr/X11R6/lib -lXmu" \
"libXaw        = -L/usr/X11R6/lib -lXaw" \
"libXm         = -L/usr/X11R6/lib -lXm -lXaw" \
$(This_macros)

SunOS-CC : This_needed
	@CONFIG=$@;export CONFIG;MAKEFLAGS='';export MAKEFLAGS;$(MAKE) -f This.mk $(t) $(target)  \
"a_or_so       = a" \
"so            = so" \
"CC            = cc" \
"CFLAGS        = $(opt)" \
"CCLD          = cc" \
"CCLDFLAGS     = " \
"CCLDEND       = " \
"CXX           = CC" \
"CXXFLAGS      = -I/usr/include -ptr${G4WORKDIR}/tmp/$${CONFIG} $(opt)" \
"CXXLD         = CC" \
"CXXLDFLAGS    = -ptr${G4WORKDIR}/tmp/$${CONFIG} $(opt)" \
"CXXLDEND      = -lsocket -lnsl" \
"F77           = f77" \
"F77FLAGS      = " \
"F77LD         = f77" \
"X11_CPPFLAGS  = " \
"Xext_CPPFLAGS = -I/usr/include/X11/extensions" \
"Xmu_CPPFLAGS  = -I/usr/include/X11" \
"Xt_CPPFLAGS   = " \
"Xaw_CPPFLAGS  = " \
"Xm_CPPFLAGS   = " \
"libc          = -lc" \
"libm          = -lm" \
"libf77        = -lF77" \
"libX11        = -L/usr/lib -lX11" \
"libXext       = -L/usr/lib -lXext" \
"libXt         = -L/usr/lib -lXt" \
"libXmu        = -L/usr/lib -lXmu" \
"libXaw        = -L/usr/lib -lXaw" \
"libXm         = -L/usr/lib -lXm" \
$(This_macros)


