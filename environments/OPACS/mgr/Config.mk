#/* +---------------------- Copyright notice -------------------------------+ */
#/* | Copyright (C) 1995, Guy Barrand, LAL Orsay. (barrand@lalcls.in2p3.fr) | */
#/* |   Permission to use, copy, modify, and distribute this software       | */
#/* |   and its documentation for any purpose and without fee is hereby     | */
#/* |   granted, provided that the above copyright notice appear in all     | */
#/* |   copies and that both that copyright notice and this permission      | */
#/* |   notice appear in supporting documentation.  This software is        | */
#/* |   provided "as is" without express or implied warranty.               | */
#/* +---------------------- Copyright notice -------------------------------+ */
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
	@echo "     OSF1"
	@echo "     HP-UX"
	@echo "     HP-UX-aCC"
	@echo "     Linux"
	@echo "     ULTRIX"
	@echo "     AIX"
	@echo "     IRIX"
	@echo "     SunOS"
	@echo "     Linux-gxx = Linux compiling with g++."
	@echo "     insure = OSF1 + insure product."
	@echo "     HP-UX-9 = HP-UX A.09.05."
	@echo "  For example : "
	@echo "     UNIX> make OSF1"
	@echo "  If none of the above 'config' matches you environment,"
	@echo " edit the Co/v3/mgr/Config.mk file."

# Perfect situation still never reached !!!
POSIX : This_needed
	@CONFIG=$@;export CONFIG;MAKEFLAGS='';export MAKEFLAGS;$(MAKE) -f This.mk $(t) $(target) \
"a_or_so       = so" \
"so            = so" \
"CC            = cc" \
"CFLAGS        = -DUNIX -DHAS_DLD -O" \
"CCLD          = cc" \
"CCLDFLAGS     = " \
"CCLDEND       = " \
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
"libf77        = -lf" \
"libX11        = -L/usr/lib -lX11" \
"libXext       = -L/usr/lib -lXext" \
"libXt         = -L/usr/lib -lXt" \
"libXmu        = -L/usr/lib -lXmu" \
"libXaw        = -L/usr/lib -lXaw" \
"libXm         = -L/usr/lib -lXm" \
$(This_macros)


# OSF1 asa V4.0 464 alpha
# CCLDFLAGS : -non_shared to compell link with .a libs
# For cxx 5.x : APP_CXXFLAGS  = -define_templates -DSOLVE_TEMPLATES
OSF1 : This_needed
	@CONFIG=$@;export CONFIG;MAKEFLAGS='';export MAKEFLAGS;$(MAKE) -f This.mk $(t) $(target) \
"a_or_so       = so" \
"so            = so" \
"CC            = cc" \
"CFLAGS        = -std1 -DUNIX -DHAS_DLD -O" \
"CCLD          = cc" \
"CCLDFLAGS     = -std1" \
"CCLDEND       = -lots -lSM -lICE -ldnet_stub" \
"CXX           = cxx" \
"CXXFLAGS      = -ptr ${OREPOSITORY} -I${RWROOT}/include -g" \
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
"JAVA_CPPFLAGS = -D_OSF1_SOURCE -I/usr/local/java/include -I/usr/local/java/include/osfport/" \
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

# HP-UX aleph B.10.20 A 9000/780 2011153535 two-user license
#"CCLDFLAGS   = -Wl,-a,archive"
#"CCLDFLAGS   = -Wl,-a,archive_shared"
# At least on 10.20, Xm and Xaw are incompatibles.
# At load, put "archive_shared" so that linking G4o apps with G4 libs works.
#"libXm = -L/usr/lib -lXm" is imcompatible with loading option "-Wl,-a,archive_shared" needed for linking with G4.
#"libXm = /usr/lib/Motif1.2_R6/libXm.a"
HP-UX : This_needed
	@CONFIG=$@;export CONFIG;MAKEFLAGS='';export MAKEFLAGS;$(MAKE) -f This.mk $(t) $(target) \
"a_or_so       = sl" \
"so            = sl" \
"CC            = c89" \
"CFLAGS        = +Z -DUNIX -DHAS_DLD -O" \
"CCLD          = c89" \
"CCLDFLAGS     = " \
"CCLDEND       = -ldld -lisamstub" \
"CXX           = CC" \
"CXXFLAGS      = +a1 +Z -pta -ptr${OREPOSITORY} -I${RWROOT}/include -O" \
"CXXLD         = CC" \
"CXXLDFLAGS    = +a1 +Z -pta -ptr${OREPOSITORY} -O" \
"CXXLDEND      = -lisamstub" \
"APP_CXXFLAGS  = -DSOLVE_TEMPLATES" \
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

# HP-UX aleph B.10.20 A 9000/780 2011153535 two-user license
HP-UX-aCC : This_needed
	@CONFIG=$@;export CONFIG;MAKEFLAGS='';export MAKEFLAGS;$(MAKE) -f This.mk $(t) $(target) \
"a_or_so       = sl" \
"so            = sl" \
"CC            = c89" \
"CFLAGS        = +Z -DUNIX -DHAS_DLD -O" \
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

# HP-UX cchp004 A.09.05 A 9000/770 2004200392 two-user license
#"CCLDFLAGS   = -Wl,-a,archive"
#"CCLDFLAGS   = -Wl,-a,archive_shared"
# At load, put "archive_shared" so that linking G4o apps with G4 libs works.
HP-UX-9 : This_needed
	@CONFIG=$@;export CONFIG;MAKEFLAGS='';export MAKEFLAGS;$(MAKE) -f This.mk $(t) $(target) \
"a_or_so       = sl" \
"so            = sl" \
"CC            = c89" \
"CFLAGS        = +Z -D_HPUX_SOURCE -DUNIX -DHAS_DLD -O" \
"CCLD          = c89" \
"CCLDFLAGS     = " \
"CCLDEND       = -ldld -lisamstub" \
"CXX           = CC" \
"CXXFLAGS      = +a1 +Z -pts -ptr${OREPOSITORY} -O" \
"CXXLD         = CC" \
"CXXLDFLAGS    = +a1 +Z -pts -ptr${OREPOSITORY} -O -Wl,-a,archive_shared" \
"CXXLDEND      = " \
"APP_CXXFLAGS  = -DSOLVE_TEMPLATES" \
"F77           = f77" \
"F77FLAGS      = -w +ppu" \
"F77LD         = fort77" \
"X11_CPPFLAGS  = -I/usr/include/X11R5" \
"Xext_CPPFLAGS = -I/usr/include/X11R5/X11/extensions" \
"Xmu_CPPFLAGS  = -I../usr/HP-UX/include/X11" \
"Xt_CPPFLAGS   = " \
"Xaw_CPPFLAGS  = -I../usr/HP-UX/include" \
"Xm_CPPFLAGS   = -I/usr/include/Motif1.2" \
"libc          = -lc" \
"libm          = -lm" \
"libf77        = -lf" \
"libX11        = -L/usr/lib/X11R5 -lX11" \
"libXext       = -L/usr/lib/X11R5 -lXext" \
"libXt         = -L/usr/lib/X11R5 -lXt" \
"libXmu        = /usr/lib/X11R4/libXmu.sl" \
"libXaw        = /usr/lib/X11R4/libXaw.sl" \
"libXm         = -L/usr/lib/Motif1.2 -lXm" \
$(This_macros)

#ULTRIX ux3 4.4 0 RISC
ULTRIX : This_needed
	@CONFIG=$@;export CONFIG;MAKEFLAGS='';export MAKEFLAGS;$(MAKE) -f This.mk $(t) $(target) \
"a_or_so       = a" \
"so            = so" \
"CC            = cc" \
"CFLAGS        = -YPOSIX -DUNIX" \
"CCLD          = cc" \
"CCLDFLAGS     = " \
"CCLDEND       = " \
"F77           = f77" \
"F77FLAGS      = -w -G 3 -static -u" \
"F77LD         = f77" \
"X11_CPPFLAGS  = " \
"Xext_CPPFLAGS = -I/usr/include/X11/extensions" \
"Xmu_CPPFLAGS  = -I/usr/include/X11" \
"Xt_CPPFLAGS   = " \
"Xaw_CPPFLAGS  = " \
"Xm_CPPFLAGS   = " \
"libc          = -lc" \
"libm          = -lm" \
"libf77        = -lfor -lutil -lUfor -lots -li" \
"libX11        = -lX11" \
"libXext       = -lXext" \
"libXt         = -lXt" \
"libXmu        = -lXmu" \
"libXaw        = -lXaw" \
"libXm         = -lXm" \
$(This_macros)

# "OSF1 asa V4.0 464 alpha" + insight
# Do not use shared lib to have source code visibility.
# In CXXFLAGS "-ptr ${OREPOSITORY}" should be put in a .insight file with : preprocessor_flag -ptr ${OREPOSITORY}
insure : This_needed
	@CONFIG=$@;export CONFIG;MAKEFLAGS='';export MAKEFLAGS;$(MAKE) -f This.mk $(t) $(target) \
"a_or_so       = so" \
"so            = so" \
"CC            = insight" \
"CFLAGS        = -std1 -DUNIX -DHAS_DLD -g" \
"CCLD          = insight" \
"CCLDFLAGS     = -std1 -g" \
"CCLDEND       = -lots -L/lal/insure/3.0.1/lib.alpha  -ltqsiiccc" \
"CXX           = insight" \
"CXXFLAGS      = -ptr ${OREPOSITORY} -g" \
"CXXLD         = insight" \
"CXXLDFLAGS    = -ptr ${OREPOSITORY} -non_shared" \
"CXXLDEND      = -lcxx -lexc -lots -lSM -lICE -ldnet_stub -L/lal/insure/3.0.1/lib.alpha  -ltqsiiccc" \
"F77           = f77" \
"F77FLAGS      = -w -g" \
"F77LD         = f77" \
"F77LDEND      = -L/lal/insure/3.0.1/lib.alpha  -linsight" \
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

# Linux pcbp.lal.in2p3.fr 2.0.25 #2 Thu Dec 5 09:56:01 MET 1996 i586
# CXXFLAGS dedicated to compile G4 include files.
#"CCLDFLAGS     = -Wl,-Bstatic"
# Xm from XiG (previous Xinside).
Linux : This_needed
	@CONFIG=$@;export CONFIG;MAKEFLAGS='';export MAKEFLAGS;$(MAKE) -f This.mk $(t) $(target)  \
"a_or_so       = so" \
"so            = so" \
"CC            = gcc" \
"CFLAGS        = -ansi -DUNIX -DHAS_DLD -O" \
"CCLD          = gcc" \
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
"libXm         = -L/usr/X11R6/lib -lXm" \
$(This_macros)

# Linux PC-Barrand.lal.in2p3.fr 2.0.30 #9 Mon Sep 22 17:53:05 MET DST 1997 i686 unknown
#"CCLDFLAGS     = -Wl,-Bstatic"
# Xm from XiG (previous Xinside).
Linux-gxx : This_needed
	@CONFIG=$@;export CONFIG;MAKEFLAGS='';export MAKEFLAGS;$(MAKE) -f This.mk $(t) $(target)  \
"a_or_so       = so" \
"so            = so" \
"CC            = g++" \
"CFLAGS        = -x c++ -Wall -ansi -DUNIX -DHAS_DLD -O" \
"CCCC          = -x c" \
"CCLD          = g++" \
"CCLDFLAGS     = " \
"CCLDEND       = " \
"CXX           = g++" \
"CXXFLAGS      = -Wall -ansi -pipe -fno-for-scope -DGNU_GCC -DG4_HAVE_BOOL -O" \
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
"JAVA_CPPFLAGS = -I/usr/local/JDK/jdk-1.1.3/include/ -I/usr/local/JDK/jdk-1.1.3/include/genunix" \
"libc          = -lc" \
"libm          = -lm" \
"libf77        = -lf" \
"libX11        = -L/usr/X11R6/lib -lX11 -lSM -lICE" \
"libXext       = -L/usr/X11R6/lib -lXext" \
"libXt         = -L/usr/X11R6/lib -lXt" \
"libXmu        = -L/usr/X11R6/lib -lXmu" \
"libXaw        = -L/usr/X11R6/lib -lXaw" \
"libXm         = -L/usr/X11R6/lib -lXm" \
$(This_macros)


# AIX sp070 1 4 000171755700
# -D_ALL_SOURCE -D_BSD so that sys/types.h includes sys/select.h and then defines fd_set.
# Pb with malloc on AIX 3.2., do a : 
#   csh> setenv  MALLOCTYPE "3.1" 
# to use older version.
# -xlf90 seems now needed in addition to -lxlf.
# On some system use -bI:/usr/lpp/xlf/lib/lowsys.exp to solve reference .__divus, .__quous, .__mulh
# Can do a find/grep over *exp files to find a symbol.
AIX : This_needed
	@CONFIG=$@;export CONFIG;MAKEFLAGS='';export MAKEFLAGS;$(MAKE) -f This.mk $(t) $(target)  \
"a_or_so       = a" \
"so            = so" \
"CC            = xlc" \
"CFLAGS        = -D_POSIX_SOURCE -D_ALL_SOURCE -D_BSD -DUNIX -DHAS_DLD -O" \
"CCLD          = xlc" \
"CCLDFLAGS     = " \
"CCLDEND       = -lld" \
"F77           = xlf" \
"F77FLAGS      = -w -qextname -qrndsngl" \
"F77LD         = xlf" \
"X11_CPPFLAGS  = " \
"Xext_CPPFLAGS = -I/usr/include/X11/extensions" \
"Xmu_CPPFLAGS  = -I/usr/include/X11" \
"Xt_CPPFLAGS   = " \
"Xaw_CPPFLAGS  = " \
"Xm_CPPFLAGS   = " \
"libc          = -lc" \
"libm          = -lm" \
"libf77        = -lxlf90 -lxlf" \
"libX11        = -L/usr/lib -lX11" \
"libXext       = -L/usr/lib -lXext" \
"libXt         = -L/usr/lib -lXt" \
"libXmu        = -L/usr/lib -lXmu" \
"libXaw        = -L/usr/lib -lXaw" \
"libXm         = -L/usr/lib -lXm" \
$(This_macros)

# IRIX fnpspa 5.3 11091810 IP12 mips 
IRIX : This_needed
	@CONFIG=$@;export CONFIG;MAKEFLAGS='';export MAKEFLAGS;$(MAKE) -f This.mk $(t) $(target)  \
"a_or_so       = so" \
"so            = so" \
"CC            = cc" \
"CFLAGS        = -w -DUNIX -DHAS_DLD -O" \
"CCLD          = cc" \
"CCLDFLAGS     = " \
"CCLDEND       = " \
"CXX           = CC" \
"CXXFLAGS      = -ptused -O" \
"CXXLD         = CC" \
"CXXLDFLAGS    = " \
"CXXLDEND      = " \
"F77           = f77" \
"F77FLAGS      = -w +ppu" \
"F77LD         = f77" \
"X11_CPPFLAGS  = " \
"Xext_CPPFLAGS = -I/usr/include/X11/extensions" \
"Xmu_CPPFLAGS  = -I/usr/include/X11" \
"Xt_CPPFLAGS   = " \
"Xaw_CPPFLAGS  = " \
"Xm_CPPFLAGS   = " \
"libc          = -lc" \
"libm          = -lm" \
"libf77        = -lftn -lm -lc" \
"libX11        = -L/usr/lib -lX11" \
"libXext       = -L/usr/lib -lXext" \
"libXt         = -L/usr/lib -lXt" \
"libXmu        = -L/usr/lib -lXmu" \
"libXaw        = -L/usr/lib -lXaw" \
"libXm         = -L/usr/lib -lXm" \
$(This_macros)

# NOT TESTED YET.
# In general SunOS is a nightmare ; cc is not ANSI by default and MOTIF is not usaully in standard place.
# We do not attempt to have shared libs !
SunOS : This_needed
	@CONFIG=$@;export CONFIG;MAKEFLAGS='';export MAKEFLAGS;$(MAKE) -f This.mk $(t) $(target)  \
"a_or_so       = a" \
"so            = so" \
"CC            = acc" \
"CFLAGS        = -DXTFUNCPROTO -DCOSPRINTF -DUNIX -DHAS_DLD -O" \
"CCLD          = acc" \
"CCLDFLAGS     = " \
"CCLDEND       = " \
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

# LynxOS virgo_vme5.lal.in2p3.fr 2.2.2 2 VGPW2
# RAND_MAX not found in /usr/include files.
# -lc_p to find mktime which is not in libc.a.
# If gcc error 11, try to set TMPDIR to a directory 
# with more place :
#   UNIX> setenv TMPDIR `pwd`
LynxOS : This_needed
	@CONFIG=$@;export CONFIG;MAKEFLAGS='';export MAKEFLAGS;$(MAKE) -f This.mk $(t) $(target)  \
"a_or_so       = a" \
"so            = so" \
"CC            = gcc" \
"CFLAGS        = -ansi -DUNIX -DHAS_DLD -DRAND_MAX=32767" \
"CCLD          = gcc" \
"CCLDFLAGS     = " \
"CCLDEND       = -lnetinet -lrpc -lc_p" \
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
"libf77        = -lf" \
"libX11        = -L/usr/X11R6/lib -lX11" \
"libXext       = -L/usr/X11R6/lib -lXext" \
"libXt         = -L/usr/X11R6/lib -lXt" \
"libXmu        = -L/usr/X11R6/lib -lXmu" \
"libXaw        = -L/usr/X11R6/lib -lXaw" \
"libXm         = -L/usr/X11R6/lib -lXm" \
$(This_macros)

