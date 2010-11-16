#!/bin/sh
cat << EoI | grep -v -e "^#"

	
Random_FAIL exam-ext-field-field04-run \
	slc4_amd64_gcc34_dbg	\
	i686-slc5-gcc43-opt
	
EXPECTED_FAIL exam-ext-g3tog4-clgeometry-run \
       osx106_ia32_gcc401      \
       osx106_ia32_gcc401_dbg
		       
EXPECTED_FAIL exam-ext-g3tog4-cltog4-run  \
       osx106_ia32_gcc401      \
       osx106_ia32_gcc401_dbg



#28 Apr back Jul 10 on slc4_ia32_gcc34
EXPECTED_FAIL test29-run \
	slc4_ia32_gcc34

	


#newly in  testing Jul 2010, fail to run for graphics so far....
EXPECTED_FAIL test202-build \
	i686-winxp-vc9-opt	\
	i686-winxp-vc9-dbg
	
#new in testing Oct 2010, ok 11-Nov
#EXPECTED_FAIL exam-ext-analysis-n03con-run

#new in testing Oct 2010
EXPECTED_FAIL exam-ext-biasing-ReverseMC01-build \
	x86_64-slc5-gcc43-opt	\
	x86_64-slc5-gcc43-dbg	
	
#new in testing Oct 2010
EXPECTED_FAIL exam-ext-biasing-ReverseMC01-run \
	slc4_amd64_gcc34	\
	slc4_amd64_gcc34_dbg	\
	slc4_ia32_gcc34		\
	slc4_ia32_gcc34_dbg	\
	x86_64-slc5-gcc41-opt	\
	x86_64-slc5-gcc41-dbg	\
	x86_64-mac106-gcc42-deb
	
	
EoI

exit

slc4_amd64_gcc34
slc4_amd64_gcc34_dbg

slc4_ia32_gcc34
slc4_ia32_gcc34_dbg

i686-slc5-gcc43-opt
i686-slc5-gcc43-dbg

x86_64-slc5-gcc43-opt
x86_64-slc5-gcc43-dbg

osx105_ia32_gcc401
osx105_ia32_gcc401_dbg

i686-winxp-vc9-dbg

	slc4_amd64_gcc34	\
	slc4_amd64_gcc34_dbg	\
	slc4_ia32_gcc34		\
	slc4_ia32_gcc34_dbg	\
	i686-slc5-gcc43-opt	\
	i686-slc5-gcc43-dbg	\
	x86_64-slc5-gcc43-opt	\
	x86_64-slc5-gcc43-dbg	\
	x86_64-slc5-gcc41-opt	\
	x86_64-slc5-gcc41-dbg	\
	x86_64-mac106-gcc42-opt	\
	x86_64-mac106-gcc42-deb
