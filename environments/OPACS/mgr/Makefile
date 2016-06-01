#/* +---------------------- Copyright notice -------------------------------+ */
#/* | Copyright (C) 1995, Guy Barrand, LAL Orsay. (barrand@lalcls.in2p3.fr) | */
#/* |   Permission to use, copy, modify, and distribute this software       | */
#/* |   and its documentation for any purpose and without fee is hereby     | */
#/* |   granted, provided that the above copyright notice appear in all     | */
#/* |   copies and that both that copyright notice and this permission      | */
#/* |   notice appear in supporting documentation.  This software is        | */
#/* |   provided "as is" without express or implied warranty.               | */
#/* +---------------------- Copyright notice -------------------------------+ */

This_macros     = \
"bin            = ../$${CONFIG}" \
"mgr            = ../mgr" \
"src            = ../src" \
"app            = ../app" \
"make_so        = sh ${COROOT}/mgr/maklibs.sh" \
"G4o_CPPFLAGS   = -I${COROOT}/include -I${GOROOT}/include -I../include -DHAS_HO -DHAS_G4O" \
"libG4o         = -L../$${CONFIG}  -lG4o" \
"RW_CPPFLAGS    = -I${RWROOT}/include" \
"librw          = -L${RWROOT}/lib -lrwtool" \
"CLHEP_CPPFLAGS = -I${CLHEPROOT}/include " \
"libCLHEP       = -L${CLHEPROOT}/lib -lCLHEP" \
"GL_CPPFLAGS    = -I${GLOROOT}/usr/include" \
"libGL          = -L${GLOROOT}/$${CONFIG} -lGLo" \
"libGLU         = -L${GLOROOT}/$${CONFIG} -lGLUo" \
"libGLX         = -L${GLOROOT}/$${CONFIG} -lGLXo" \
"G4_CPPFLAGS = -I${G4INSTALL}/source/global/management/include -I${G4INSTALL}/source/global/HEPRandom/include -I${G4INSTALL}/source/global/HEPGeometry/include -I${G4INSTALL}/source/run/include -I${G4INSTALL}/source/materials/include -I${G4INSTALL}/source/particles/management/include -I${G4INSTALL}/source/particles/bosons/include -I${G4INSTALL}/source/particles/leptons/include -I${G4INSTALL}/source/particles/hadrons/barions/include -I${G4INSTALL}/source/particles/hadrons/ions/include -I${G4INSTALL}/source/particles/hadrons/mesons/include -I${G4INSTALL}/source/processes/management/include -I${G4INSTALL}/source/processes/general/include -I${G4INSTALL}/source/processes/decay/include -I${G4INSTALL}/source/processes/optical/include -I${G4INSTALL}/source/processes/transportation/include -I${G4INSTALL}/source/processes/parameterisation/include -I${G4INSTALL}/source/processes/electromagnetic/standard/include -I${G4INSTALL}/source/processes/electromagnetic/muons/include -I${G4INSTALL}/source/processes/electromagnetic/utils/include -I${G4INSTALL}/source/processes/electromagnetic/xrays/include -I${G4INSTALL}/source/processes/hadronic/management/include -I${G4INSTALL}/source/processes/hadronic/processes/include -I${G4INSTALL}/source/processes/hadronic/util/include -I${G4INSTALL}/source/processes/hadronic/cross_sections/include -I${G4INSTALL}/source/processes/hadronic/models/generator/include -I${G4INSTALL}/source/processes/hadronic/models/high_energy/include -I${G4INSTALL}/source/processes/hadronic/models/low_energy/include -I${G4INSTALL}/source/processes/hadronic/stopping/include -I${G4INSTALL}/source/geometry/management/include -I${G4INSTALL}/source/geometry/solids/CSG/include -I${G4INSTALL}/source/geometry/solids/BREPS/include -I${G4INSTALL}/source/geometry/solids/STEP/include -I${G4INSTALL}/source/geometry/volumes/include -I${G4INSTALL}/source/geometry/magneticfield/include -I${G4INSTALL}/source/track/include -I${G4INSTALL}/source/tracking/include -I${G4INSTALL}/source/digits+hits/detector/include -I${G4INSTALL}/source/digits+hits/hits/include -I${G4INSTALL}/source/digits+hits/digits/include -I${G4INSTALL}/source/readout/include -I${G4INSTALL}/source/event/include -I${G4INSTALL}/source/graphics_reps/include -I${G4INSTALL}/source/persistency/include -I${G4INSTALL}/source/intercoms/include -I${G4INSTALL}/source/visualization/modeling/include" \
"libg4 = -L${G4INSTALL}/lib/$${G4SYSTEM} -lG4modeling -lG4readout -lG4run -lG4event -lG4tracking -lG4particles -lG4processes -lG4track -lG4materials -lG4digits+hits -lG4geometry -lG4graphics_reps -lG4intercoms -lG4global" \
"G4TEST_CPPFLAGS = -I${G4INSTALL}/source/visualization/test/test_utils/include" \
"libg4test       = -L${G4INSTALL}/lib/${G4SYSTEM} -lG4vis_utils" \
"G4ATLAS_CPPFLAGS = -I${G4ATLAS}//include" \
"libg4atlas      = -L${G4INSTALL}/lib/${G4SYSTEM} -lG4ATLAS" \
"Wo_CPPFLAGS     = -I${COROOT}/include -I${XXROOT}/include -I${WOROOT}/include" \
"libCo           = -L${COROOT}/$${CONFIG} -lCo" \
"libXx           = -L${XXROOT}/$${CONFIG} -lXx" \
"libWo           = -L${WOROOT}/$${CONFIG} -lWo" \
"libWoGL         = -L${WOROOT}/$${CONFIG} -lWoGL" \
"WoXm_CPPFLAGS   = -DHAS_XM" \
"libWoXm         = -L${WOROOT}/$${CONFIG} -lWoXm" \
"Xo_CPPFLAGS     = -DHAS_XO -DHAS_HTML -I${XOROOT}/include -I${GOROOT}/include -I${HOROOT}/include" \
"libXo           = -L${XOROOT}/$${CONFIG}  -lXo" \
"libHo           = -L${HOROOT}/$${CONFIG} -lHo" \
"libGo           = -L${GOROOT}/$${CONFIG}  -lGo" \
"libwww          = -L${W3OROOT}/$${CONFIG} -lW3o" \
"libhtmlw        = -L${HTMLOROOT}/$${CONFIG} -lHTMLo"

include Config.mk


This_needed :
	@if test "${COROOT}" = "" ; then echo "COROOT not set" ; exit 1 ; else exit 0 ; fi
	@if test "${GOROOT}" = "" ; then echo "GOROOT not set"; exit 1 ; else exit 0 ; fi
	@if test "${G4INSTALL}" = "" ; then echo "G4INSTALL not set" exit 1 ; else exit 0 ; fi
	@
	@if test "${GLOROOT}" = "" ; then echo "GLOROOT not set" ; else exit 0 ; fi


