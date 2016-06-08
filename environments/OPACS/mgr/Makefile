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
"RW_CPPFLAGS    = ${RWFLAGS}" \
"librw          = ${RWLIBS}" \
"CLHEP_CPPFLAGS = -I${CLHEPROOT}/include " \
"libCLHEP       = -L${CLHEPROOT}/lib -lCLHEP" \
"GL_CPPFLAGS    = -I${GLOROOT}/usr/include" \
"libGL          = -L${GLOROOT}/${OCONFIG} -lGLo" \
"libGLU         = -L${GLOROOT}/${OCONFIG} -lGLUo" \
"libGLX         = -L${GLOROOT}/${OCONFIG} -lGLXo" \
"Wo_CPPFLAGS     = -I${COROOT}/include -I${XXROOT}/include -I${WOROOT}/include" \
"libCo           = -L${COROOT}/${OCONFIG} -lCo" \
"libXx           = -L${XXROOT}/${OCONFIG} -lXx" \
"libWo           = -L${WOROOT}/${OCONFIG} -lWo" \
"libWoGL         = -L${WOROOT}/${OCONFIG} -lWoGL" \
"WoXm_CPPFLAGS   = -DHAS_XM" \
"libWoXm         = -L${WOROOT}/${OCONFIG} -lWoXm" \
"Xo_CPPFLAGS     = -DHAS_XO -DHAS_HTML -I${XOROOT}/include -I${GOROOT}/include -I${HOROOT}/include" \
"libXo           = -L${XOROOT}/${OCONFIG}  -lXo" \
"libHo           = -L${HOROOT}/${OCONFIG} -lHo" \
"libGo           = -L${GOROOT}/${OCONFIG}  -lGo" \
"libwww          = -L${W3OROOT}/${OCONFIG} -lW3o" \
"libhtmlw        = -L${HTMLOROOT}/${OCONFIG} -lHTMLo" \
"G4_CPPFLAGS = -DG4USE_STL -DG4USE_NAMESPACE -DG4USE_STD_NAMESPACE -I${G4INSTALL}/source/global/management/include -I${G4INSTALL}/source/global/HEPRandom/include -I${G4INSTALL}/source/global/HEPGeometry/include -I${G4INSTALL}/source/run/include -I${G4INSTALL}/source/materials/include -I${G4INSTALL}/source/particles/management/include -I${G4INSTALL}/source/particles/bosons/include -I${G4INSTALL}/source/particles/leptons/include -I${G4INSTALL}/source/particles/hadrons/barions/include -I${G4INSTALL}/source/particles/hadrons/ions/include -I${G4INSTALL}/source/particles/hadrons/mesons/include -I${G4INSTALL}/source/processes/management/include -I${G4INSTALL}/source/processes/general/include -I${G4INSTALL}/source/processes/decay/include -I${G4INSTALL}/source/processes/optical/include -I${G4INSTALL}/source/processes/transportation/include -I${G4INSTALL}/source/processes/parameterisation/include -I${G4INSTALL}/source/processes/electromagnetic/standard/include -I${G4INSTALL}/source/processes/electromagnetic/muons/include -I${G4INSTALL}/source/processes/electromagnetic/utils/include -I${G4INSTALL}/source/processes/electromagnetic/xrays/include -I${G4INSTALL}/source/processes/hadronic/management/include -I${G4INSTALL}/source/processes/hadronic/processes/include -I${G4INSTALL}/source/processes/hadronic/util/include -I${G4INSTALL}/source/processes/hadronic/cross_sections/include -I${G4INSTALL}/source/processes/hadronic/models/generator/include -I${G4INSTALL}/source/processes/hadronic/models/high_energy/include -I${G4INSTALL}/source/processes/hadronic/models/low_energy/include -I${G4INSTALL}/source/processes/hadronic/stopping/include -I${G4INSTALL}/source/geometry/management/include -I${G4INSTALL}/source/geometry/solids/CSG/include -I${G4INSTALL}/source/geometry/solids/specific/include -I${G4INSTALL}/source/geometry/solids/BREPS/include -I${G4INSTALL}/source/geometry/solids/STEP/include -I${G4INSTALL}/source/geometry/volumes/include -I${G4INSTALL}/source/geometry/magneticfield/include -I${G4INSTALL}/source/track/include -I${G4INSTALL}/source/tracking/include -I${G4INSTALL}/source/digits+hits/detector/include -I${G4INSTALL}/source/digits+hits/hits/include -I${G4INSTALL}/source/digits+hits/digits/include -I${G4INSTALL}/source/readout/include -I${G4INSTALL}/source/event/include -I${G4INSTALL}/source/graphics_reps/include -I${G4INSTALL}/source/persistency/include -I${G4INSTALL}/source/intercoms/include -I${G4INSTALL}/source/visualization/modeling/include" \
"libg4 = -L${G4LIB}/$${G4SYSTEM} -lG4modeling -lG4emstandard -lG4decay -lG4UIbasic -lG4muons -lG4ions -lG4mesons -lG4leptons -lG4baryons -lG4run -lG4event -lG4bosons -lG4tracking -lG4transportation -lG4detector -lG4emutils -lG4hits -lG4procman -lG4digits -lG4track -lG4UIcommon -lG4partman -lG4volumes -lG4magneticfield -lG4geomBoolean -lG4csg -lG4geometrymng -lG4materials -lG4intercoms -lG4graphics_reps -lG4globman"

# With Geant4 global libs use :
#"libg4 = -L${G4LIB}/$${G4SYSTEM} -lG4modeling -lG4run -lG4event -lG4tracking -lG4processes -lG4particles -lG4digits+hits -lG4track -lG4materials -lG4geometry -lG4graphics_reps -lG4intercoms -lG4global"

include Config.mk

This_needed :
	@if test "${COROOT}" = "" ; then echo "COROOT not set" ; exit 1 ; else exit 0 ; fi
	@if test "${GOROOT}" = "" ; then echo "GOROOT not set"; exit 1 ; else exit 0 ; fi
	@if test "${G4INSTALL}" = "" ; then echo "G4INSTALL not set" exit 1 ; else exit 0 ; fi
	@if test "${G4LIB}" = "" ; then echo "G4LIB not set" exit 1 ; else exit 0 ; fi
	@
	@if test "${GLOROOT}" = "" ; then echo "GLOROOT not set" ; else exit 0 ; fi


