G4LIB_BUILD_SHARED = 1

name := Brachy

extrabuild:
	g++ -shared -o${G4WORKDIR}/tmp/${G4SYSTEM}/${name}/libDIANE_${name}.so -L${G4WORKDIR}/tmp/${G4SYSTEM}/${name} -l${name} ${LDFLAGS} ${LDLIBS} 

include $(G4INSTALL)/config/binmake.gmk
