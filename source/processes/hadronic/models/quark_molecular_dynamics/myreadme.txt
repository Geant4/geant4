0 - To update the checked out version of Geant4, use first
      cvs log GNUmakefile | less
    on the main Makefile (G4INSTAL/source/GNUmakefile)
    to get the status of the development.
    Chose the latest released version, e.g.
      geant4-03-00: 1.25
    and check it out by 
      cvs update -r geant4-03-00
    in G4INSTAL/source and G4INSTAL/config
    Important:
    save current development with
      cvs commit -m "current status"
    restore developping tree in relevant directories as
    geant4/source/processes/hadronic/models/generator/
    geant4/source/processes/hadronic/models/quark_molecular_dynamics
    with 
      cvs -n update -A 
    (-A: head revision - latest changes; -n: do nothing)


1 - make sure that variables 
     G4INSTALL=/afs/cern.ch/user/s/sscherer/sungeant/geant4
     G4SYSTEM=SUN-CC5
    are set (set in .bashrc on sungeant)

2 - compile "util" and "body" by calling gmake in these directries
    this creates object and library files in 
      geant4/tmp/SUN-CC/G4hadronic_quark_md_body 
    and 
      geant4/tmp/SUN-CC/G4hadronic_quark_md_util
    and finally creates object files 
      geant4/lib/SUN-CC/libG4hadronic_quark_md_body.a
    and 
      geant4/lib/SUN-CC/libG4hadronic_quark_md_util.a

    flag G4DEBUG may be necessary to set: export G4DEBUG=1 

3 - set variable TESTTARGET to source code in "test", eg. 
      export TESTTARGET=qmd_fromG4String
    clear the variable with export -n TESTTARGET or export -d TESTTARGET

4 - call gmake in directory "test"

5 - changing some files in "body" or "util" may yield complete 
    recompilation only after "gmake clean"

6 - If linking fails, checking the existence of all used libraries may
    solve the problem, 
    e.g. calling gmake in source/global/management creates the library 
      lib/SUN-CC/libG4globman.a
    names of libaries are defined in the "name := theName" statement of
    the GNUmakefile, so if the linking produvces the error
      ld: fatal: library -lG4baryons: not found 
    the library libG4baryons.a is missing. It can be created by starting
    gmake in the directory found by calling
      find . -name GNUmakefile | xargs grep "name := G4baryons"
    in $G4INSTALL.

7 - If there are "undefined symbols first referenced in..."
    check whether
      gmake libmap
    has been done after compiling quark_md_util/body
    If this does not help...
    recompile the whole stuff...

8 - binary will be found in geant4/bin/SUN-CC

