#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#  File produced by the Co/omake tool
# using file Make.odb.
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#  The soft linked packages has been taken
# from the property 'has' of the 
# 'Make' object of file Make.odb.
# It had the value :
#    X11 Xext Xt Xmu Xm GL RW CLHEP G4 G4TEST G4ATLAS Wo WoXm Xo
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
HAS_CPPFLAGS = \
 $(G4o_CPPFLAGS)\
 $(X11_CPPFLAGS)\
 $(Xext_CPPFLAGS)\
 $(Xt_CPPFLAGS)\
 $(Xmu_CPPFLAGS)\
 $(Xm_CPPFLAGS)\
 $(GL_CPPFLAGS)\
 $(RW_CPPFLAGS)\
 $(CLHEP_CPPFLAGS)\
 $(G4_CPPFLAGS)\
 $(G4TEST_CPPFLAGS)\
 $(G4ATLAS_CPPFLAGS)\
 $(Wo_CPPFLAGS)\
 $(WoXm_CPPFLAGS)\
 $(Xo_CPPFLAGS)
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
all   :: mkdir

clean ::
	/bin/rm -f $(bin)/*.o;/bin/rm -f $(bin)/*.exe;/bin/rm -f $(bin)/*.class;/bin/rm -f $(bin)/*.a;/bin/rm -f $(bin)/*.so;/bin/rm -f $(bin)/*.sl
rmlib ::
	/bin/rm -f $(bin)/*.a
rmo   ::
	/bin/rm -f $(bin)/*.o
rmexe ::
	/bin/rm -f $(bin)/*.exe;/bin/rm -f $(bin)/*.class

mkdir :
	@if test -d $(bin) ; then exit ; else mkdir $(bin) ; echo "$(bin) created." ; fi

libG4o_target = $(bin)/libG4o.a
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
all :: mkdir \
$(libG4o_target) \
$(bin)/EXPO.exe 
	@echo "G4o : all ok." 

libs :: mkdir \
$(libG4o_target) 
	@echo "G4o : libs ok." 

apps :: mkdir \
$(bin)/EXPO.exe \
$(bin)/TEST.exe 
	@echo "G4o : apps ok." 
#--------------------------------------------
rmexeo :
	/bin/rm -f $(bin)/EXPO.o
	/bin/rm -f $(bin)/TEST.o
#--------------------------------------------
EXPO : $(bin)/EXPO.exe
	@echo "G4o : EXPO ok."
TEST : $(bin)/TEST.exe
	@echo "G4o : TEST ok."
libG4o : $(libG4o_target)
	@echo "G4o : libG4o ok."
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#--------------------------------------------
$(bin)/EXPO.exe : $(bin)/EXPO.o \
$(libG4o_target) 
	$(CXXLD) $(CXXLDFLAGS) $(HAS_CPPFLAGS) -o $(bin)/EXPO.exe $(bin)/EXPO.o \
$(libG4o) \
$(libg4) \
$(libCLHEP) \
$(librw) \
$(libWoGL) \
$(libXo) \
$(libhtmlw) \
$(libWoXm) \
$(libXm) \
$(libWo) \
$(libXx) \
$(libXmu) \
$(libXt) \
$(libGo) \
$(libGLX) \
$(libGLU) \
$(libGL) \
$(libXext) \
$(libX11) \
$(libHo) \
$(libCo) \
$(libwww) \
$(libm) \
$(CXXLDEND) 

$(bin)/EXPO.o   : $(app)/EXPO.cc
	$(CXX) $(CXXFLAGS) $(APP_CXXFLAGS) $(HAS_CPPFLAGS) -c -o $(bin)/EXPO.o $(app)/EXPO.cc
#--------------------------------------------
$(bin)/TEST.exe : $(bin)/TEST.o \
$(libG4o_target) 
	$(CXXLD) $(CXXLDFLAGS) $(HAS_CPPFLAGS) -o $(bin)/TEST.exe $(bin)/TEST.o \
$(libg4test) \
$(libg4atlas) \
$(libG4o) \
$(libg4) \
$(libCLHEP) \
$(librw) \
$(libWoGL) \
$(libXo) \
$(libhtmlw) \
$(libWoXm) \
$(libXm) \
$(libWo) \
$(libXx) \
$(libXmu) \
$(libXt) \
$(libGo) \
$(libGLX) \
$(libGLU) \
$(libGL) \
$(libXext) \
$(libX11) \
$(libHo) \
$(libCo) \
$(libwww) \
$(libm) \
$(CXXLDEND) 

$(bin)/TEST.o   : $(app)/TEST.cc
	$(CXX) $(CXXFLAGS) $(APP_CXXFLAGS) $(HAS_CPPFLAGS) -c -o $(bin)/TEST.o $(app)/TEST.cc
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# libraries 
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
.PRECIOUS : $(bin)/libG4o.a 
#--------------------------------------------

$(bin)/libG4o.a : \
$(bin)/libG4o.a(G4o.o) \
$(bin)/libG4o.a(G4oScene.o) \
$(bin)/libG4o.a(G4oDrawer.o) \
$(bin)/libG4o.a(G4oState.o) 
	@cd $(bin) ; if [ -f /bin/ranlib ] ; then /bin/ranlib libG4o.a ; fi ; cd $(mgr)

#--------------------------------------------
# libG4o dependencies 
#--------------------------------------------
$(bin)/libG4o.a(G4o.o) : $(src)/G4o.cc
	$(CXX) $(CXXFLAGS) $(HAS_CPPFLAGS) -c -o $(bin)/G4o.o $(src)/G4o.cc
	ar cr $(bin)/libG4o.a $(bin)/G4o.o
	/bin/rm -f $(bin)/G4o.o

$(bin)/libG4o.a(G4oScene.o) : $(src)/G4oScene.cc
	$(CXX) $(CXXFLAGS) $(HAS_CPPFLAGS) -c -o $(bin)/G4oScene.o $(src)/G4oScene.cc
	ar cr $(bin)/libG4o.a $(bin)/G4oScene.o
	/bin/rm -f $(bin)/G4oScene.o

$(bin)/libG4o.a(G4oDrawer.o) : $(src)/G4oDrawer.cc
	$(CXX) $(CXXFLAGS) $(HAS_CPPFLAGS) -c -o $(bin)/G4oDrawer.o $(src)/G4oDrawer.cc
	ar cr $(bin)/libG4o.a $(bin)/G4oDrawer.o
	/bin/rm -f $(bin)/G4oDrawer.o

$(bin)/libG4o.a(G4oState.o) : $(src)/G4oState.cc
	$(CXX) $(CXXFLAGS) $(HAS_CPPFLAGS) -c -o $(bin)/G4oState.o $(src)/G4oState.cc
	ar cr $(bin)/libG4o.a $(bin)/G4oState.o
	/bin/rm -f $(bin)/G4oState.o


