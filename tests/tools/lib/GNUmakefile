ifndef G4INSTALL
  G4INSTALL = ../../..
endif
include $(G4INSTALL)/config/architecture.gmk

CPPFLAGS+=-I../include

SRCS:=$(wildcard ../src/*.cc)
OBJS:=$(patsubst ../src/%.cc,$(G4SYSTEM)/%.o,$(SRCS))
dependencies := $(patsubst ../src/%.cc,$(G4SYSTEM)/%.d,$(SRCS))

############### RULES ##########################################
$(G4SYSTEM)/libG4TestTool.a: $(G4SYSTEM)/obj.last
	@echo Creating/replacing object files in $(G4SYSTEM)/libG4TestTool.a
	@$(AR) $(OUT_LIB)$@ $(G4SYSTEM)/*.o

# Touch the versioning file
$(G4SYSTEM)/obj.last: $(OBJS)
	@touch $@

$(G4SYSTEM)/%.o: ../src/%.cc
	@echo Compiling $*.cc ...
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $< $(OUT_OBJ)$(G4SYSTEM)/$*.o

$(G4SYSTEM)/%.d: ../src/%.cc
	@if [ ! -d $(G4SYSTEM) ] ; then mkdir $(G4SYSTEM) ;fi
	@echo Making dependency for file $< ...
	@($(ECHO) $(G4SYSTEM)/\\c ; \
	g++ -MM -I../include $< ) | sed 's!$(G4SYSTEM)/$*.o!& $@!' >$@
-include $(dependencies)
###########################################################################

.PHONY: all obj lib

obj: $(G4SYSTEM)/obj.last
lib: $(G4SYSTEM)/libG4TestTool.a

clean:
	@echo Cleaning up ...
	@rm -f $(G4SYSTEM)/libG4TestTool.a
	@rm -rf $(G4SYSTEM)

clean_libs:
	@if [ -f $(G4SYSTEM)/libG4TestTool.a ] ; then \
	$(ECHO) Removing library libG4TestTool.a ... ; \
	$(RM) -f $(G4LIBDIR)/libG4TestTool.a ; fi
