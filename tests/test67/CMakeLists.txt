cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(test67)
find_package(Geant4 REQUIRED)
include(${Geant4_USE_FILE})

find_package(ROOT QUIET)
if(ROOT_FOUND)
  #message(STATUS "G4 TESTS: ROOT package found --> test67 with analysis enabled") 
  add_definitions(-DG4_USE_ROOT)
  include_directories(${ROOT_INCLUDE_DIR} ${PROJECT_SOURCE_DIR}/include)  
  # copy .root files to build area
  FILE(GLOB _refs RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.root test67*.in)
  FILE(COPY ${_refs} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
else()
  message(STATUS "G4 TESTS: ROOT package not found. --> test67 Analysis disabled")
  include_directories(${PROJECT_SOURCE_DIR}/include)
endif() 

if(GEANT4_ENABLE_TESTING)
	if(ROOT_FOUND)
	 GEANT4_EXECUTABLE(test67 test67.cc src/*.cc LIBRARIES ${ROOT_LIBRARIES})
	else()
	 GEANT4_EXECUTABLE(test67 test67.cc src/*.cc LIBRARIES)
	endif()
	
	#---Test definitions----------------------------------------------------------------------------------
	# For all groups
	GEANT4_ADD_TEST( test67-build 
			 	         BUILD test67 
					 LABELS Nightly Continuous PhysicsChecks
	  				 ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
	GEANT4_ADD_TEST(test67 COMMAND test67 ${CMAKE_CURRENT_SOURCE_DIR}/test67.in
	                       DEPENDS test67-build
			       LABELS Nightly Continuous PhysicsChecks
	                       ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
	
	#Add special tests for PhysicsList group
	if(ROOT_FOUND)
	   #Define a test that builds the executable: needed independent since tests after will be
	   #executed in parallel. This ensure correct dependency treatment
	   #Loop on all different test types
	   foreach(_confs test67_livermore test67_stdopt3 test67_stdopt0)
		#Use for testname same string as macro name (test names need to be unique)
	  	#Note we need a different working directory for each test so each test can write its own output
		#w/o interference 
		#set( _wd ${CMAKE_BINARY_DIR}/${_confs} )
		#file(MAKE_DIRECTORY ${_wd} )
		GEANT4_ADD_TEST( ${_confs} COMMAND test67 ${CMAKE_CURRENT_SOURCE_DIR}/${_confs}.in
		#			   WORKING_DIRECTORY ${_wd}
					   DEPENDS test67-build
					   LABELS PhysicsChecks
					   ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
					   
	   endforeach()
	
	endif()
else()
	  # allow a build without testing enabled....
	  # Locate sources and headers for this project
	  include_directories(${PROJECT_SOURCE_DIR}/include 
	                  	 ${Geant4_INCLUDE_DIR})
	  file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
	  file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)
	
	  #----------------------------------------------------------------------------
	  # Add the executable, and link it to the Geant4 libraries
	  #
	  add_executable(test67 test67.cc ${sources} ${headers})
	  if(ROOT_FOUND)
	      target_link_libraries(test67 ${Geant4_LIBRARIES} ${ROOT_LIBRARIES})
	  else()
	      target_link_libraries(test67 ${Geant4_LIBRARIES})
	  endif()    
	  install(TARGETS test67 DESTINATION bin)

endif()
