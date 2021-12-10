
# start from nothing, so repeated inclusion of this into CMake context doesnt repeat the flags 
set(CMAKE_CXX_FLAGS)

if(WIN32)

  # need to detect compiler not os?
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -W4") # overall warning level 4
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -wd4996")   # disable  C4996: 'strdup': The POSIX name for this item is deprecated.
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DNOMINMAX")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_SCL_SECURE_NO_WARNINGS")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_CRT_SECURE_NO_WARNINGS")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_USE_MATH_DEFINES")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_ITERATOR_DEBUG_LEVEL=0")


else(WIN32)

  ## c++11 forced by AsioZMQ : AsioZMQ not used here, but expect best to use same compiler options as far as possible
  #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -std=c++0x")  ## huh nvcc compilation fails with this ???
  if (${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang")
     # set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -std=c++11 -stdlib=libc++")
     set(CMAKE_CXX_STANDARD 14)
     set(CMAKE_CXX_STANDARD_REQUIRED on)
  else ()
      #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
      #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -std=c++0x")
     # set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -std=c++11")   #needed for numpyserver- on Linux ?
     set(CMAKE_CXX_STANDARD 14)
     set(CMAKE_CXX_STANDARD_REQUIRED on)

  endif ()

  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fvisibility=hidden")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fvisibility-inlines-hidden") ## avoid boostrap visibility warning at link 
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fdiagnostics-show-option") 
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-function")

  if (${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
     set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-comment")
     set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-deprecated")
     set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-shadow")
  else()
     set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-private-field")
     set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-shadow")
  endif()

endif(WIN32)


if(FLAGS_VERBOSE)
   # https://cmake.org/Wiki/CMake_Useful_Variables
   message(STATUS "OpticksCompilationFlags.cmake : CMAKE_BUILD_TYPE = ${CMAKE_BUILD_TYPE}")
   message(STATUS "OpticksCompilationFlags.cmake : CMAKE_CXX_FLAGS = ${CMAKE_CXX_FLAGS}")
   message(STATUS "OpticksCompilationFlags.cmake : CMAKE_CXX_FLAGS_DEBUG = ${CMAKE_CXX_FLAGS_DEBUG}")
   message(STATUS "OpticksCompilationFlags.cmake : CMAKE_CXX_FLAGS_RELEASE = ${CMAKE_CXX_FLAGS_RELEASE}")
   message(STATUS "OpticksCompilationFlags.cmake : CMAKE_CXX_FLAGS_RELWITHDEBINFO= ${CMAKE_CXX_FLAGS_RELWITHDEBINFO}")
   message(STATUS "OpticksCompilationFlags.cmake : CMAKE_CXX_STANDARD : ${CMAKE_CXX_STANDARD} " )
   message(STATUS "OpticksCompilationFlags.cmake : CMAKE_CXX_STANDARD_REQUIRED : ${CMAKE_CXX_STANDARD_REQUIRED} " )
endif()



