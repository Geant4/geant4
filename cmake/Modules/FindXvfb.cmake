# Check if Xvfb is installed.
#  Xvfb is an X server that can run on machines with no display hardware 
#  and no physical input devices.  
#  It emulates a dumb framebuffer using virtual memory.
# It is used to test graphic examples.

# Usage of this module is as follows
#
#   find_package( Xvfb )
#   if (XVFB_FOUND)
#     ...
# sets XVFB_EXECUTABLE


set(XVFB_FOUND FALSE)

if(UNIX)
   find_program(XVFB_EXECUTABLE Xvfb
	             PATHS /usr/bin
					 DOC "Xvfb dumb framebuffer X server")
endif()

# handle the QUIETLY and REQUIRED arguments and set PYTHIA6_FOUND to TRUE if
# all listed variables are TRUE
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(Xvfb DEFAULT_MSG XVFB_EXECUTABLE)
				 
mark_as_advanced(XVFB_FOUND XVFB_EXECUTABLE)
