#.rst:
# G4GitUtilities
# --------------
#
# .. code-block::cmake
#
#   include(G4GitUtilities)
#
# CMake functions and macros for working with Git where needed in
# Geant4's build/packaging systems.
#

#-----------------------------------------------------------------
# License and Disclaimer
#
# The  Geant4 software  is  copyright of the Copyright Holders  of
# the Geant4 Collaboration.  It is provided  under  the terms  and
# conditions of the Geant4 Software License,  included in the file
# LICENSE and available at  http://cern.ch/geant4/license .  These
# include a list of copyright holders.
#
# Neither the authors of this software system, nor their employing
# institutes,nor the agencies providing financial support for this
# work  make  any representation or  warranty, express or implied,
# regarding  this  software system or assume any liability for its
# use.  Please see the license in the file  LICENSE  and URL above
# for the full disclaimer and the limitation of liability.
#
# This  code  implementation is the result of  the  scientific and
# technical work of the GEANT4 collaboration.
# By using,  copying,  modifying or  distributing the software (or
# any work based  on the software)  you  agree  to acknowledge its
# use  in  resulting  scientific  publications,  and indicate your
# acceptance of all terms of the Geant4 Software license.
#
#-----------------------------------------------------------------

include_guard(GLOBAL)

#-----------------------------------------------------------------------
#.rst:
# .. cmake:command:: geant4_git_find_dirty
#
#   .. code-block:: cmake
#
#     geant4_git_find_dirty(<basedir> <modified> <untracked>)
#
#   Runs ``git status`` in the ``<basedir>``, writing lists of modified
#   and untracked/ignored files to ``<modified>`` and ``<untracked>``
#   respectively.
#
#   ``<basedir>`` must be the root of the repository queried, i.e. containing
#   the ``.git/`` directory.
#
#   Nothing is returned in the output variables if ``<basedir>`` is not a
#   git repository or if a ``git`` executable cannot be located.
#
function(geant4_git_find_dirty _basedir _output_modified _output_untracked)
  if(NOT EXISTS "${_basedir}/.git")
    return()
  endif()

  find_package(Git QUIET)
  if(NOT Git_FOUND)
    return()
  endif()

  execute_process(COMMAND ${GIT_EXECUTABLE} status -s --ignored
    WORKING_DIRECTORY "${_basedir}"
    OUTPUT_VARIABLE GEANT4_GIT_UNTRACKED
    ERROR_VARIABLE GEANT4_GIT_UNTRACKED_ERROR)
  
  if(GEANT4_GIT_UNTRACKED_ERROR)
    message(FATAL_ERROR "geant4_git_find_dirty: Calling `git status` failed with error: '${GEANT4_GIT_UNTRACKED_ERROR}'")
  endif()

  if(GEANT4_GIT_UNTRACKED)
    set(_modded_files)
    set(_untracked_files)

    string(REPLACE "\n" ";" GEANT4_GIT_UNTRACKED "${GEANT4_GIT_UNTRACKED}")
    foreach(item ${GEANT4_GIT_UNTRACKED})
      if(item MATCHES "^ *(M|A|R|D|C)")
        string(REGEX REPLACE "^ *[MARDC][MARDC] *" "" _modded_path "${item}")
        list(APPEND _modded_files ${_modded_path})
      endif()

      if(item MATCHES "^ *(\\!|\\?)")
        string(REGEX REPLACE "^ *(\\!\\!|\\?\\?) *" "" _untracked_path "${item}")
        list(APPEND _untracked_files "${_untracked_path}")
      endif()
    endforeach()

    set(${_output_modified} ${_modded_files} PARENT_SCOPE)
    set(${_output_untracked} ${_untracked_files} PARENT_SCOPE)
  endif()
endfunction()
