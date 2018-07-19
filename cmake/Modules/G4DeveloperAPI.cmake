#.rst:
# G4DeveloperAPI
# --------------
#
# .. code-block::cmake
#
#   include(G4DeveloperAPI)
#
# CMake functions and macros for declaring and working with build
# products of Geant4.
#
# WIP: At present, this simply forwards to the old implementation
#      Eventually will allow switch between the old and the new,
#      then just the new once transition complete.

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

if(NOT __G4DEVELOPERAPI_INCLUDED)
  set(__G4DEVELOPERAPI_INCLUDED TRUE)
else()
  return()
endif()

include(G4DeveloperAPI_OLD)
