//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
//
//---------------------------------------------------------------------------
//
// ClassName:   G4OpticalProcessIndex
//
// Author:      P.Gumplinger 29.09.2011
//              (based on code from I. Hrivnacova)
//
// Modified:
//
//----------------------------------------------------------------------------
//
// Enumeration for processes defined in G4OpticalPhysics
//
// The enumeration constants are used as the indices of the instatianted
// processes in the vector of processes; needed to configure G4OpticalPhysics
// (PhysicsList) according to selected process choices.
//

#ifndef G4OpticalProcessIndex_h
#define G4OpticalProcessIndex_h 1

#include "globals.hh"

enum G4OpticalProcessIndex {
  kCerenkov,      ///< Cerenkov process index
  kScintillation, ///< Scintillation process index
  kAbsorption,    ///< Absorption process index
  kRayleigh,      ///< Rayleigh scattering process index
  kMieHG,         ///< Mie scattering process index
  kBoundary,      ///< Boundary process index
  kWLS,           ///< Wave Length Shifting process index
  kNoProcess      ///< Number of processes, no selected process
};

/// Return the name for a given optical process index
G4String G4OpticalProcessName(G4int );

////////////////////
// Inline methods
////////////////////

inline
G4String G4OpticalProcessName(G4int processNumber)
{
  switch ( processNumber ) {
    case kCerenkov:      return "Cerenkov";
    case kScintillation: return "Scintillation";
    case kAbsorption:    return "OpAbsorption";
    case kRayleigh:      return "OpRayleigh";
    case kMieHG:         return "OpMieHG";
    case kBoundary:      return "OpBoundary";
    case kWLS:           return "OpWLS";
    default:             return "NoProcess";
  }
}

#endif // G4OpticalProcessIndex_h
