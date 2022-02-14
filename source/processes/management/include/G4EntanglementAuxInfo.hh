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
// --------------------------------------------------------------------
// GEANT4 class header file
//
// Class Description:
// Manages a clipboard for communicatimg between quantum entangled tracks.
//
// ------------------ G4EntanglementAuxInfo ------------------
//
// Author: J.Allison, May 2017
//
// --------------------------------------------------------------------

// Usage:
//
// This is a concrete class based on G4VAuxiliaryTrackInformation. It is
// designed to hold a shared pointer to a "clip board".  The idea is to
// provide an object - the clip board - that is common across all tracks
// that have this track information and is deleted only when all tracks
// have been tracked.
//
// See G4VEntanglementClipBoard.hh for more infromation.

#ifndef G4EntanglementAuxInfo_hh
#define G4EntanglementAuxInfo_hh

#include "G4VAuxiliaryTrackInformation.hh"

#include <memory>

#include "G4VEntanglementClipBoard.hh"

class G4EntanglementAuxInfo : public G4VAuxiliaryTrackInformation {

public:

  G4EntanglementAuxInfo
  (const std::shared_ptr<G4VEntanglementClipBoard>& clipBoard)
  : fEntanglementClipBoard(clipBoard) {}

  ~G4EntanglementAuxInfo() {}

  G4VEntanglementClipBoard* GetEntanglementClipBoard() const
  {return fEntanglementClipBoard.get();}

private:

  std::shared_ptr<G4VEntanglementClipBoard> fEntanglementClipBoard;

};

#endif
