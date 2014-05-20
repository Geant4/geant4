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
// $Id: G4TouchableDumpScene.hh 66773 2013-01-12 14:48:08Z allison $
//
// 
// John Allison  10th August 1998.
// An artificial scene to find physical volumes.

#ifndef G4TOUCHABLEDUMPSCENE_HH
#define G4TOUCHABLEDUMPSCENE_HH

#include "G4PseudoScene.hh"
#include "G4VisExtent.hh"
#include "G4ModelingParameters.hh"
#include "G4PhysicalVolumeModel.hh"

#include <iostream>

class G4TouchableDumpScene: public G4PseudoScene {

public:

  G4TouchableDumpScene
  (std::ostream& os,
   G4PhysicalVolumeModel* pPVModel,
   const G4ModelingParameters::PVNameCopyNoPath& requiredTouchable);

  virtual ~G4TouchableDumpScene ();

  G4bool IsFound() {return fFound;}

private:

  void ProcessVolume(const G4VSolid&);

  std::ostream& fos;
  const G4PhysicalVolumeModel* fpPVModel;
  G4ModelingParameters::PVNameCopyNoPath fRequiredTouchable;
  G4bool fFound;
};

#endif
