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
//
// John Allison  22nd August 2018.
// An artificial scene to find the properties of a touchable.

#ifndef G4TOUCHABLEPROPERTIESSCENE_HH
#define G4TOUCHABLEPROPERTIESSCENE_HH

#include "G4PseudoScene.hh"
#include "G4ModelingParameters.hh"
#include "G4Transform3D.hh"
#include "G4PhysicalVolumeModel.hh"

#include <iostream>

class G4TouchablePropertiesScene: public G4PseudoScene {

public:

  G4TouchablePropertiesScene
  (G4PhysicalVolumeModel* pSearchPVModel,
   const G4ModelingParameters::PVNameCopyNoPath& requiredTouchable);

  virtual ~G4TouchablePropertiesScene ();

  const G4PhysicalVolumeModel::TouchableProperties&
  GetFoundTouchableProperties() const
  {return fFoundTouchableProperties;}

private:

  void ProcessVolume(const G4VSolid&);

  const G4PhysicalVolumeModel*           fpSearchPVModel;
  G4ModelingParameters::PVNameCopyNoPath fRequiredTouchable;
  G4PhysicalVolumeModel::TouchableProperties fFoundTouchableProperties;
};

#endif
