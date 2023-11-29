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
// Authors: M. Large, S. Guatelli -  University of Wollongong, Australia
// This class developed by John Allison, March 2021
//

#ifndef ICRP110PhantomPseudoScene_h
#define ICRP110PhantomPseudoScene_h 1

#include "G4PseudoScene.hh"

#include "G4Colour.hh"
#include "G4ThreeVector.hh"

#include <map>

class G4PhysicalVolumeModel;
class G4Material;

class ICRP110PhantomPseudoScene: public G4PseudoScene
{
public:
  // Given a world volume, fill a multimap provided by the instantiator
  ICRP110PhantomPseudoScene
  (G4PhysicalVolumeModel*  // input...the following are outputs
   ,std::multimap<const G4Material*,G4ThreeVector>& positionByMaterial
   ,std::map     <const G4Material*,G4Colour>&      colourByMaterial);

  ~ICRP110PhantomPseudoScene();

  using G4PseudoScene::AddSolid;  // except for...
  void AddSolid(const G4Box&);
  void ProcessVolume(const G4VSolid&);

  const std::multimap<const G4Material*,G4ThreeVector>& GetPositionByMaterial()
  {return fPositionByMaterial;}

  const std::map     <const G4Material*,G4Colour>&      GetColourByMaterial()
  {return fColourByMaterial;}

private:
  G4PhysicalVolumeModel* fpPVModel;
  std::multimap<const G4Material*,G4ThreeVector>& fPositionByMaterial;
  std::map     <const G4Material*,G4Colour>&      fColourByMaterial;
};
#endif

