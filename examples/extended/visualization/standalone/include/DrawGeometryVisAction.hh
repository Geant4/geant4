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
/// \file visualization/standalone/include/DrawGeometryVisAction.hh
/// \brief Definition of the DrawGeometryVisAction class
//
//

#ifndef DRAWGEOMETRYVISACTION_HH
#define DRAWGEOMETRYVISACTION_HH

#include "G4VUserVisAction.hh"

namespace B1 {  // Taken from example/basic/B1
  class DetectorConstruction;
}
class G4VPhysicalVolume;
#include "G4Transform3D.hh"
#include "G4VisExtent.hh"

class DrawGeometryVisAction: public G4VUserVisAction {
public:
  DrawGeometryVisAction();
  ~DrawGeometryVisAction();
  void Draw() override;
  const G4VisExtent& GetVisxtent() const {return fExtent;}
private:
  B1::DetectorConstruction* fDetectorConstruction;
  G4VPhysicalVolume* fPhysicalVolume;
  G4Transform3D fTransform;
  G4VisExtent fExtent;
};

#endif

