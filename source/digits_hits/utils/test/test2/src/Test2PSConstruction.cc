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

#include "Test2PSConstruction.hh"

#include "G4LogicalVolume.hh"

Test2PSConstruction::Test2PSConstruction(const G4String& name, G4int Segment[3])
  :fname(name)
{ 
  nSegment[0] = Segment[0];
  nSegment[1] = Segment[1];
  nSegment[2] = Segment[2];
}

Test2PSConstruction::~Test2PSConstruction()
{;}


#include "G4SDManager.hh"


#include "G4SDManager.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4PSEnergyDeposit3D.hh"
#include "G4PSTrackLength3D.hh"
#include "G4PSPassageTrackLength3D.hh"
#include "G4PSDoseDeposit3D.hh"
#include "G4PSFlatSurfaceCurrent3D.hh"
#include "G4PSSphereSurfaceCurrent3D.hh"
#include "G4PSPassageCellCurrent3D.hh"
#include "G4PSFlatSurfaceFlux3D.hh"
#include "G4PSCellFlux3D.hh"
#include "G4PSPassageCellFlux3D.hh"
#include "G4PSMinKinEAtGeneration3D.hh"
#include "G4PSNofSecondary3D.hh"
#include "G4PSCellCharge3D.hh"
#include "G4PSNofStep3D.hh"
#include "G4SDParticleFilter.hh"

void Test2PSConstruction::SetupSensitivity(G4LogicalVolume* logVol) {
  //
  // sensitive detectors
  //
  G4SDManager * sdManager = G4SDManager::GetSDMpointer();
  sdManager->SetVerboseLevel(1);

  G4MultiFunctionalDetector* det = new G4MultiFunctionalDetector(fname);

  // filters
  G4String filterName, particleName;
  G4SDParticleFilter* gammaFilter
    = new G4SDParticleFilter(filterName="gammaFilter",particleName="gamma");
  G4SDParticleFilter* electronFilter
    = new G4SDParticleFilter(filterName="electronFilter",particleName="e-");
  G4SDParticleFilter* positronFilter
    = new G4SDParticleFilter(filterName="positronFilter",particleName="e+");

  // primitive scorers
  G4VPrimitiveScorer* primitive;
  primitive = 
    new G4PSEnergyDeposit3D("eDep",nSegment[0],nSegment[1],nSegment[2]);
  det->RegisterPrimitive(primitive);
  primitive = 
    new G4PSTrackLength3D("trackLengthGamma",nSegment[0],nSegment[1],nSegment[2]);
  primitive->SetFilter(gammaFilter);
  det->RegisterPrimitive(primitive);
  primitive = 
    new G4PSTrackLength3D("trackLengthElec",nSegment[0],nSegment[1],nSegment[2]);
  primitive->SetFilter(electronFilter);
  det->RegisterPrimitive(primitive);
  primitive = 
    new G4PSTrackLength3D("trackLengthPosi",nSegment[0],nSegment[1],nSegment[2]);
  primitive->SetFilter(positronFilter);
  det->RegisterPrimitive(primitive);
  primitive = 
    new G4PSNofStep3D("nStepGamma",nSegment[0],nSegment[1],nSegment[2]);
  primitive->SetFilter(gammaFilter);
  det->RegisterPrimitive(primitive);
  primitive = 
    new G4PSNofStep3D("nStepElec",nSegment[0],nSegment[1],nSegment[2]);
  primitive->SetFilter(electronFilter);
  det->RegisterPrimitive(primitive);
  primitive = 
    new G4PSNofStep3D("nStepPosi",nSegment[0],nSegment[1],nSegment[2]);
  primitive->SetFilter(positronFilter);
  det->RegisterPrimitive(primitive);

  //-----
  primitive = new G4PSPassageTrackLength3D("PassageTrackLength",nSegment[0],nSegment[1],nSegment[2]);
  det->RegisterPrimitive(primitive);

  primitive = new G4PSDoseDeposit3D("DoseDeposit",nSegment[0],nSegment[1],nSegment[2]);
  det->RegisterPrimitive(primitive);

  primitive = new G4PSFlatSurfaceCurrent3D("FlatSurfaceCurrent",fCurrent_In);
  det->RegisterPrimitive(primitive);

  primitive = new G4PSPassageCellCurrent3D("PassageCellCurrent",nSegment[0],nSegment[1],nSegment[2]);
  det->RegisterPrimitive(primitive);

  primitive = new G4PSFlatSurfaceFlux3D("FlatSurfaceFlux",fFlux_In);
  det->RegisterPrimitive(primitive);

  primitive = new G4PSCellFlux3D("CellFlux",nSegment[0],nSegment[1],nSegment[2]);
  det->RegisterPrimitive(primitive);

  primitive = new G4PSPassageCellFlux3D("PassageCellFlux",nSegment[0],nSegment[1],nSegment[2]);
  det->RegisterPrimitive(primitive);

  primitive = new G4PSNofSecondary3D("NofSecondary",nSegment[0],nSegment[1],nSegment[2]);
  det->RegisterPrimitive(primitive);

  primitive = new G4PSCellCharge3D("CellCharge",nSegment[0],nSegment[1],nSegment[2]);
  det->RegisterPrimitive(primitive);

  sdManager->AddNewDetector(det);
  logVol->SetSensitiveDetector(det);
  sdManager->SetVerboseLevel(0);

}
