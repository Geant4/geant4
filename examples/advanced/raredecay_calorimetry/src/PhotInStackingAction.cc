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

//#define debug
//#define errors

#include "PhotInStackingAction.hh"

// Can be moved to the PhotInConstants.hh
G4int PhotInStackingAction::PhotInNGamma     [PhotInDiNSections] = {0,0,0,0,0,0};
G4int PhotInStackingAction::PhotInNElectron  [PhotInDiNSections] = {0,0,0,0,0,0};
G4int PhotInStackingAction::PhotInNPositron  [PhotInDiNSections] = {0,0,0,0,0,0};
G4int PhotInStackingAction::PhotInNNeutrons  [PhotInDiNSections] = {0,0,0,0,0,0};
G4int PhotInStackingAction::PhotInNProtons   [PhotInDiNSections] = {0,0,0,0,0,0};
G4int PhotInStackingAction::PhotInNDeuterons [PhotInDiNSections] = {0,0,0,0,0,0};
G4int PhotInStackingAction::PhotInNTritons   [PhotInDiNSections] = {0,0,0,0,0,0};
G4int PhotInStackingAction::PhotInNHe3s      [PhotInDiNSections] = {0,0,0,0,0,0};
G4int PhotInStackingAction::PhotInNAlphas    [PhotInDiNSections] = {0,0,0,0,0,0};
G4int PhotInStackingAction::PhotInNLambdas   [PhotInDiNSections] = {0,0,0,0,0,0};
G4int PhotInStackingAction::PhotInNHeavyFrags[PhotInDiNSections] = {0,0,0,0,0,0};
G4int PhotInStackingAction::PhotInNMesons    [PhotInDiNSections] = {0,0,0,0,0,0};
G4double PhotInStackingAction::PhotInEMinGamma     [PhotInDiNSections] = {DBL_MAX,DBL_MAX,DBL_MAX,DBL_MAX,DBL_MAX,DBL_MAX};
G4double PhotInStackingAction::PhotInEMinElectron  [PhotInDiNSections] = {DBL_MAX,DBL_MAX,DBL_MAX,DBL_MAX,DBL_MAX,DBL_MAX};
G4double PhotInStackingAction::PhotInEMinPositron  [PhotInDiNSections] = {DBL_MAX,DBL_MAX,DBL_MAX,DBL_MAX,DBL_MAX,DBL_MAX};
G4double PhotInStackingAction::PhotInEMinNeutrons  [PhotInDiNSections] = {DBL_MAX,DBL_MAX,DBL_MAX,DBL_MAX,DBL_MAX,DBL_MAX};
G4double PhotInStackingAction::PhotInEMinProtons   [PhotInDiNSections] = {DBL_MAX,DBL_MAX,DBL_MAX,DBL_MAX,DBL_MAX,DBL_MAX};
G4double PhotInStackingAction::PhotInEMinDeuterons [PhotInDiNSections] = {DBL_MAX,DBL_MAX,DBL_MAX,DBL_MAX,DBL_MAX,DBL_MAX};
G4double PhotInStackingAction::PhotInEMinTritons   [PhotInDiNSections] = {DBL_MAX,DBL_MAX,DBL_MAX,DBL_MAX,DBL_MAX,DBL_MAX};
G4double PhotInStackingAction::PhotInEMinHe3s      [PhotInDiNSections] = {DBL_MAX,DBL_MAX,DBL_MAX,DBL_MAX,DBL_MAX,DBL_MAX};
G4double PhotInStackingAction::PhotInEMinAlphas    [PhotInDiNSections] = {DBL_MAX,DBL_MAX,DBL_MAX,DBL_MAX,DBL_MAX,DBL_MAX};
G4double PhotInStackingAction::PhotInEMinLambdas   [PhotInDiNSections] = {DBL_MAX,DBL_MAX,DBL_MAX,DBL_MAX,DBL_MAX,DBL_MAX};
G4double PhotInStackingAction::PhotInEMinHeavyFrags[PhotInDiNSections] = {DBL_MAX,DBL_MAX,DBL_MAX,DBL_MAX,DBL_MAX,DBL_MAX};
G4double PhotInStackingAction::PhotInEMinMesons    [PhotInDiNSections] = {DBL_MAX,DBL_MAX,DBL_MAX,DBL_MAX,DBL_MAX,DBL_MAX};
G4double PhotInStackingAction::PhotInEMaxGamma     [PhotInDiNSections]={0.,0.,0.,0.,0.,0.};
G4double PhotInStackingAction::PhotInEMaxElectron  [PhotInDiNSections]={0.,0.,0.,0.,0.,0.};
G4double PhotInStackingAction::PhotInEMaxPositron  [PhotInDiNSections]={0.,0.,0.,0.,0.,0.};
G4double PhotInStackingAction::PhotInEMaxNeutrons  [PhotInDiNSections]={0.,0.,0.,0.,0.,0.};
G4double PhotInStackingAction::PhotInEMaxProtons   [PhotInDiNSections]={0.,0.,0.,0.,0.,0.};
G4double PhotInStackingAction::PhotInEMaxDeuterons [PhotInDiNSections]={0.,0.,0.,0.,0.,0.};
G4double PhotInStackingAction::PhotInEMaxTritons   [PhotInDiNSections]={0.,0.,0.,0.,0.,0.};
G4double PhotInStackingAction::PhotInEMaxHe3s      [PhotInDiNSections]={0.,0.,0.,0.,0.,0.};
G4double PhotInStackingAction::PhotInEMaxAlphas    [PhotInDiNSections]={0.,0.,0.,0.,0.,0.};
G4double PhotInStackingAction::PhotInEMaxLambdas   [PhotInDiNSections]={0.,0.,0.,0.,0.,0.};
G4double PhotInStackingAction::PhotInEMaxHeavyFrags[PhotInDiNSections]={0.,0.,0.,0.,0.,0.};
G4double PhotInStackingAction::PhotInEMaxMesons    [PhotInDiNSections]={0.,0.,0.,0.,0.,0.};

PhotInStackingAction::PhotInStackingAction() {}

PhotInStackingAction::~PhotInStackingAction(){}

G4ClassificationOfNewTrack PhotInStackingAction::ClassifyNewTrack(const G4Track* aTrack)
{
  G4TouchableHistory* theTouchable = (G4TouchableHistory*)(aTrack->GetTouchable());
  G4int depth=0;                  // Depth of geometrical hierarchy of volumes
  if(theTouchable)depth=theTouchable->GetHistoryDepth();
#ifdef debug
  G4cout<<"PhotInStackingAction::ClassifyNewTrack: --------->>>>>>> depth="<<depth<<G4endl;
#endif
#ifdef debug
		G4cout<<"PhotInStackingAction::ClassifyNewTrack: XYZ="<<aTrack->GetPosition()<<G4endl;
#endif
  if(depth>1)                     // It is in one of the Sections, not in the Mother volume
  {
    G4VPhysicalVolume* calPhys = theTouchable->GetVolume(depth-1);// PhysVolume of Sections
    if(calPhys)
    {
      G4String sectName=calPhys->GetName();
#ifdef debug
      G4cout<<"PhotInStackingAction::ClassifyNewTrack: SectionName="<<sectName<<G4endl;
#endif
      G4int i = -1;                // 2 * Section_number
      if     (sectName==PhotInCalName[0]) { i = 0; }
      else if(sectName==PhotInCalName[1]) { i = 2; }
      else if(sectName==PhotInCalName[2]) { i = 4; }
#ifdef errors
      if(i<0) G4cerr<<"*PhotInStackingAction::ClassifyNewTrack: sensitive i="<<i<<G4endl;
#endif
#ifdef debug
      G4int layerNumber=theTouchable->GetReplicaNumber(depth-2);//WhatIfLayers notFill Sect
      G4cout<<"PhotInStackingAction::ClassifyNewTrack: Layer#"<<layerNumber<<G4endl;
#endif
      G4int slabNumber=-1;        // Prototype of the slab number (-1 means absorber)
      if(depth>2) slabNumber=theTouchable->GetReplicaNumber(depth-3);
      else i++;                   // Collect information in the absorber
#ifdef debug
      G4cout<<"PhotInStackingAction::ClassifyNewTrack:Slab#"<<slabNumber<<",i="<<i<<G4endl;
#endif
      G4ParticleDefinition* particleType = aTrack->GetDefinition();
      if(particleType==G4Gamma::GammaDefinition())
      {
        PhotInNGamma[i]++;
        if(PhotInEMinGamma[i]>aTrack->GetKineticEnergy())
          PhotInEMinGamma[i]=aTrack->GetKineticEnergy();
        if(PhotInEMaxGamma[i]<aTrack->GetKineticEnergy())
          PhotInEMaxGamma[i]=aTrack->GetKineticEnergy();
      }
      else if(particleType==G4Electron::ElectronDefinition())
      {
        PhotInNElectron[i]++;
        if(PhotInEMinElectron[i]>aTrack->GetKineticEnergy())
          PhotInEMinElectron[i]=aTrack->GetKineticEnergy();
        if(PhotInEMaxElectron[i]<aTrack->GetKineticEnergy())
          PhotInEMaxElectron[i]=aTrack->GetKineticEnergy();
      }
      else if(particleType==G4Positron::PositronDefinition())
      {
        PhotInNPositron[i]++;
        if(PhotInEMinPositron[i]>aTrack->GetKineticEnergy())
          PhotInEMinPositron[i]=aTrack->GetKineticEnergy();
        if(PhotInEMaxPositron[i]<aTrack->GetKineticEnergy())
          PhotInEMaxPositron[i]=aTrack->GetKineticEnergy();
      }
      else if(particleType==G4Neutron::NeutronDefinition())
      {
        PhotInNNeutrons[i]++;
        if(PhotInEMinNeutrons[i]>aTrack->GetKineticEnergy())
          PhotInEMinNeutrons[i]=aTrack->GetKineticEnergy();
        if(PhotInEMaxNeutrons[i]<aTrack->GetKineticEnergy())
          PhotInEMaxNeutrons[i]=aTrack->GetKineticEnergy();
      }
      else if(particleType==G4Proton::ProtonDefinition())
      {
        PhotInNProtons[i]++;
        if(PhotInEMinProtons[i]>aTrack->GetKineticEnergy())
          PhotInEMinProtons[i]=aTrack->GetKineticEnergy();
        if(PhotInEMaxProtons[i]<aTrack->GetKineticEnergy())
          PhotInEMaxProtons[i]=aTrack->GetKineticEnergy();
      }
      else if(particleType==G4Deuteron::DeuteronDefinition())
      {
        PhotInNDeuterons[i]++;
        if(PhotInEMinDeuterons[i]>aTrack->GetKineticEnergy())
          PhotInEMinDeuterons[i]=aTrack->GetKineticEnergy();
        if(PhotInEMaxDeuterons[i]<aTrack->GetKineticEnergy())
          PhotInEMaxDeuterons[i]=aTrack->GetKineticEnergy();
      }
      else if(particleType==G4Triton::TritonDefinition())
      {
        PhotInNTritons[i]++;
        if(PhotInEMinTritons[i]>aTrack->GetKineticEnergy())
          PhotInEMinTritons[i]=aTrack->GetKineticEnergy();
        if(PhotInEMaxTritons[i]<aTrack->GetKineticEnergy())
          PhotInEMaxTritons[i]=aTrack->GetKineticEnergy();
      }
      else if(particleType==G4He3::He3Definition())
      {
        PhotInNHe3s[i]++;
        if(PhotInEMinHe3s[i]>aTrack->GetKineticEnergy())
          PhotInEMinHe3s[i]=aTrack->GetKineticEnergy();
        if(PhotInEMaxHe3s[i]<aTrack->GetKineticEnergy())
          PhotInEMaxHe3s[i]=aTrack->GetKineticEnergy();
      }
      else if(particleType==G4Alpha::AlphaDefinition())
      {
        PhotInNAlphas[i]++;
        if(PhotInEMinAlphas[i]>aTrack->GetKineticEnergy())
          PhotInEMinAlphas[i]=aTrack->GetKineticEnergy();
        if(PhotInEMaxAlphas[i]<aTrack->GetKineticEnergy())
          PhotInEMaxAlphas[i]=aTrack->GetKineticEnergy();
      }
      else if(particleType==G4Lambda::LambdaDefinition())
      {
        PhotInNLambdas[i]++;
        if(PhotInEMinLambdas[i]>aTrack->GetKineticEnergy())
          PhotInEMinLambdas[i]=aTrack->GetKineticEnergy();
        if(PhotInEMaxLambdas[i]<aTrack->GetKineticEnergy())
          PhotInEMaxLambdas[i]=aTrack->GetKineticEnergy();
      }
      else if(particleType->GetBaryonNumber()>4)
      {
        PhotInNHeavyFrags[i]++;
        if(PhotInEMinHeavyFrags[i]>aTrack->GetKineticEnergy())
          PhotInEMinHeavyFrags[i]=aTrack->GetKineticEnergy();
        if(PhotInEMaxHeavyFrags[i]<aTrack->GetKineticEnergy())
          PhotInEMaxHeavyFrags[i]=aTrack->GetKineticEnergy();
      }
      else if(!particleType->GetBaryonNumber())
      {
        PhotInNMesons[i]++;
        if(PhotInEMinMesons[i]>aTrack->GetKineticEnergy())
          PhotInEMinMesons[i]=aTrack->GetKineticEnergy();
        if(PhotInEMaxMesons[i]<aTrack->GetKineticEnergy())
          PhotInEMaxMesons[i]=aTrack->GetKineticEnergy();
      }
    }
  }
  return fUrgent;
}

void PhotInStackingAction::PrepareNewEvent()
{ 
  for(int i=0; i<PhotInDiNSections ;i++)
  {
    PhotInNGamma     [i] = 0; 
    PhotInNElectron  [i] = 0;
    PhotInNPositron  [i] = 0;
    PhotInNNeutrons  [i] = 0; 
    PhotInNProtons   [i] = 0;
    PhotInNDeuterons [i] = 0;
    PhotInNTritons   [i] = 0; 
    PhotInNHe3s      [i] = 0;
    PhotInNAlphas    [i] = 0;
    PhotInNLambdas   [i] = 0; 
    PhotInNHeavyFrags[i] = 0;
    PhotInNMesons    [i] = 0;

    PhotInEMinGamma     [i] = DBL_MAX;
    PhotInEMinElectron  [i] = DBL_MAX;
    PhotInEMinPositron  [i] = DBL_MAX;
    PhotInEMinNeutrons  [i] = DBL_MAX;
    PhotInEMinProtons   [i] = DBL_MAX;
    PhotInEMinDeuterons [i] = DBL_MAX;
    PhotInEMinTritons   [i] = DBL_MAX;
    PhotInEMinHe3s      [i] = DBL_MAX;
    PhotInEMinAlphas    [i] = DBL_MAX;
    PhotInEMinLambdas   [i] = DBL_MAX;
    PhotInEMinHeavyFrags[i] = DBL_MAX;
    PhotInEMinMesons    [i] = DBL_MAX;

    PhotInEMaxGamma     [i] = 0.;
    PhotInEMaxElectron  [i] = 0.;
    PhotInEMaxPositron  [i] = 0.;
    PhotInEMaxNeutrons  [i] = 0.;
    PhotInEMaxProtons   [i] = 0.;
    PhotInEMaxDeuterons [i] = 0.;
    PhotInEMaxTritons   [i] = 0.;
    PhotInEMaxHe3s      [i] = 0.;
    PhotInEMaxAlphas    [i] = 0.;
    PhotInEMaxLambdas   [i] = 0.;
    PhotInEMaxHeavyFrags[i] = 0.;
    PhotInEMaxMesons    [i] = 0.;
  }
}


