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
// -------------------------------------------------------------
//
//
//      ---------- test31SD -------------
//              
//  Modified: 05.04.01 Vladimir Ivanchenko new design of test31 
// 
// -------------------------------------------------------------

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "test31SD.hh"

#include "globals.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4RunManager.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4Track.hh"
#include "G4Positron.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

test31SD::test31SD(G4String name)
 :G4VSensitiveDetector(name),
  theHisto(test31Histo::GetPointer()),
  evno(0),
  evnOld(-1),
  trIDold(-1),
  delta(1.0e-6*mm)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

test31SD::~test31SD()
{
  energy.clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void test31SD::Initialize(G4HCofThisEvent*)
{
  evno++;
  numAbs = theHisto->GetNumberOfAbsorbers();
  if(0 < theHisto->GetVerbose()) 
    G4cout << "test31SD: Begin Of Event # " 
           << evno
           << "  numAbs= " << numAbs 
           << G4endl;
 
  energy.resize(numAbs);
  for(G4int i=0; i<numAbs; i++) { energy[i] = 0.0; }
  backEnergy = 0.0;
  leakEnergy = 0.0;
  G4double gap = theHisto->GetGap();
  zmax = (theHisto->GetAbsorberThickness() + gap)*numAbs - gap - delta;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool test31SD::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{
  //G4cout << "SD: =====" << G4endl;
  G4double edep = aStep->GetTotalEnergyDeposit();

  theHisto->AddStep();
  theHisto->AddTrackLength(aStep->GetStepLength());

  G4int j = aStep->GetPreStepPoint()->GetPhysicalVolume()->GetCopyNo();

  //G4cout << "SD: edep= " << edep << "  j= " << j << G4endl;
   
  if(j < 0 || j >= numAbs) {
      G4cout << "Warning!!! test31SD: cannot add " << edep/MeV
             << " MeV to the slice # " << j << G4endl;

  } else {
      energy[j] += edep;
  }

  const G4Track* track = aStep->GetTrack();
  G4int trIDnow  = track->GetTrackID();

  if(1 < theHisto->GetVerbose() && edep > 0.0) {
      G4cout << "test31SD: energy = " << edep/MeV
             << " MeV is deposited at " << j
             << "-th absorber slice "
             << " TrackID= " << trIDnow 
             << G4endl;
  }

  G4double tkin  = track->GetKineticEnergy(); 
  G4double theta = (track->GetMomentumDirection()).theta();
  G4double zend  = aStep->GetPostStepPoint()->GetPosition().z();

  G4bool stop = false;
  G4bool primary = false;
  G4bool outAbs = false;

  //G4cout << "Next vol " << track->GetNextVolume() << G4endl;

  if(track->GetNextVolume() != track->GetVolume()
     && (zend <= delta || zend >= zmax)  ) { outAbs = true; }

  if(tkin == 0.0) { stop = true; }
  if(0 == aStep->GetTrack()->GetParentID()) { primary = true; }

  // new particle
  if(trIDnow != trIDold || evno != evnOld) {
    trIDold = trIDnow;
    evnOld  = evno;
    tkinold = aStep->GetPreStepPoint()->GetKineticEnergy();
    part_is_out = true;
  }

  // Primary particle stop 
  //G4cout << "primary= " << primary << " stop= " << stop 
  //<< " outAbs= " << outAbs << G4endl;
   
  if(primary && (stop || outAbs)) {

      G4double xend = aStep->GetPostStepPoint()->GetPosition().x();
      G4double yend = aStep->GetPostStepPoint()->GetPosition().y();
      theHisto->SaveToTuple("xend",xend/mm);      
      theHisto->SaveToTuple("yend",yend/mm);      
      theHisto->SaveToTuple("zend",zend/mm);      
      theHisto->SaveToTuple("ltpk",(theHisto->GetTrackLength())/mm);      
      theHisto->SaveToTuple("tend",tkin);
      theHisto->SaveToTuple("teta",theta);      
      theHisto->SaveToTuple("loss",(tkinold-tkin)/MeV);      
      theHisto->SaveToTuple("dedx",(tkinold-tkin)*mm/(zmax*MeV));      

    // exclude cases when track return back

      if(theHisto->GetTrackLength() < 2.0*zend) theHisto->AddEndPoint(zend);

  }

  // After step in absorber

  if(outAbs && part_is_out) {
    G4double e = tkin;
    
    if(track->GetDefinition() == G4Positron::Positron()) 
      e += 2.*electron_mass_c2;

    if(zend > zmax-delta)    {
       leakEnergy += e;
       theHisto->AddParticleLeak(track);       
    } else if(zend < delta)  {
       backEnergy += e;
       theHisto->AddParticleBack(track);       
    }
    part_is_out = false;
  }
  //G4cout << " ================ " << G4endl;
  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void test31SD::EndOfEvent(G4HCofThisEvent*)
{
  G4int i, j;

  theHisto->SaveToTuple(G4String("back"),backEnergy);      
  theHisto->SaveToTuple(G4String("leak"),leakEnergy);      
  G4double etot = 0.0;

  // The histogramm on the energy deposition profile
  if(numAbs > 0) {
    G4double sss = theHisto->GetAbsorberThickness();
    G4double z = -0.5 * sss;
    for(i=0; i<numAbs; i++) {
      z += sss; 
      etot += energy[i];
      //      G4cout << "i= " << i << " e= " << energy[i] << G4endl;
      theHisto->AddEnergy(energy[i], z);
    }
  }
  theHisto->SaveToTuple(G4String("edep"),etot);      
  //theHisto->SaveToTuple(G4String("Evt"),G4double(evno));      

  // Integrated energy deposition to nTuple
  G4int nMax = 60;
  G4double EE[60];
  /*
  G4String eSlice[60]={
      "S0", "S1", "S2", "S3", "S4", "S5", "S6", "S7", "S8", "S9", 
      "S10", "S11", "S12", "S13", "S14", "S15", "S16", "S17", "S18", "S19", 
      "S20", "S21", "S22", "S23", "S24", "S25", "S26", "S27", "S28", "S29", 
      "S30", "S31", "S32", "S33", "S34", "S35", "S36", "S37", "S38", "S39", 
      "S40", "S41", "S42", "S43", "S44", "S45", "S46", "S47", "S48", "S49", 
      "S50", "S51", "S52", "S53", "S54", "S55", "S56", "S57", "S58", "S59"};
  G4String eInteg[60]={
      "E0", "E1", "E2", "E3", "E4", "E5", "E6", "E7", "E8", "E9", 
      "E10", "E11", "E12", "E13", "E14", "E15", "E16", "E17", "E18", "E19", 
      "E20", "E21", "E22", "E23", "E24", "E25", "E26", "E27", "E28", "E29", 
      "E30", "E31", "E32", "E33", "E34", "E35", "E36", "E37", "E38", "E39", 
      "E40", "E41", "E42", "E43", "E44", "E45", "E46", "E47", "E48", "E49", 
      "E50", "E51", "E52", "E53", "E54", "E55", "E56", "E57", "E58", "E59"};
  */
  G4int k = theHisto->GetNumAbsorbersSaved();
  if (nMax > k) nMax = k;

  if(nMax > 1) {
    for(i=0; i<nMax; i++){
      EE[i]=0.0;
      for(j=0; j<i+1; j++) {
        EE[i] += energy[j];
      }  
      //      theHisto->SaveToTuple(eSlice[i],energy[i]);      
      //      theHisto->SaveToTuple(eInteg[i],EE[i]);      
    }
  }

  // Dump information about this event 
  theHisto->SaveEvent();

  if(0 < theHisto->GetVerbose()) { 
    G4cout << "test31SD: Event #" << evno << " ended" << G4endl;
    G4cout << "Edep(MeV)= " << etot/MeV 
           << "; backE(MeV)= " << backEnergy/MeV
           << "; leakE(MeV)= " << leakEnergy/MeV
           << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void test31SD::clear()
{} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


void test31SD::PrintAll()
{} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....







