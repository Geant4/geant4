//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
// -------------------------------------------------------------
//
//
// -------------------------------------------------------------
//      GEANT4 test31
//
//      History: based on object model of
//      2nd December 1995, G.Cosmo
//      ---------- test31SD -------------
//              
//  Modified: 05.04.01 Vladimir Ivanchenko new design of test31 
// 
// -------------------------------------------------------------

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "test31SD.hh"

#include "G4RunManager.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4Track.hh"
#include "globals.hh"

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
  if(0 < theHisto->GetVerbose()) 
    G4cout << "test31SD: Begin Of Event # " << evno << G4endl;

  numAbs = theHisto->GetNumberOfAbsorbers();
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
  G4double edep = aStep->GetTotalEnergyDeposit();
  //  if(0.0 == edep) return true;

  theHisto->AddStep();
  theHisto->AddTrackLength(aStep->GetStepLength());

  G4int j = aStep->GetPreStepPoint()->GetPhysicalVolume()->GetCopyNo();
   
  if(j < 0 || j >= numAbs) {
      G4cout << "Warning!!! test31SD: cannot add " << edep/MeV
             << " MeV to the slice # " << j << G4endl;

  } else {
      energy[j] += edep;
  }

  if(1 < theHisto->GetVerbose()) {
      G4cout << "test31SD: energy = " << edep/MeV
             << " MeV is deposited at " << j
             << "-th absorber slice " << G4endl;
  }

  const G4Track* track = aStep->GetTrack();
  G4int trIDnow  = track->GetTrackID();
  G4double tkin  = track->GetKineticEnergy(); 
  G4double theta = (track->GetMomentumDirection()).theta();
  G4double zend  = aStep->GetPostStepPoint()->GetPosition().z();
  //  G4double zstart= aStep->GetPreStepPoint()->GetPosition().z();

  G4bool stop = false;
  G4bool primary = false;
  G4bool outAbs = false;
  if(track->GetNextVolume()->GetName() == "World"
     && (zend <= delta || zend >= zmax)  ) outAbs = true;

  if(tkin == 0.0) stop = true;
  if(0 == aStep->GetTrack()->GetParentID()) primary = true;

  // new particle
  if(trIDnow != trIDold || evno != evnOld) {
    trIDold = trIDnow;
    evnOld  = evno;
    tkinold = aStep->GetPreStepPoint()->GetKineticEnergy();
  }

  // Primary particle stop 
   
  if(primary && (stop || outAbs)) {

      G4double xend = aStep->GetPostStepPoint()->GetPosition().x();
      G4double yend = aStep->GetPostStepPoint()->GetPosition().y();
      theHisto->SaveToTuple("XEND",xend/mm);      
      theHisto->SaveToTuple("YEND",yend/mm);      
      theHisto->SaveToTuple("ZEND",zend/mm);      
      theHisto->SaveToTuple("LTRK",(theHisto->GetTrackLength())/mm);      
      theHisto->SaveToTuple("TEND",tkin);
      theHisto->SaveToTuple("TETA",theta);      
      theHisto->SaveToTuple("LOSS",(tkinold-tkin)/MeV);      
      theHisto->SaveToTuple("DEDX",(tkinold-tkin)*mm/(zmax*MeV));      

    // exclude cases when track return back

      if(theHisto->GetTrackLength() < 2.0*zend) theHisto->AddEndPoint(zend);

  }

  // After step in absorber

  if(outAbs) {
      if(zend >= zmax)       leakEnergy += tkin;
      else                   backEnergy += tkin;
  }
  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void test31SD::EndOfEvent(G4HCofThisEvent*)
{
  G4int i, j;

  theHisto->SaveToTuple(G4String("backE"),backEnergy);      
  theHisto->SaveToTuple(G4String("leakE"),leakEnergy);      
  G4double etot = 0.0;

  // The histogramm on the energy deposition profile
  if(numAbs > 0) {
    G4double s = theHisto->GetAbsorberThickness();
    G4double z = -0.5 * s;
    for(i=0; i<numAbs; i++) {
      z += s; 
      etot += energy[i];
      theHisto->AddEnergy(energy[i], z);
    }
  }
  theHisto->SaveToTuple(G4String("Edep"),etot);      
  theHisto->SaveToTuple(G4String("Evt"),G4double(evno));      

  // Integrated energy deposition to nTuple
  G4int nMax = 60;
  G4double EE[60];
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

  G4int k = theHisto->GetNumAbsorbersSaved();
  if (nMax > k) nMax = k;

  if(nMax > 1) {
    for(i=0; i<nMax; i++){
      EE[i]=0.0;
      for(j=0; j<i+1; j++) {
        EE[i] += energy[j];
      }  
      theHisto->SaveToTuple(eSlice[i],energy[i]);      
      theHisto->SaveToTuple(eInteg[i],EE[i]);      
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







