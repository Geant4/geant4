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
//---------------------------------------------------------------------------
//
// ClassName:   hTestHisto
//  
//
// Author:      V.Ivanchenko 30/01/01
//
//----------------------------------------------------------------------------
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "hTestHisto.hh"
#include "g4std/iomanip"

#include "CLHEP/Hist/HBookFile.h"
#include <assert.h>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

hTestHisto* hTestHisto::fManager = 0;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

hTestHisto* hTestHisto::GetPointer()
{
  if(!fManager) {
    fManager = new hTestHisto();
  }
  return fManager;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

hTestHisto::hTestHisto() 
{
  verbose = 0;
  histName = G4String("histo.hbook");
  hbookManager = 0;
  ntup = 0;
  nHisto = 1;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

hTestHisto::~hTestHisto() 
{
  delete hbookManager;
  histo.clear(); 
  G4cout << "hTestHisto: Histograms are deleted for " << theName << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void hTestHisto::BeginOfHisto(G4int num)
{  
  if(0 < verbose) G4cout << "hTestHisto # " << num << " started " << G4endl;
  zend     = 0.0;
  zend2    = 0.0;
  zEvt     = 0.0;
  
  if(0 < nHisto) bookHisto();

  if(verbose > 0) {
    G4cout << "hTestHisto: Histograms are booked and run has been started" 
           << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void hTestHisto::EndOfHisto()
{

  G4cout << "hTestHisto: End of run actions are started" << G4endl;

  // Zend average
  if(zEvt > 0.0) {
    zend  /= zEvt;
    zend2 /= zEvt;
    G4double sig = sqrt(zend2 - zend*zend);
    zend2 = sig / sqrt(zEvt);
    G4cout<<"========================================================"<<G4endl;
    G4cout << setprecision(4) << "Range(mm)= " << zend/mm 
           << "; Stragling(mm)= " << sig/mm 
           << setprecision(2) << " +- " << zend2/mm << G4endl;
    G4cout<<"========================================================"<<G4endl;
  }  

   // Write histogram file
  if(0 < nHisto) {
    hbookManager->write();
    G4cout << "Histograms and Ntuples are saved" << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void hTestHisto::SaveEvent()
{
  if(0 < nHisto) ntup->dumpData();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void hTestHisto::SaveToTuple(const G4String& parname, G4double val)
{
  if(0 < nHisto) ntup->column(parname,val);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void hTestHisto::SaveToTuple(const G4String& parname,G4double val,G4double defval)
{
  if(0 < nHisto) ntup->column(parname,val,defval);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void hTestHisto::bookHisto()
{
  G4cout << "hTestHisto: Histograms will be saved to the file <" 
         << histName << ">" << G4endl;

  // init hbook
  hbookManager = new HBookFile(histName, 68);

  // book histograms

  histo.resize(nHisto);

  G4double zmax = AbsorberThickness * (G4double)NumberOfAbsorbers / mm;

  if(0 < nHisto) histo[0] = hbookManager->histogram(
    "Energy deposit (MeV) in absorber (mm)",NumberOfAbsorbers,0.0,zmax);

  if(1 < nHisto) histo[1] = hbookManager->histogram(
    "Energy (MeV) of delta-electrons",50,0.0,maxEnergy/MeV);

  if(2 < nHisto) histo[2] = hbookManager->histogram(
    "Theta (degrees) of delta-electrons",36,0.0,180.);

  if(3 < nHisto) histo[3] = hbookManager->histogram(
    "Energy (MeV) of secondary gamma",50,0.0,maxEnergy/MeV);

  if(4 < nHisto) histo[4] = hbookManager->histogram(
    "Theta (degrees) of secondary gamma",36,0.0,180.);

  // book ntuple
  ntup = hbookManager->ntuple("Range/Energy");

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void hTestHisto::AddEnergy(G4double edep, G4double z)
{
  if(0 < nHisto) histo[0]->accumulate(z/mm, edep/MeV);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void hTestHisto::AddEndPoint(G4double z)
{
  zend  += z;
  zend2 += z*z;
  zEvt  += 1.0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void hTestHisto::AddDeltaElectron(const G4DynamicParticle* elec)
{
  if(1 < nHisto) histo[1]->accumulate(elec->GetKineticEnergy()/MeV,1.0);
  if(2 < nHisto)
     histo[2]->accumulate((elec->GetMomentumDirection()).theta()/deg,1.0);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void hTestHisto::AddPhoton(const G4DynamicParticle* ph)
{
  if(3 < nHisto) histo[3]->accumulate(ph->GetKineticEnergy()/MeV,1.0);
  if(4 < nHisto)
     histo[4]->accumulate((ph->GetMomentumDirection()).theta()/deg,1.0);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

