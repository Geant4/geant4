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
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "Test17RunAction.hh"

#include "G4Run.hh"
#include "G4UImanager.hh"
#include "G4IonC12.hh"
#include "G4ios.hh"
#include "g4std/iomanip"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Test17RunAction::Test17RunAction():
   theProton (G4Proton::Proton()),
   theElectron (G4Electron::Electron())
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Test17RunAction::~Test17RunAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Test17RunAction::BeginOfRunAction(const G4Run* aRun)
{  
  G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;
  edepTot = 0.;
  length  = 0.;
  length2 = 0.;
  xend    = 0.;
  xend2   = 0.;
  nEvents = 0;
  nCharged= 0;
  nNeutral= 0;
  kinEnergy0 = 0.;
  part0   = 0;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Test17RunAction::EndOfRunAction(const G4Run* aRun)
{
  G4double nev = (G4double)nEvents;
  if(nev <= 0.) nev = 1.0;
  edepTot /= nev;
  length  /= nev;
  length2 /= nev;
  xend    /= nev;
  xend2   /= nev;
  G4double sigl = length2 - length*length;
  if(sigl > 0.0) sigl = sqrt(sigl/nev);
  G4double sigx = xend2 - xend*xend;
  if(sigl > 0.0) sigx = sqrt(sigx/nev);
  G4double nc = (G4double)nCharged / nev;
  G4double nn = (G4double)nNeutral / nev;
  G4bool icru = false;

  G4double protR = 0.0;
  G4double protL = 0.0;

  if(abs(kinEnergy0 - 500.0*keV)<0.1*keV && part0 == theProton) { 
    icru  = true;
    protL = 0.009059*mm;
    protR = 0.008869*mm;
  }

  if(abs(kinEnergy0 - MeV)<0.1*keV && part0->GetParticleName() == "pi-") { 
    icru  = true;
    protL = 0.009158*cm*0.9059/0.8869;
    protR = 0.009158*cm;
  }


  G4cout << " ================== run summary =====================" << G4endl;

  G4String name = "";
  if(part0) name = part0->GetParticleName();

  G4cout << G4endl;
  //  G4int prec = G4cout.precision(6);
  G4cout << " end of Run TotNbofEvents = " <<  nEvents 
         << " for " <<  name
         << " with Ekin = " << kinEnergy0/MeV << " MeV" << G4endl ;
  G4cout << "    Track Length in absorber = " <<
          length/micrometer     << " +- " << sigl/micrometer   <<
          "  microns " << G4endl; 
  G4cout << "    CSDA  Range  in absorber = " <<
          xend/micrometer     << " +- " << sigx/micrometer   <<
          "  microns " << G4endl; 
  G4cout << G4endl;
  G4cout << "    Energy deposit in absorber = " <<
           edepTot/MeV << "  MeV " << G4endl ;
  G4cout << G4endl;

  if(icru) {
    G4cout << "### Comparison with the ICRU49 data: " << G4endl;
    G4cout << "    Track Length (G4 - ICRU49) = " 
           << (length - protL)/micrometer
           << " +- " << sigl/micrometer
           << " microns " << G4endl; 
    G4cout << "    CSDA  Range  (G4 - ICRU49) = " 
           << (xend - protR)/micrometer 
           << " +- " << sigx/micrometer
           << " microns " << G4endl;
    G4cout << G4endl ;
  }

  G4cout << "Average Number electrons per event = " << nc << G4endl;
  G4cout << "Average Number photons per event   = " << nn << G4endl;
  G4cout << G4endl; 
      
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....













