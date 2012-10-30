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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "Test17RunAction.hh"

#include "G4SystemOfUnits.hh"
#include "G4Run.hh"
#include "G4UImanager.hh"
#include "G4ios.hh"
#include "G4Proton.hh"
#include "G4Electron.hh"
#include <iomanip>

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
  run     = aRun;
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

void Test17RunAction::EndOfRunAction(const G4Run*)
{
  G4double nev = (G4double)nEvents;
  if(nev <= 0.) nev = 1.0;
  edepTot /= nev;
  length  /= nev;
  length2 /= nev;
  xend    /= nev;
  xend2   /= nev;

  G4double sigl = length2 - length*length;
  if(sigl > 0.0) sigl = std::sqrt(sigl/nev);
  else           sigl = 0.0;

  G4double sigx = xend2 - xend*xend;
  if(sigl > 0.0) sigx = std::sqrt(sigx/nev);
  else           sigx = 0.0;

  G4double nc = (G4double)nCharged / nev;
  G4double nn = (G4double)nNeutral / nev;
  G4bool icru = false;

  G4double protR = 0.0;
  G4double protL = 0.0;

  if(std::abs(kinEnergy0 - 500.0*keV)<0.1*keV && part0 == theProton) { 

    icru  = true;
    protL = 0.009059*mm;
    protR = 0.008869*mm;
  }

  if(std::abs(kinEnergy0 - MeV)<0.1*keV && part0->GetParticleName() == "pi-") { 

    icru  = true;
    protL = 0.009158*cm*0.9059/0.8869;
    protR = 0.009158*cm;
  }


  G4cout << " ================== run summary =====================" << G4endl;

  G4String name = "";
  G4double mass = 0.0;
  if(part0) {
    name = part0->GetParticleName();
    mass = part0->GetPDGMass();
  }
  G4double p = std::sqrt(kinEnergy0*(kinEnergy0 + 2.0*mass))/MeV;
  G4cout << G4endl;
  //  G4int prec = G4cout.precision(6);
  G4cout << " end of Run TotNbofEvents = " <<  nEvents
         << " for " <<  name
         << " with Ekin = " << kinEnergy0/MeV << " MeV" << " p= " << p << " MeV/c" << G4endl;
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













