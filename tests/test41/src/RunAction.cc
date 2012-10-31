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
#include <fstream>
#include <iomanip>
#include <vector>

#include "RunAction.hh"
#include "globals.hh"
#include "G4ios.hh"
#include "G4SystemOfUnits.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction()
{ 
  ptot = 172*MeV;
  sigma = 0;
  fname   = "test41.out";
  matname = "unknown";
  const G4double t[11] = {0.00269, 0.00895, 0.0162, 0.0248, 0.0347, 0.0463,
                          0.0597, 0.0754, 0.0938, 0.1151, 3.15};
  for(G4int i=0; i<11; i++) {tbins[i] = t[i];}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run* aRun)
{
  G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;
  nevt = 0;
  for(G4int i=0; i<11; i++) {hist[i] = 0;}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run* aRun)
{
  if(nevt == 0) {
    G4cout << "### Run " << aRun->GetRunID() << " ended with zero events" 
	   << G4endl;
    return;
  }
  G4double x[11];
  G4double dx[11];
  G4double norm = 1.0/G4double(nevt);
  G4double w, y, z;
  G4int i;
  for(i=0; i<11; i++) {
    if(i == 0) w = tbins[0];
    else       w = tbins[i] - tbins[i - 1];
    y = G4double(hist[i]);
    x[i] = 0.5*y*norm/w;
    z    = std::sqrt(y);
    if(z <= 1.0) z = 1.0;   
    dx[i]= x[i]/z;
  }

  std::ofstream fout(fname, std::ios::out|std::ios::trunc);
  if(!fout.is_open()) {
    G4cout << "### Error in open of file " << fname  
	   << G4endl;
    return;
  }    
  G4cout << "###  " << nevt << "  events of mu+  with p(MeV/c)= "
	 << ptot/MeV << " scattered off " 
	 << width/mm << " mm of "
	 << matname << G4endl;
  G4cout << "###    Beam momentum spread " << sigma << " MeV/c" << G4endl;
  fout << matname << "  " << ptot << G4endl;
  for(i=0; i<11; i++) {
    fout << tbins[i] << "  " << x[i] << "  " << dx[i] << G4endl;
  }
  fout.close();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::AddProjectileTheta(G4double theta)
{
  for(G4int i=0; i<11; i++) {
    if(theta <= tbins[i]) {
      hist[i] += 1;
      nevt++;
      break;
    }
  }
  return;
}
