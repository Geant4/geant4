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
// $Id: Tst26RunAction.cc,v 1.3 2003-02-06 11:53:27 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// 
/////////////////////////////////////////////////////////////////////////
//
// test26: Cut per region physics
//
// Created: 31.01.03 V.Ivanchenko
//
// Modified:
//
////////////////////////////////////////////////////////////////////////
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "Tst26RunAction.hh"

#include "Tst26PrimaryGeneratorAction.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4ios.hh"
#include "G4UnitsTable.hh"
#include "Randomize.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Tst26RunAction::Tst26RunAction(Tst26PrimaryGeneratorAction* kin)
  :Tst26Kin(kin),
   trancFactor(0.8)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Tst26RunAction::~Tst26RunAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Tst26RunAction::BeginOfRunAction(const G4Run* aRun)
{
  beamE = Tst26Kin->GetParticleGun()->GetParticleEnergy();   
  G4cout << "### Run " << aRun->GetRunID() << " start  beam E= " 
         << beamE/MeV << " MeV"
         << G4endl;
  evtNo = 0;
  calNo = 0; 
  for (G4int i=0; i<8; i++) {
    e[i] = 0.0;
    e2[i] = 0.0;
  }
  sgam  = 0;
  sgamV = 0;
  sgamM = 0;
  sel   = 0;
  selV  = 0;
  selM  = 0;
  spos  = 0;
  sposV = 0;
  sposM = 0;

  // save Rndm status
  G4RunManager::GetRunManager()->SetRandomNumberStore(true);
  HepRandom::showEngineStatus();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Tst26RunAction::EndOfRunAction(const G4Run* aRun)
{
  const G4String pn = Tst26Kin->GetParticleGun()
                      ->GetParticleDefinition()->GetParticleName();   
  G4cout.setf(G4std::ios::fixed,G4std::ios::floatfield);
  G4cout.precision(5);

  G4double x, y, z;     
  G4double ntot = (G4double)evtNo;
  if(ntot == 0.0) ntot = 1.0;

  G4double norm = 1./ntot;
  for(G4int i=0; i<8; i++) {
    e[i]  *= norm;
    e2[i] *= norm;
    e2[i] -= e[i]*e[i];
    e2[i]  = sqrt(e2[i]);
    s[i]   = e2[i]*sqrt(norm);
  }

  G4cout << G4endl;
  G4cout << "============   SUMMARY   ============" << G4endl;
  G4cout << G4endl;
  G4cout << "   " << evtNo << " Events were generated" << G4endl; 
  G4cout << "   Incident particle " << pn << " Ekin(meV) = " << beamE/MeV << G4endl; 
  G4cout << G4endl;
  G4cout << "    Average number of secondaries" << G4endl;
  G4cout << G4endl;
  G4cout << "             World         Vertex Detector       Muon Detector" << G4endl;
  G4cout << G4endl;
  x = ((G4double)sgam)*norm;
  y = ((G4double)sgamV)*norm;
  z = ((G4double)sgamM)*norm;
  G4cout << "  Gammas:    " << x << "        " << y << "        " << z << G4endl;
  x = ((G4double)sel)*norm;
  y = ((G4double)selV)*norm;
  z = ((G4double)selM)*norm;
  G4cout << "  Electrons: " << x << "        " << y << "        " << z << G4endl;
  x = ((G4double)spos)*norm;
  y = ((G4double)sposV)*norm;
  z = ((G4double)selM)*norm;
  G4cout << "  Positrons:  " << x << "        " << y << "        " << z << G4endl;
  G4cout << G4endl;
  G4cout << "    Normalised energy deposition in the calorimeter" 
         << G4endl;
  G4cout << G4endl;
  G4cout << "  E1 = " << e[0] << " +- " << s[0] 
	 << " RMS(E1) = " << e2[0] <<  " +- " << s[0] << G4endl;
  G4cout << "  E9 = " << e[1] << " +- " << s[1] 
	 << " RMS(E9) = " << e2[1] <<  " +- " << s[1] << G4endl;
  G4cout << "  E25= " << e[2] << " +- " << s[2] 
	 << " RMS(E25)= " << e2[2] <<  " +- " << s[2] << G4endl;
  G4cout << G4endl;
  G4cout << "    Energy deposition in the absorbers" 
         << G4endl;
  G4cout << G4endl;
  G4cout << "  Eabs1(MeV)= " << e[3]/MeV << " +- " << s[3]/MeV 
	 << " RMS(Eabs1)= " << e2[3]/MeV <<  " +- " << s[3]/MeV << G4endl;
  G4cout << "  Eabs2(MeV)= " << e[4]/MeV << " +- " << s[4]/MeV 
	 << " RMS(Eabs2)= " << e2[4]/MeV <<  " +- " << s[4]/MeV << G4endl;
  G4cout << "  Eabs3(MeV)= " << e[5]/MeV << " +- " << s[5]/MeV 
	 << " RMS(Eabs3)= " << e2[5]/MeV <<  " +- " << s[5]/MeV << G4endl;
  G4cout << "  Eabs4(MeV)= " << e[6]/MeV << " +- " << s[6]/MeV 
	 << " RMS(Eabs4)= " << e2[6]/MeV <<  " +- " << s[6]/MeV << G4endl;

  G4cout << G4endl;
  G4cout << "    Average number of hits in the vertex detector" 
         << G4endl;
  G4cout << G4endl;
  G4cout << "  Nhit= " << e[7] << " +- " << s[7] 
	 << G4endl;
                   

  // show Rndm status
  HepRandom::showEngineStatus();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Tst26RunAction::AddParticle(G4int partIndex, G4int regionIndex)
{
  if(partIndex == 0) {
    if(regionIndex == 0)      sgam++;
    else if(regionIndex == 1) sgamV++;
    else if(regionIndex == 2) sgamM++;
  } else if(partIndex == 1) {
    if(regionIndex == 0)      sel++;
    else if(regionIndex == 1) selV++;
    else if(regionIndex == 2) selM++;
  } else if(partIndex == 2) {
    if(regionIndex == 0)      spos++;
    else if(regionIndex == 1) sposV++;
    else if(regionIndex == 2) sposM++;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Tst26RunAction::AddEvent(G4double E1, G4double E9, G4double E25,
                              G4double Eabs1, G4double Eabs2, 
                              G4double Eabs3, G4double Eabs4, 
                              G4int nPad)
{
  evtNo++;
  E1  /= beamE;
  E9  /= beamE;
  E25 /= beamE;
  e[0]  += E1;
  e2[0] += E1*E1;
  e[1]  += E9;
  e2[1] += E9*E9;
  e[2]  += E25;
  e2[2] += E25*E25;
 
  e[3]  += Eabs1;
  e2[3] += Eabs1*Eabs1;
  e[4]  += Eabs2;
  e2[4] += Eabs2*Eabs2;
  e[5]  += Eabs3;
  e2[5] += Eabs3*Eabs3;
  e[6]  += Eabs4;
  e2[6] += Eabs4*Eabs4;

  e[7]  += nPad;
  e2[7] += nPad*nPad;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

