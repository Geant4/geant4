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
#include "globals.hh"
#include "G4ios.hh"

#include "HistoEPTest47.hh"
#include "G4VParticleChange.hh"
#include "G4UnitsTable.hh"

HistoEPTest47::HistoEPTest47() : hiPx(0), hiPy(0), hiPz(0), hiE(0), hiKE(0) {

  G4cout << "HistoEPTest47:: Initialized" << G4endl;

}

HistoEPTest47::~HistoEPTest47() {
  if (hiPx > 0) delete hiPx;
  if (hiPy > 0) delete hiPy;
  if (hiPz > 0) delete hiPz;
  if (hiE  > 0) delete hiE;
  if (hiKE > 0) delete hiKE;
}

void HistoEPTest47::fill(G4VParticleChange* aChange, G4LorentzVector pinit) {

  G4LorentzVector labv(pinit), fm;
  G4int n = aChange->GetNumberOfSecondaries();
  const G4DynamicParticle* sec = 0;
  G4ThreeVector  mom;
  G4ParticleDefinition* pd;
  G4double  keBalance=labv.e()-labv.m();
  for (G4int j=0; j<n; j++) {
    sec = aChange->GetSecondary(j)->GetDynamicParticle();
    pd  = sec->GetDefinition();
    mom = sec->GetMomentumDirection();
    G4double ke    = (sec->GetKineticEnergy())/GeV;
    if (ke < 0.0) ke = 0.0;
    G4double m     = (pd->GetPDGMass())/GeV;
    G4double p     = std::sqrt(ke*(ke + 2.0*m));
    G4double ee    = ke + m;
    mom           *= p;
    fm             = G4LorentzVector(mom, ee);
    labv          -= fm;
    keBalance     -= ke;
  }

  //  G4cout << "Balance " << labv << " KE " << keBalance << G4endl;
  hiPx->Fill(labv.px());
  hiPy->Fill(labv.py());
  hiPz->Fill(labv.pz());
  hiE->Fill(labv.e());
  hiKE->Fill(keBalance);
}

void HistoEPTest47::write() {

  hiPx->Write(); hiPy->Write(); hiPz->Write();
  hiE->Write(); hiKE->Write(); 
}

void HistoEPTest47::initialize(std::string namePart, std::string nameMat, 
			       G4double momentum, std::string nameGen) {


  char name[100], title[160];
  if (hiPx > 0) delete hiPx;
  if (hiPy > 0) delete hiPy;
  if (hiPz > 0) delete hiPz;
  if (hiE  > 0) delete hiE;
  if (hiKE > 0) delete hiKE;
  sprintf (name, "PxBalance");
  sprintf (title, "Px Balance in %s + %s at %6.2f GeV/c (%s)",namePart.c_str(),
	   nameMat.c_str(), momentum, nameGen.c_str());
  hiPx = new TH1F (name, title, 5000, -5.0, 5.0);
  sprintf (name, "PyBalance");
  sprintf (title, "Py Balance in %s + %s at %6.2f GeV/c (%s)",namePart.c_str(),
	   nameMat.c_str(), momentum, nameGen.c_str());
  hiPy = new TH1F (name, title, 5000, -5.0, 5.0);
  sprintf (name, "PzBalance");
  sprintf (title, "Pz Balance in %s + %s at %6.2f GeV/c (%s)",namePart.c_str(),
	   nameMat.c_str(), momentum, nameGen.c_str());
  hiPz = new TH1F (name, title, 5000, -5.0, 5.0);
  sprintf (name, "EBalance");
  sprintf (title, "E Balance in %s + %s at %6.2f GeV/c (%s)", namePart.c_str(),
	   nameMat.c_str(), momentum, nameGen.c_str());
  hiE  = new TH1F (name, title, 5000, -5.0, 195.0);
  sprintf (name, "KEBalance");
  sprintf (title, "KE Balance in %s + %s at %6.2f GeV/c (%s)",namePart.c_str(),
	   nameMat.c_str(), momentum, nameGen.c_str());
  hiKE = new TH1F (name, title, 5000, -25.0, 25.0);
}
