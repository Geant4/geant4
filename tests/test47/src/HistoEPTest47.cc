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
#include "CLHEP/Units/PhysicalConstants.h"

HistoEPTest47::HistoEPTest47() : pin(0), debug(false), hiPx(0), hiPy(0), 
				 hiPz(0), hiE(0), hiKE(0), hiX(0), hiY(0),
				 hiDpByp(0), hiR(0), hiZ(0), hiAngX(0), 
				 hiAngY(0) {

  G4cout << "HistoEPTest47:: Initialized" << G4endl;

}

HistoEPTest47::~HistoEPTest47() {

  if (hiPx)    delete hiPx;
  if (hiPy)    delete hiPy;
  if (hiPz)    delete hiPz;
  if (hiE)     delete hiE;
  if (hiKE)    delete hiKE;
  if (hiDpByp) delete hiDpByp;
  if (hiR)     delete hiR;
  if (hiX)     delete hiX;
  if (hiY)     delete hiY;
  if (hiZ)     delete hiZ;
  if (hiAngX)  delete hiAngX;
  if (hiAngY)  delete hiAngY;

}

void HistoEPTest47::fill(G4VParticleChange* aChange, G4LorentzVector& pinit,
			 G4ThreeVector& aPosition) {

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
    G4double ke    = (sec->GetKineticEnergy())/CLHEP::GeV;
    if (ke < 0.0) ke = 0.0;
    G4double m     = (pd->GetPDGMass())/CLHEP::GeV;
    G4double p     = std::sqrt(ke*(ke + 2.0*m));
    G4double ee    = ke + m;
    mom           *= p;
    fm             = G4LorentzVector(mom, ee);
    labv          -= fm;
    keBalance     -= ke;
  }

  if (debug) G4cout << "Balance " << labv << " KE " << keBalance << G4endl;
  hiPx->Fill(labv.px());
  hiPy->Fill(labv.py());
  hiPz->Fill(labv.pz());
  hiE->Fill(labv.e());
  hiKE->Fill(keBalance);

  G4double dpbyp = 0;
  if (pin > 0) dpbyp = pinit.rho()/pin;
  G4double x     = (aPosition.x()/CLHEP::cm);
  G4double y     = (aPosition.y()/CLHEP::cm);
  G4double z     = (aPosition.z()/CLHEP::cm);
  G4double r     = std::sqrt(x*x + y*y);
  G4double angx  = std::asin(pinit.y()/pinit.rho());
  G4double angy  =-std::atan2(pinit.x(),pinit.z());
  if (debug) G4cout << "HistoEPTest47:: DpByP " << dpbyp << " R:Z " << r << ":"
		    << z << " Angle(X:Y) " << angx << ":" << angy << G4endl;
  hiDpByp->Fill(dpbyp);
  hiR    ->Fill(r);
  hiX    ->Fill(x);
  hiY    ->Fill(y);
  hiZ    ->Fill(z);
  hiAngX ->Fill(angx);
  hiAngY ->Fill(angy);
}

void HistoEPTest47::write() {

  hiPx->Write();    hiPy->Write(); hiPz->Write();
  hiE->Write();     hiKE->Write(); hiX->Write(); hiY->Write();
  hiDpByp->Write(); hiR->Write();  hiZ->Write();
  hiAngX->Write();  hiAngY->Write();
}

void HistoEPTest47::initialize(std::string namePart, std::string nameMat, 
			       G4double momentum, std::string nameGen) {

  pin = momentum;
  char name[100], title[160];
  if (hiPx)    delete hiPx;
  if (hiPy)    delete hiPy;
  if (hiPz)    delete hiPz;
  if (hiE)     delete hiE;
  if (hiKE)    delete hiKE;
  if (hiDpByp) delete hiDpByp;
  if (hiR)     delete hiR;
  if (hiX)     delete hiX;
  if (hiY)     delete hiY;
  if (hiZ)     delete hiZ;
  if (hiAngX)  delete hiAngX;
  if (hiAngY)  delete hiAngY;
  sprintf (name, "PxBalance");
  sprintf (title, "Px Balance in %s + %s at %6.2f GeV/c (%s)",namePart.c_str(),
	   nameMat.c_str(), momentum, nameGen.c_str());
  hiPx = new TH1F (name, title, 5000, -5.0, 5.0);
  hiPx->Sumw2();
  sprintf (name, "PyBalance");
  sprintf (title, "Py Balance in %s + %s at %6.2f GeV/c (%s)",namePart.c_str(),
	   nameMat.c_str(), momentum, nameGen.c_str());
  hiPy = new TH1F (name, title, 5000, -5.0, 5.0);
  hiPy->Sumw2();
  sprintf (name, "PzBalance");
  sprintf (title, "Pz Balance in %s + %s at %6.2f GeV/c (%s)",namePart.c_str(),
	   nameMat.c_str(), momentum, nameGen.c_str());
  hiPz = new TH1F (name, title, 5000, -5.0, 5.0);
  hiPz->Sumw2();
  sprintf (name, "EBalance");
  sprintf (title, "E Balance in %s + %s at %6.2f GeV/c (%s)", namePart.c_str(),
	   nameMat.c_str(), momentum, nameGen.c_str());
  hiE  = new TH1F (name, title, 5000, -5.0, 195.0);
  hiE->Sumw2();
  sprintf (name, "KEBalance");
  sprintf (title, "KE Balance in %s + %s at %6.2f GeV/c (%s)",namePart.c_str(),
	   nameMat.c_str(), momentum, nameGen.c_str());
  hiKE = new TH1F (name, title, 5000, -25.0, 25.0);
  hiKE->Sumw2();

  hiDpByp = new TH1F ("DpByp", "#delta p/p of beam", 1600, 0.8, 1.2);
  hiDpByp->Sumw2();
  hiR     = new TH1F ("RPos",  "R of interaction vertex", 100, 0, 10.0);
  hiR->Sumw2();
  hiX     = new TH1F ("XPos",  "X of interaction vertex", 400, -10.0, 10.0);
  hiX->Sumw2();
  hiY     = new TH1F ("YPos",  "Y of interaction vertex", 400, -10.0, 10.0);
  hiY->Sumw2();
  hiZ     = new TH1F ("ZPos",  "Z of interaction vertex", 400, -10.0, 10.0);
  hiZ->Sumw2();
  hiAngX  = new TH1F ("AngleX","Angle of the beam in YZ-plane", 400, -0.1,0.1);
  hiAngX->Sumw2();
  hiAngY  = new TH1F ("AngleY","Angle of the beam in XZ plane", 400, -0.1,0.1);
  hiAngY->Sumw2();

}
