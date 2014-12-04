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
#ifdef TEST69_HAS_ROOT

#include "Tst69INCLXXTallyROOT.hh"
#include "G4SystemOfUnits.hh"
#include "G4HadSecondary.hh"
#include "G4DynamicParticle.hh"
#include "CLHEP/Units/PhysicalConstants.h"
#include <algorithm>

Tst69INCLXXTallyROOT::Tst69INCLXXTallyROOT(TString const &physList) :
  file(NULL),
  reac(NULL),
  Ap(0),
  Zp(0),
  Ep(0.),
  At(0),
  Zt(0),
  nParticles(0)
{
  filename = "test69_rtally";
  filename += physList;
  filename += ".root";
  std::fill_n(A, maxNParticles, 0);
  std::fill_n(Z, maxNParticles, 0);
  std::fill_n(EKin, maxNParticles, 0.);
  std::fill_n(px, maxNParticles, 0.);
  std::fill_n(py, maxNParticles, 0.);
  std::fill_n(pz, maxNParticles, 0.);
  std::fill_n(theta, maxNParticles, 0.);
  std::fill_n(phi, maxNParticles, 0.);
}

Tst69INCLXXTallyROOT::~Tst69INCLXXTallyROOT() {
  Close();
}

void Tst69INCLXXTallyROOT::Open() {
  // Open an output file
  file = new TFile(filename, "RECREATE");
  if (!file->IsOpen()) {
    G4cout << "\n---> Tst69INCLXXTallyROOT::Open(): cannot open " << filename << G4endl;
    delete file;
    file = NULL;
    return;
  }

  // Create tree
  reac = new TTree("reac", "INCL++ event tally");

  // Branches for theEventTree
  reac->Branch("Ap", &(Ap), "Ap/S");
  reac->Branch("Zp", &(Zp), "Zp/S");
  reac->Branch("Ep", &(Ep), "Ep/F");
  reac->Branch("At", &(At), "At/S");
  reac->Branch("Zt", &(Zt), "Zt/S");
  reac->Branch("nParticles", &(nParticles), "nParticles/S");
  reac->Branch("A", A, "A[nParticles]/S");
  reac->Branch("Z", Z, "Z[nParticles]/S");
  reac->Branch("EKin", EKin, "EKin[nParticles]/F");
  reac->Branch("px", px, "px[nParticles]/F");
  reac->Branch("py", py, "py[nParticles]/F");
  reac->Branch("pz", pz, "pz[nParticles]/F");
  reac->Branch("theta", theta, "theta[nParticles]/F");
  reac->Branch("phi", phi, "phi[nParticles]/F");

  G4cout << "\n----> Tree is opened in " << filename << G4endl;
}

void Tst69INCLXXTallyROOT::Close() {
  if(file) {
    file->Write();
    file->Close();  
    file = NULL;
    reac = NULL;
    G4cout << "\n----> Tree is saved in " << filename << G4endl;
  }                    
}

void Tst69INCLXXTallyROOT::Tally(G4HadProjectile const &aTrack, G4Nucleus const &theNucleus, G4HadFinalState const &result) {
  if(!reac)
    return;

  G4ParticleDefinition const * const trackDefinition = aTrack.GetDefinition();

  // Fill ntuple
  Ap = trackDefinition->GetAtomicMass();
  Zp = (Short_t) trackDefinition->GetPDGCharge();
  Ep = aTrack.GetKineticEnergy()/MeV;
  At = theNucleus.GetA_asInt();
  Zt = theNucleus.GetZ_asInt();
  nParticles = (Short_t) result.GetNumberOfSecondaries();
  for(Short_t i=0; i<nParticles; ++i) {
    G4DynamicParticle const * const p = result.GetSecondary(i)->GetParticle();
    G4ParticleDefinition const * const pdef = p->GetDefinition();
    A[i] = pdef->GetAtomicMass();
    Z[i] = pdef->GetPDGCharge();
    EKin[i] = p->GetKineticEnergy()/MeV;
    G4ThreeVector const pmom = p->GetMomentum();
    px[i] = pmom.x();
    py[i] = pmom.y();
    pz[i] = pmom.z();
    theta[i] = 180.*pmom.theta()/CLHEP::pi;
    phi[i] = 180.*pmom.phi()/CLHEP::pi;
  }
  reac->Fill();
  
}

#endif // TEST69_HAS_ROOT
