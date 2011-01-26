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
// -------------------------------------------------------------------
//      GEANT 4 class file --- Copyright CERN 1998
//      CERN Geneva Switzerland
//
//
//      File name:     G4PhotonEvaporationTest.cc 
//
//      Author:        Maria Grazia Pia (pia@genova.infn.it), 
// 
//      Creation date: 27 October 1998
//
//      Modifications: 
//
// -------------------------------------------------------------------

#include "globals.hh"

#include "G4ios.hh"
#include <assert.h>
#include <fstream>
#include <iomanip>
#include <iostream>

#include "CLHEP/Hist/TupleManager.h"
#include "CLHEP/Hist/HBookFile.h"
#include "CLHEP/Hist/Histogram.h"
#include "CLHEP/Hist/Tuple.h"

#include "G4PhotonEvaporation.hh"
#include "G4VGammaDeexcitation.hh"
#include "G4DataVector.hh"
#include "G4ContinuumGammaDeexcitation.hh"
#include "G4DiscreteGammaDeexcitation.hh"
#include "G4LorentzVector.hh"
#include "G4ThreeVector.hh"
#include "G4NucleiProperties.hh"
#include "G4NucleiPropertiesTable.hh"
#include "G4Fragment.hh"
#include "Randomize.hh"
#include "G4FragmentVector.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleDefinition.hh"

int main()
{
  // MGP ---- HBOOK initialization
  HepTupleManager* hbookManager;
  hbookManager = new HBookFile("mgcarlo.hbook", 58);
  assert (hbookManager != 0);

  // MGP ---- Book histograms

  HepHistogram* hGammaE;
  hGammaE = hbookManager->histogram("Gamma energy", 100,0.,10.);
  assert (hGammaE != 0);  

  HepHistogram* hNProducts;
  hNProducts = hbookManager->histogram("Number of products", 20,0.,20.);
  assert (hNProducts != 0);  

  HepHistogram* hProb;
  hProb = hbookManager->histogram("Probability * 1.e25", 100,0.,10.);
  assert (hProb != 0);  

  // MGP ---- Book a ntuple
  HepTuple* ntuple;
  ntuple = hbookManager->ntuple("G4PhotonEvaporation");
  assert (ntuple != 0);

  HepTuple* ntupleCons;
  ntupleCons = hbookManager->ntuple("Conservation");
  assert (ntupleCons != 0);

  G4int Z;
  G4int A;

  G4cout << "Enter Z and A" << G4endl;
  G4cin >> Z >> A;

  assert (Z > 0);
  assert (A > 0);
  assert (A > Z);
  
 
  G4int verbose;
  G4cout << "Enter verbose level " << G4endl;
  G4cin >> verbose;

  G4int mode = 0;
  G4cout << "Enter continuum (0) or discrete (1) deexcitation mode" << G4endl;
  G4cin >> mode;

  G4int iter = 1;
  G4cout << "Enter number of iterations " << G4endl;
  G4cin >> iter;
  if (iter <1) iter = 1;

  G4double excMin;
  G4double excMax;    
  G4cout << "Enter initial min an max excitation energy" << G4endl;
  G4cin >> excMin >> excMax;
  assert (excMin >= 0.);
  assert (excMax > 0.);
  assert (excMax >= excMin);
 

  G4ParticleDefinition* gamma = G4Gamma::GammaDefinition();

  G4ParticleDefinition* proton = G4Proton::ProtonDefinition();
  G4ParticleDefinition* antiProton = G4AntiProton::AntiProtonDefinition();
  G4ParticleDefinition* neutron = G4Neutron::NeutronDefinition();   
  G4ParticleDefinition* antiNeutron = G4AntiNeutron::AntiNeutronDefinition();   

  G4ParticleDefinition* pionPlus = G4PionPlus::PionPlusDefinition();
  G4ParticleDefinition* pionMinus = G4PionMinus::PionMinusDefinition();
  G4ParticleDefinition* pionZero = G4PionZero::PionZeroDefinition();

  G4ParticleDefinition* kaonPlus = G4KaonPlus::KaonPlusDefinition();
  G4ParticleDefinition* kaonMinus = G4KaonMinus::KaonMinusDefinition();
   
  G4ParticleDefinition* lambda = G4Lambda::LambdaDefinition();

  G4ParticleDefinition* theMuonPlus = G4MuonPlus::MuonPlusDefinition();
  G4ParticleDefinition* theMuonMinus = G4MuonMinus::MuonMinusDefinition();

  G4ParticleDefinition* theNeutrinoMu = G4NeutrinoMu::NeutrinoMuDefinition();
  G4ParticleDefinition* theAntiNeutrinoMu = G4AntiNeutrinoMu::AntiNeutrinoMuDefinition();

  // G4PhotonEvaporation for this (Z,A) material
  G4PhotonEvaporation* photonEvaporation = new G4PhotonEvaporation;
  G4cout << "TEST ---- G4PhotonEvaporation created ----" << G4endl;

  photonEvaporation->SetVerboseLevel(verbose);

  G4int i;
  for (i=0; i<iter; i++)
    {
      G4double excitation = excMin + G4UniformRand() * (excMax - excMin);

	  G4cout << G4endl << "TEST >>>>>>>>> Iteration " << i 
		 << " <<<<<<<<< Initial excitation " << excitation << " ";

      if (mode > 0) 
	{
	  // Transform excitation energy into nearest level energy 

	  G4NuclearLevelManager levelManager;
	  levelManager.SetNucleus(Z,A);
          const G4NuclearLevel* level = levelManager.NearestLevel(excitation);
	  if (level != 0)
	    {
	      excitation = level->Energy();
	      G4cout << "Transformed into excitation " << excitation << G4endl;
	    }
	}
      else
	{
	  G4cout << G4endl;
	}

      G4LorentzVector p4(0.,0.,0.,G4NucleiProperties::GetNuclearMass(A,Z)+excitation);
      G4Fragment nucleus(A,Z,p4);
      //      nucleus.SetExcitationEnergy(excitation);
      G4double initialE = nucleus.GetMomentum().e();
         
      G4FragmentVector* products = photonEvaporation->BreakItUp(nucleus);
      if(verbose >= 0) G4cout << "TEST: BreakItUp done" << G4endl;

      G4double scale = 1./millibarn;
      G4double prob = photonEvaporation->GetEmissionProbability() * scale;
      if(verbose >= 0) G4cout << "TEST: probability " << prob << G4endl;
 
      // Fill histograms 

      G4int nProducts = 0;
      if (products !=0) nProducts = products->entries();
      hNProducts->accumulate(nProducts);

      hProb->accumulate(prob);

      if(verbose >= 0) G4cout << "TEST: " << nProducts << " products generated: gammaE ";

      G4double currentExcitation = excitation;
      G4double sumPx = 0.;
      G4double sumPy = 0.;
      G4double sumPz = 0.;
      G4double sumE = 0.;
      G4double sumEgamma = 0.;

      G4int ig = 0;
      for (ig=0; ig<nProducts; ig++)
	{
	  G4double productE = products->at(ig)->GetMomentum().e();
	  G4ThreeVector pProd(products->at(ig)->GetMomentum());
	  sumPx = sumPx + pProd.x();
	  sumPy = sumPy + pProd.y();
	  sumPz = sumPz + pProd.z();
	  sumE = sumE + productE;
	  if (products->at(ig)->GetA() < 1)
	    {
	      if(verbose >= 0) G4cout << productE << " ";
	      sumEgamma = sumEgamma + productE;
	      hGammaE->accumulate(productE);
	      ntuple->column("exc",currentExcitation);
              ntuple->column("egamma",productE);
              ntuple->column("prob",prob);
	      ntuple->dumpData();
	      // ntuple->column("px",pProd.x());
	      // ntuple->column("py",pProd.y());
	      // ntuple->column("pz",pProd.z());
	      currentExcitation = currentExcitation - productE;
	    }
	  else
	    { if(verbose >= 0) G4cout << "(" << products->at(ig)->GetZ() << "," << products->at(ig)->GetA() << ") "; }
	}
      if(verbose >= 0) G4cout << G4endl;

      if(verbose >= 0) G4cout << "TEST: Sum of product energies = " << sumE 
	     << ", Sum of gamma energies = " << sumEgamma << G4endl;

      if (products != 0) 
	{
	  products->clearAndDestroy();

	  ntupleCons->column("px",sumPx);
	  ntupleCons->column("py",sumPy);
	  ntupleCons->column("pz",sumPz);
	  ntupleCons->column("eprod",sumE);
	  ntupleCons->column("exc",excitation);
	  ntupleCons->column("egammas",sumEgamma);
	  ntupleCons->column("ein",initialE);
	  ntupleCons->dumpData();
	}
      delete products;
    }

  hbookManager->write();

  delete photonEvaporation;
  
  return EXIT_SUCCESS;
}



















