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
//      File name:     G4ExcitationHandlerTest.cc 
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
#include <fstream>
#include <iomanip>
#include <iostream>

#include "CLHEP/Hist/TupleManager.h"
#include "CLHEP/Hist/HBookFile.h"
#include "CLHEP/Hist/Histogram.h"
#include "CLHEP/Hist/Tuple.h"
#include "CLHEP/String/Strings.h"
#include "G4VEvaporationChannel.hh"
#include "G4CompetitiveFission.hh"

#include "G4ReactionProductVector.hh"
#include "G4ReactionProduct.hh"
#include "G4DynamicParticleVector.hh"
#include "G4DynamicParticle.hh"
#include "G4ParticleDefinition.hh"

#include "G4ExcitationHandler.hh"
#include "G4EvaporationChannel.hh"
#include "G4PhotonEvaporation.hh"
#include "G4DataVector.hh"
#include "G4LorentzVector.hh"
#include "G4ThreeVector.hh"
#include "G4NucleiProperties.hh"
#include "G4Fragment.hh"
#include "Randomize.hh"
#include "G4Material.hh"
#include "G4Element.hh"
#include "G4FragmentVector.hh"
#include "G4ParticleDefinition.hh"
#include "G4IonTable.hh"
#include "G4Electron.hh"
#include "G4Proton.hh"
#include "G4Neutron.hh"
#include "G4ParticleTypes.hh"
#include "G4IonTable.hh"
#include "G4GenericIon.hh"
#include "G4NucleiPropertiesTable.hh"

#include "g4rw/tvvector.h"

int main()
{
  // MGP ---- HBOOK initialization
  HepTupleManager* hbookManager;
  hbookManager = new HBookFile("handler.hbook", 58);
  assert (hbookManager != 0);

  // MGP ---- Book histograms

  HepHistogram* hGammaE;
  hGammaE = hbookManager->histogram("Gamma energy", 200,0.,20.);
  assert (hGammaE != 0);  

  HepHistogram* hNProducts;
  hNProducts = hbookManager->histogram("Number of products", 100,0.,100.);
  assert (hNProducts != 0);  

  HepHistogram* hProb;
  hProb = hbookManager->histogram("Probability * 1.e25", 100,0.,10.);
  assert (hProb != 0);  

  // MGP ---- Book a ntuple
  HepTuple* ntuple;
  ntuple = hbookManager->ntuple("G4ExcitationHandler");
  assert (ntuple != 0);

  G4int Z;
  G4int A;

  G4cout << "Enter Z and A" << G4endl;
  G4cin >> Z >> A;

  assert (Z > 0);
  assert (A > 0);
  assert (A > Z);
  
  G4int iter;
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

  G4ParticleDefinition *theElectron = G4Electron::ElectronDefinition();
  G4ParticleDefinition* ion  = G4GenericIon::GenericIonDefinition();
  G4ParticleTable* theTableOfParticles;
  theTableOfParticles = G4ParticleTable::GetParticleTable();
  G4IonTable* theTable = new G4IonTable();

  // G4ExcitationHandler for this (Z,A) material
  // Excitation energy levels for each channel

  G4int i;
  for (i=0; i<iter; i++)
    {
      G4double excitation = excMin + G4UniformRand() * (excMax - excMin);

      G4cout << G4endl << "TEST >>>>>>>>> Iteration " << i 
	     << " <<<<<<<<< Initial excitation " << excitation << G4endl;

      G4double px = 0.;
      G4double py = 0.;
      G4double pz = 0.;

      G4ThreeVector p(px,py,pz);
      G4double nuclearMass = G4NucleiProperties::GetNuclearMass(A,Z);
      G4LorentzVector p4(p,std::sqrt((nuclearMass+excitation)*(nuclearMass+excitation)));

      G4Fragment nucleus(A,Z,p4);

      G4ExcitationHandler* handler = new G4ExcitationHandler;
      //      G4DynamicParticleVector* products = handler->BreakItUp(nucleus);
      G4ReactionProductVector* products = handler->BreakItUp(nucleus);

      // Fill histograms 

      G4int nProducts = 0;
      if (products !=0) nProducts = products->entries();

      G4cout << nProducts << " products: gammaE ";
      
      G4int ig = 0;
      for (ig=0; ig<nProducts; ig++)
	{
	  if (products->at(ig)->GetDefinition() == G4Gamma::GammaDefinition())
	    {
	      //	      G4double productE = products->at(ig)->Get4Momentum().e();
	      G4double productE = products->at(ig)->GetTotalEnergy();
	      G4cout << productE << " ";
	      hGammaE->accumulate(productE);
              ntuple->column("egamma",productE);
	      ntuple->dumpData();
	    }
	}
      G4cout << G4endl;
      
      delete handler;
      products->clearAndDestroy();
      delete products;
    }

  hbookManager->write();

  
  return EXIT_SUCCESS;
}



















