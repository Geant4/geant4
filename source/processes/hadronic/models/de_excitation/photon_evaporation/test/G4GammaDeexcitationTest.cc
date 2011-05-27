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
//      File name:     G4GammaDeexcitationTest.cc 
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
#include <assert.h>

#include "CLHEP/Hist/TupleManager.h"
#include "CLHEP/Hist/HBookFile.h"
#include "CLHEP/Hist/Histogram.h"
#include "CLHEP/Hist/Tuple.h"

#include "G4VGammaDeexcitation.hh"
#include "G4DataVector.hh"
#include "G4ContinuumGammaDeexcitation.hh"
#include "G4DiscreteGammaDeexcitation.hh"
#include "G4LorentzVector.hh"
#include "G4NucleiProperties.hh"
#include "G4NucleiPropertiesTable.hh"
#include "G4Fragment.hh"
#include "Randomize.hh"
#include "G4FragmentVector.hh"

int main()
{
  // MGP ---- HBOOK initialization
  HepTupleManager* hbookManager;
  hbookManager = new HBookFile("gammadeex.hbook", 58);
  assert (hbookManager != 0);

  // MGP ---- Book histograms

  HepHistogram* hGammaE;
  hGammaE = hbookManager->histogram("Gamma energy", 100,0.,10.);
  assert (hGammaE != 0);  

  HepHistogram* hNGammas;
  hNGammas = hbookManager->histogram("Number of gammas", 20,0.,20.);
  assert (hNGammas != 0);  

  // MGP ---- Book a ntuple
  HepTuple* ntuple;
  ntuple = hbookManager->ntuple("G4GammaDeexcitation ntuple");
  assert (ntuple != 0);

  G4int Z;
  G4int A;

  G4cout << "Enter Z and A" << G4endl;
  G4cin >> Z >> A;

  assert (Z > 0);
  assert (A > 0);
  assert (A > Z);
  
  G4NuclearLevelManager levelManager(Z,A);

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

  G4int i;
  for (i=0; i<iter; i++)
    {
      G4double excitation = excMin + G4UniformRand() * (excMax - excMin);
      G4VGammaDeexcitation* deexcitation;

      G4cout << G4endl << "TEST >>>>>>>>> Iteration " << i 
	     << " <<<<<<<<< Initial excitation " << excitation << " ";

      // Discrete deexcitation

      if (mode > 0) 
	{
	  // Transform excitation energy into nearest level energy 
	  const G4NuclearLevel* level = levelManager.NearestLevel(excitation);
	  if (level != 0)
	    {
	      excitation = level->Energy();
	      G4cout << "Transformed into excitation " << excitation << G4endl;
	    }

	  deexcitation = new G4DiscreteGammaDeexcitation();
	  G4cout << "TEST ---- G4DiscreteGammaDeexcitation created ----" << G4endl;
	}

      // Continuum deexcitation
      
      else
	{
	  deexcitation = new G4ContinuumGammaDeexcitation();
	  G4cout << G4endl << "TEST ---- G4ContinuumGammaDeexcitation created ----" << G4endl;
	}
      
      deexcitation->SetVerboseLevel(verbose);

      G4LorentzVector p4(0.,0.,0.,G4NucleiProperties::GetNuclearMass(A,Z));
      G4Fragment nucleus(A,Z,p4);
      nucleus.SetExcitationEnergy(excitation);
      deexcitation->SetNucleus(nucleus);
      
      G4FragmentVector* gammas = deexcitation->DoChain();

      // Fill histograms 
      G4int nGammas = 0;
      if (gammas !=0) nGammas = gammas->entries();
      hNGammas->accumulate(nGammas);
      G4cout << "TEST: " << nGammas << " photons generated: ";
      G4int ig = 0;
      for (ig=0; ig<nGammas; ig++)
	{
	  G4double gammaE = gammas->at(ig)->GetMomentum().e();
	  G4cout << gammaE << " ";
          hGammaE->accumulate(gammaE);
	}
      G4cout << G4endl;

      if (gammas != 0) gammas->clearAndDestroy();
      delete gammas;
      gammas = 0;
      delete deexcitation;
      deexcitation = 0;
    }

  hbookManager->write();
  
  return EXIT_SUCCESS;
}





