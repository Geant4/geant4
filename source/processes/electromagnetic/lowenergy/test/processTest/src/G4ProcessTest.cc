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
// $Id: G4ProcessTest.cc,v 1.2 2001-10-28 18:00:34 pia Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Author: Maria Grazia Pia (Maria.Grazia.Pia@cern.ch)
//
// History:
// -----------
// 07 Oct 2001   MGP        Created
//
// -------------------------------------------------------------------
// Class description:
// Test DoIt method of physics processes
// Further documentation available from http://www.ge.infn.it/geant4/lowE/index.html

// -------------------------------------------------------------------

#include "globals.hh"
#include "G4ProcessTest.hh"
#include "G4VProcess.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4VParticleChange.hh"
#include "G4ParticleChange.hh"
#include "SystemOfUnits.h"

G4ProcessTest::G4ProcessTest() : process(0), ioni(0), brem(0)
{ }

G4ProcessTest:: ~G4ProcessTest()
{
  delete process;
  process = 0;
  delete ioni;
  ioni = 0;
  delete brem;
  brem = 0;
 }

void G4ProcessTest::init(const G4String& type, G4bool isPolarised = false)
{ 
  process = createProcess(const G4String& type, G4bool isPolarised = false);
  G4cout << process->GetName() << " created" << G4endl;
  ioni = createElectronIonisation(const G4String& type);
  if (ioni != 0)  G4cout << ioni->GetName() << " created" << G4endl;
  brem = createBremsstrahlung(const G4String& type);
  if (brem != 0)  G4cout << brem->GetName() << " created" << G4endl;
}
 
void G4ProcessTest::postStepTest(const G4Track& track,
				 const G4Step& step) const
{
  G4ParticleChange* particleChange = (G4ParticleChange*) process->PostStepDoIt(track, step);
  
  // Primary physical quantities 
  
  G4double energyChange = particleChange->GetEnergyChange();
  G4double initEnergy = track.GetKineticEnergy();
  G4double dedx = initEnergy - energyChange ;
  G4double dedxNow = dedx / (step.GetStepLength());
  
  G4ThreeVector eChange = particleChange->CalcMomentum(energyChange,
						       (*particleChange->GetMomentumChange()),
						       particleChange->GetMassChange());
  
  G4double pxChange  = eChange.x();
  G4double pyChange  = eChange.y();
  G4double pzChange  = eChange.z();
  G4double pChange   = sqrt(pxChange*pxChange + pyChange*pyChange + pzChange*pzChange);
  
  G4double xChange = particleChange->GetPositionChange()->x();
  G4double yChange = particleChange->GetPositionChange()->y();
  G4double zChange = particleChange->GetPositionChange()->z();
  G4double thetaChange = particleChange->GetPositionChange()->theta();
  
  G4cout << "---- Primary after the step ---- " << G4endl;
  
  //      G4cout << "Position (x,y,z) = " 
  //	     << xChange << "  " 
  //	     << yChange << "   " 
  //	     << zChange << "   " 
  //	     << G4endl;
  
  G4cout << "---- Energy: " << energyChange/MeV << " MeV,  " 
	 << "(px,py,pz): ("
	 << pxChange/MeV << ","
	 << pyChange/MeV << "," 
	 << pzChange/MeV << ") MeV"
	 << G4endl;
  
  G4cout << "---- Energy loss (dE) = " << dedx/keV << " keV" << G4endl;
  //      G4cout << "Stopping power (dE/dx)=" << dedxNow << G4endl;
  
  
  G4int nElectrons = 0;
  G4int nPositrons = 0;
  G4int nPhotons = 0;
  
  for (G4int i = 0; i < (particleChange->GetNumberOfSecondaries()); i++) 
    {
      // The following two items should be filled per event, not
      // per secondary; filled here just for convenience, to avoid
      // complicated logic to dump ntuple when there are no secondaries
      
      G4Track* finalParticle = particleChange->GetSecondary(i) ;
      
      G4double e    = finalParticle->GetTotalEnergy();
      G4double eKin = finalParticle->GetKineticEnergy();
      G4double px   = (finalParticle->GetMomentum()).x();
      G4double py   = (finalParticle->GetMomentum()).y();
      G4double pz   = (finalParticle->GetMomentum()).z();
      G4double theta   = (finalParticle->GetMomentum()).theta();
      G4double p   = sqrt(px*px+py*py+pz*pz);
      
      if (e > initEnergy)
	{
	  G4cout << "WARNING: eFinal > eInit " << G4endl;
	  //	     << e
	  //		     << " > " initEnergy 
	  
	}
      
      G4String particleName = finalParticle->GetDefinition()->GetParticleName();
      G4cout  << "==== Final " 
	      <<  particleName  <<  " "  
	      << "energy: " <<  e/MeV  <<  " MeV,  " 
	      << "eKin: " <<  eKin/MeV  <<  " MeV, " 
	      << "(px,py,pz): ("
	      <<  px/MeV  <<  "," 
	      <<  py/MeV  <<  ","
	      <<  pz/MeV  << ") MeV "
	      <<  G4endl;   
      G4int partType = 0;
      
      if (particleName == "e-") 
	{
	  partType = 1;
	  nElectrons++;
	}
      else if (particleName == "e+") 
	{
	  partType = 2;
	  nPositrons++;
	}
      else if (particleName == "gamma") 
	{
	  partType = 3;
	  nPhotons++;
	}
      
      delete particleChange->GetSecondary(i);
    }
  
  particleChange->Clear();
}

void G4ProcessTest::alongStepTest(const G4Track& track,
				  const G4Step& step) const 
{
  // To be implemented
}
