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
// $Id: G4ProcessTestAnalysis.cc,v 1.3 2001-11-08 23:31:24 pia Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Author:  A. Pfeiffer (Andreas.Pfeiffer@cern.ch) 
//         (copy of his UserAnalyser class)
//
// History:
// -----------
//  5 Nov 2001   MGP        Implementation of process test analysis
//
// -------------------------------------------------------------------

#include "G4ProcessTestAnalysis.hh"
#include "globals.hh"
#include "G4ParticleChange.hh"
#include "G4Track.hh"
#include "PhysicalConstants.h"
#include "SystemOfUnits.h"

G4ProcessTestAnalysis* G4ProcessTestAnalysis::instance = 0;

G4ProcessTestAnalysis::G4ProcessTestAnalysis()
{
  histoManager = createIHistoManager();
  ntFactory = Lizard::createNTupleFactory();
}

G4ProcessTestAnalysis::~G4ProcessTestAnalysis()
{ 
  delete ntFactory;
  ntFactory = 0;

  //  G4std::map< G4String,IHistogram1D*,G4std::less<G4String> >::iterator pos1;
  //  for (pos1 = histo1D.begin(); pos1 != histo1D.end(); ++pos1)
  //    {
  //    IHistogram* h = pos1->second;
  //   delete h;
  //  }

  delete histoManager;
  histoManager = 0;
}

void G4ProcessTestAnalysis::book(const G4String& storeName)
{
  const char* nameStore = (storeName + ".hbook").c_str();
  histoManager->selectStore(nameStore);

  // Book histograms

  IHistogram1D* hNSec = histoManager->create1D("1","Number of secondaries", 10,0.,10.);
  //histo1D["nSec"] = hNSec;

  IHistogram1D* hDeposit = histoManager->create1D("2","Local energy deposit", 100,0.,10.);
  //histo1D["eDeposit"] = hDeposit;

  IHistogram1D* hEKin = histoManager->create1D("3","Kinetic Energy", 100,0.,10.);
  //histo1D["eKin"] = hEKin;

  IHistogram1D* hTheta = histoManager->create1D("4","Theta", 100,0.,pi);
  //histo1D["theta"] = hTheta;

  IHistogram1D* hPhi = histoManager->create1D("5","Phi", 100,-pi,pi);
  //histo1D["phi"] = hPhi;

 // Book ntuples

  //  G4String ntFileName = storeName + "1" + ".hbook::1";
  G4String ntFileName = storeName + ".hbook::1";
  const char* name = ntFileName.c_str();
  Lizard::NTuple* ntuple1 = ntFactory->createC(name);

  //  Add and bind the attributes to the general final state nTuple

  if( !( ntuple1->addAndBind( "eprimary"  , initialEnergy) &&
	 ntuple1->addAndBind( "pxch"      , pxChange     ) &&
	 ntuple1->addAndBind( "pych"      , pyChange     ) &&
	 ntuple1->addAndBind( "pzch"      , pzChange     ) &&
	 ntuple1->addAndBind( "thetach"   , thetaChange  ) &&
	 ntuple1->addAndBind( "phich"     , phiChange    ) &&
	 ntuple1->addAndBind( "deposit"   , eDeposit     ) &&
	 ntuple1->addAndBind( "eminus"    , nElectrons   ) &&
	 ntuple1->addAndBind( "eplus"     , nPositrons   ) &&
	 ntuple1->addAndBind( "photons"   , nPhotons     ) ) ) 
    {
      G4cerr << "Error: unable to add attribute to nTuple1." << G4endl;
      // Must be cleaned up properly before any exit.
      delete ntuple1;
      G4Exception("Could not addAndBind ntuple1");
    }
  ntuples["primary"] = ntuple1;

  //  ntFileName = storeName + "2" + ".hbook::2";
  ntFileName = storeName + ".hbook::2";
  name = ntFileName.c_str();
  Lizard::NTuple* ntuple2 = ntFactory->createC(name);

  //  Add and bind the attributes to the secondary nTuple
  if ( !( ntuple2->addAndBind( "e0"      , initialEnergy) &&
	  ntuple2->addAndBind( "px"	 , px        ) &&
	  ntuple2->addAndBind( "py"	 , py        ) &&
	  ntuple2->addAndBind( "pz"	 , pz        ) &&
	  ntuple2->addAndBind( "p" 	 , p         ) &&
	  ntuple2->addAndBind( "e" 	 , e         ) &&
	  ntuple2->addAndBind( "ekin"    , eKin      ) &&
	  ntuple2->addAndBind( "theta"   , theta     ) &&
	  ntuple2->addAndBind( "phi"     , phi       ) &&
	  ntuple2->addAndBind( "type"    , partType  )  ) )
    {
      G4cerr << "Error: unable to add attribute to nTuple2" << G4endl;
      // Must be cleaned up properly before any exit
      delete ntuple2;
      G4Exception("Could not addAndBind ntuple2");
    }
  ntuples["secondaries"] = ntuple2;
}

void G4ProcessTestAnalysis::finish()
{
  // Because of a Lizard feature, ntuples must be deleted at this stage, 
  // not in the destructor (otherwise the ntuples are not stored)

  G4std::map< G4String,Lizard::NTuple*,G4std::less<G4String> >::iterator pos;
  Lizard::NTuple* ntuple;
  pos = ntuples.find("primary");
  if (pos != ntuples.end()) 
    {
      ntuple = pos->second;
      delete ntuple;
      ntuples["primary"] = 0;
    }

  pos = ntuples.find("secondaries");
  if (pos != ntuples.end()) 
    {
      ntuple = pos->second;
      delete ntuple;
      ntuples["secondary"] = 0;
    }

  histoManager->store("1");
  histoManager->store("2");
  histoManager->store("3");
  histoManager->store("4");
  histoManager->store("5");
}

void G4ProcessTestAnalysis::analyseSecondaries(const G4ParticleChange* particleChange)
{
  G4int nSecondaries = particleChange->GetNumberOfSecondaries();
  for (G4int i = 0; i < nSecondaries; i++) 
    {
      G4Track* finalParticle = particleChange->GetSecondary(i) ;

      // The quantities bound to ntuple2      
      px    = (finalParticle->GetMomentum()).x();
      py    = (finalParticle->GetMomentum()).y();
      pz    = (finalParticle->GetMomentum()).z();
      p     = sqrt(px*px + py*py + pz*pz);
      e     = finalParticle->GetTotalEnergy();
      eKin  = finalParticle->GetKineticEnergy();
      theta = (finalParticle->GetMomentum()).theta();
      phi   = (finalParticle->GetMomentum()).phi();

      partType = 0;
      G4String particleName = finalParticle->GetDefinition()->GetParticleName();

      if (particleName == "e-") partType = 1;
      else if (particleName == "e+") partType = 2;
      else if (particleName == "gamma") partType = 3;

      G4std::map< G4String,Lizard::NTuple*,G4std::less<G4String> >::iterator pos;
      pos = ntuples.find("secondaries");
      if (pos != ntuples.end()) 
	{
	  Lizard::NTuple* ntuple2 = pos->second;
	  ntuple2->addRow();
	}
            
      // Fill histograms
      
      IHistogram1D* h = histoManager->retrieveHisto1D("3");
      h->fill(eKin);
      h = histoManager->retrieveHisto1D("4");
      h->fill(theta);
      h = histoManager->retrieveHisto1D("5");
      h->fill(phi);

      /*      
      G4std::map< G4String,IHistogram1D*,G4std::less<G4String> >::iterator pos1;

      pos1 = histo1D.find("eKin");
      if (pos1 != histo1D.end()) 
	{
	  IHistogram1D* h = pos1->second;
	  h->fill(eKin);
	}

      pos1 = histo1D.find("theta");
      if (pos1 != histo1D.end()) 
	{
	  IHistogram1D* h = pos1->second;
	  h->fill(theta);
	}

      pos1 = histo1D.find("phi");
      if (pos1 != histo1D.end()) 
	{
	  IHistogram1D* h = pos1->second;
	  h->fill(phi);
	}
      */ 
    
      G4cout  << "==== Final " 
	      <<  particleName  <<  " "  
	      << "energy: " <<  e/MeV  <<  " MeV,  " 
	      << "eKin: " <<  eKin/MeV  <<  " MeV, " 
	      << "(px,py,pz): ("
	      <<  px/MeV  <<  "," 
	      <<  py/MeV  <<  ","
	      <<  pz/MeV  << ") MeV "
	      <<  G4endl;   
   }
}

void G4ProcessTestAnalysis::analyseGeneral(const G4Track& track,
					   const G4ParticleChange* particleChange)
{
  // Primary physical quantities 
  
  energyChange = particleChange->GetEnergyChange();
  initialEnergy = track.GetKineticEnergy();
  
  G4ThreeVector eChange = particleChange->CalcMomentum(energyChange,
						       (*particleChange->GetMomentumChange()),
						       particleChange->GetMassChange());
  
  pxChange  = eChange.x();
  pyChange  = eChange.y();
  pzChange  = eChange.z();
  thetaChange = particleChange->GetPositionChange()->theta();
  phiChange = particleChange->GetPositionChange()->phi();
  eDeposit = particleChange->GetLocalEnergyDeposit();

  nElectrons = 0;
  nPositrons = 0;
  nPhotons = 0;

  G4int nSecondaries = particleChange->GetNumberOfSecondaries();

  for (G4int i = 0; i < nSecondaries; i++) 
    {  
      G4Track* finalParticle = particleChange->GetSecondary(i) ;
      G4String particleName = finalParticle->GetDefinition()->GetParticleName();

      if (particleName == "e-") nElectrons++;
      else if (particleName == "e+") nPositrons++;
      else if (particleName == "gamma") nPhotons++;
    }

  G4std::map< G4String,Lizard::NTuple*,G4std::less<G4String> >::iterator pos;
  pos = ntuples.find("primary");
  if (pos != ntuples.end()) 
    {
      Lizard::NTuple* ntuple1 = pos->second;
      ntuple1->addRow();
    }

  // Fill histograms
  IHistogram1D* h = histoManager->retrieveHisto1D("1");
  h ->fill(float(nSecondaries));
  h = histoManager->retrieveHisto1D("2");
  h->fill(eDeposit);

  /*
  G4std::map< G4String,IHistogram1D*,G4std::less<G4String> >::iterator pos1;
  pos1 = histo1D.find("nSec");
  if (pos1 != histo1D.end()) 
    {
       IHistogram1D* h2 = pos1->second;
       h2->fill(float(nSecondaries));
    }
  
  pos1 = histo1D.find("eDeposit");
  if (pos1 != histo1D.end()) 
    {
      IHistogram1D* h2 = pos1->second; 
      h2->fill(eDeposit);
    }
  */
}


G4ProcessTestAnalysis* G4ProcessTestAnalysis::getInstance()
{
  if (instance == 0) instance = new G4ProcessTestAnalysis;
  return instance;
}
