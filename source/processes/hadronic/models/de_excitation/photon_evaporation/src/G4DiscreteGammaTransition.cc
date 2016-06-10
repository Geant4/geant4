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
// $Id: G4DiscreteGammaTransition.cc 88987 2015-03-17 10:39:50Z gcosmo $
//
// -------------------------------------------------------------------
//      GEANT 4 class file 
//
//      CERN, Geneva, Switzerland
//
//      File name:     G4DiscreteGammaTransition
//
//      Author:        Maria Grazia Pia   (pia@genova.infn.it)
// 
//      Creation date: 23 October 1998
//
//      Modifications:
//        15 April 1999, Alessandro Brunengo (Alessandro.Brunengo@ge.infn.it)
//              Added creation time evaluation for products of evaporation
//
//        21 Nov. 2001, Fan Lei (flei@space.qinetiq.com)
//              i) added G4int _nucleusZ initialise it through the constructor
//              ii) modified SelectGamma() to allow the generation of 
//                  conversion e-
//              iii) added #include G4AtomicShells.hh
//      
//        09 Sep. 2002, Fan Lei  (flei@space.qinetiq.com)
//              Added renormalization to determine whether transition leads to
//              electron or gamma in SelectGamma()
//
//        19 April 2010, J. M. Quesada. 
//              Corrections added for taking into account mismatch between 
//              tabulated gamma energies and level energy differences 
//              (fake photons eliminated) 
//
//        9 May 2010, V.Ivanchenko
//              Removed unphysical corretions of gamma energy; fixed default 
//              particle as gamma; do not subtract bounding energy in case of 
//              electron emmision
//
//	  3 November 2011, L. Desorgher
//		Extend the use of the code for Z>100 by not calling 
//              G4AtomicShells::GetBindingEnergy for Z>100
//		For Z>100 the binding energy is set to 0 and atomic 
//              relaxation is not simulated in G4RadDecay
//
// -------------------------------------------------------------------

#include "G4DiscreteGammaTransition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4AtomicShells.hh"
#include "G4NuclearLevel.hh"
#include "G4NuclearLevelStore.hh"
#include "G4Pow.hh"
#include "G4Log.hh"

static const G4double tolerance = 10*CLHEP::keV;

G4DiscreteGammaTransition::G4DiscreteGammaTransition(
  const G4NuclearLevel* level, G4int Z, G4int verb)
  : orbitE(-1), bondE(0.),  gammaEnergy(0.),  excitation(0.), 
    gammaCreationTime(0.), aGamma(true), icm(false)
{
  verbose = verb;
  Update(level, Z);
}

G4DiscreteGammaTransition::~G4DiscreteGammaTransition() 
{}

void G4DiscreteGammaTransition::SelectGamma()
{
  // default gamma 
  aGamma = true;    
  gammaEnergy = 0.;
  
  G4int nGammas = aLevel->NumberOfGammas();
  if (nGammas > 0) {
    G4int iGamma = 0;
    if(1 < nGammas) {
      G4double random = G4UniformRand();
      
      //G4cout << "G4DiscreteGammaTransition::SelectGamma  N= " 
      //	     << nGammas << " rand= " << random << G4endl;
      for(iGamma=0; iGamma<nGammas; ++iGamma) {
	//G4cout << iGamma << "  prob= " 
	//   << (aLevel->GammaCumulativeProbabilities())[iGamma] << G4endl;
	if(random <= (aLevel->GammaCumulativeProbabilities())[iGamma])
	  { break; }
      }
    }
    /*
    G4cout << "Elevel(MeV)= " << aLevel->Energy()
	   << " Etran(MeV)= " << (aLevel->GammaEnergies())[iGamma]
	   << " Eexc(MeV)= " << excitation << G4endl;
    */
    // VI 2014: initial excitation energy may be not exactly energy of the level
    //          final excitation energy is always energy of some level or zero
    //          transition to the ground state should be always equal to
    //          the excitation energy 
    gammaEnergy = (aLevel->GammaEnergies())[iGamma] 
      + excitation - aLevel->Energy();

    // this check is needed to remove cases when nucleaus is left in 
    // slightly excited state which will require very low energy
    // gamma emission 
    if(excitation <= gammaEnergy + tolerance) { gammaEnergy = excitation; }
    //JMQ: 
    //1)If chosen gamma energy is close enough to excitation energy, 
    //  the later is used instead for gamma dacey to gs (it guarantees 
    //  energy conservation)
    //2)For energy conservation, level energy differences instead of  
    //  tabulated gamma energies must be used (origin of final fake photons)
      
    // VI: remove fake photons - applied only for the last transition
    //     do not applied on each transition
    //if(std::fabs(excitation - gammaEnergy) < tolerance) { 
    //  gammaEnergy = excitation;
    //}

    //  JMQ: Warning: the following check is needed to avoid loops:
    //  Due essentially to missing nuclear levels in data files, it is
    //  possible that gammaEnergy is so low as the nucleus doesn't change
    //  its level after the transition.
    //  When such case is found, force the full deexcitation of the nucleus.
    //
    //    NOTE: you should force the transition to the next lower level,
    //          but this change needs a more complex revision of actual 
    //          design.
    //          I leave this for a later revision.

    // VI: the check is needed to remove very low-energy gamma
    if (gammaEnergy < tolerance) { gammaEnergy = excitation; }
    /*    
    G4cout << "G4DiscreteGammaTransition::SelectGamma: " << gammaEnergy 
	   << " Eexc= " << excitation
	   << " icm: " << icm << G4endl;
    */
    // now decide whether Internal Coversion electron should be emitted instead
    if (icm) {
      G4double random = G4UniformRand();
      if ( random <= (aLevel->TotalConvertionProbabilities())[iGamma]
	   *(aLevel->GammaWeights())[iGamma]
	   /((aLevel->TotalConvertionProbabilities())[iGamma]
	     *(aLevel->GammaWeights())[iGamma]
	     +(aLevel->GammaWeights())[iGamma])) 
	{
	  G4int iShell = 9;
	  random = G4UniformRand() ;
	  if ( random <= (aLevel->KConvertionProbabilities())[iGamma]) 
	    { iShell = 0;}
	  else if ( random <= (aLevel->L1ConvertionProbabilities())[iGamma]) 
	    { iShell = 1;}
	  else if ( random <= (aLevel->L2ConvertionProbabilities())[iGamma]) 
	    { iShell = 2;}
	  else if ( random <= (aLevel->L3ConvertionProbabilities())[iGamma]) 
	    { iShell = 3;}	
	  else if ( random <= (aLevel->M1ConvertionProbabilities())[iGamma]) 
	    { iShell = 4;}
	  else if ( random <= (aLevel->M2ConvertionProbabilities())[iGamma]) 
	    { iShell = 5;}
	  else if ( random <= (aLevel->M3ConvertionProbabilities())[iGamma]) 
	    { iShell = 6;}
	  else if ( random <= (aLevel->M4ConvertionProbabilities())[iGamma]) 
	    { iShell = 7;}
	  else if ( random <= (aLevel->M5ConvertionProbabilities())[iGamma]) 
	    { iShell = 8;}
	  // the following is needed to match the ishell to that used in  
	  // G4AtomicShells
	  if ( iShell == 9) {
	    if ( (nucleusZ < 28) && (nucleusZ > 20)) {
	      iShell--;
	    } else if ( nucleusZ == 20 || nucleusZ == 19 ) {
	      iShell = iShell -2;
	    }
	  }
	  //L.Desorgher 02/11/2011
	  //Atomic shell information is available in Geant4 only up top Z=100
	  //To extend the photo evaporation code to Z>100  the call 
	  // to G4AtomicShells::GetBindingEnergy should be forbidden for Z>100
	  bondE = 0.;
	  if (nucleusZ <=100) {
	    bondE = G4AtomicShells::GetBindingEnergy(nucleusZ, iShell);
	  }
	  if (verbose > 1) {
	    G4cout << "G4DiscreteGammaTransition: nucleusZ = " <<nucleusZ 
		   << " , iShell = " << iShell  
		   << " , Shell binding energy = " << bondE/keV
		   << " keV " << G4endl;
	  }

	  // last check on energy
	  if(gammaEnergy >  bondE + tolerance) {
	    orbitE = iShell;	  
	    aGamma = false ;   // emitted is not a gamma now 
	    gammaEnergy -= bondE; 
	  }
	  //G4cout << "gammaEnergy = " << gammaEnergy << G4endl;
	}
    }
      
    G4double tau = aLevel->HalfLife() / G4Pow::GetInstance()->logZ(2);

    //09.05.2010 VI rewrite samling of decay time 
    //              assuming ordinary exponential low
    gammaCreationTime = 0.;      
    if(tau > 0.0) {  gammaCreationTime = -tau*G4Log(G4UniformRand()); }
  }
  //G4cout << "G4DiscreteGammaTransition end nGamma= " << nGammas
  //	 << "  Egamma= " << gammaEnergy << G4endl;
}

G4double G4DiscreteGammaTransition::GetGammaEnergy()
{
  return gammaEnergy;
}

G4double G4DiscreteGammaTransition::GetGammaCreationTime()
{
  return gammaCreationTime;
}

void G4DiscreteGammaTransition::SetEnergyFrom(G4double energy)
{
  excitation = energy;
}






