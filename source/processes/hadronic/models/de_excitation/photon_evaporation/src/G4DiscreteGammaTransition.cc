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
// $Id$
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
//              ii) modified SelectGamma() to allow the generation of conversion electrons    
//              iii) added #include G4AtomicShells.hh
//      
//        09 Sep. 2002, Fan Lei  (flei@space.qinetiq.com)
//              Added renormalization to determine whether transition leads to
//              electron or gamma in SelectGamma()
//
//        19 April 2010, J. M. Quesada. 
//              Corrections added for taking into account mismatch between tabulated 
//              gamma energies and level energy differences (fake photons eliminated) 
//
//        9 May 2010, V.Ivanchenko
//              Removed unphysical corretions of gamma energy; fixed default particle 
//              as gamma; do not subtract bounding energy in case of electron emmision
//
//		  03 November 2011, L. Desorgher
//				Extend the use of the code for Z>100 by not calling G4AtomicShells::GetBindingEnergy for Z>100
//				For Z>100 the binding energy is set to 0, the atomic relaxation is not simulated in G4RadDecay
//
// -------------------------------------------------------------------

#include "G4DiscreteGammaTransition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4RandGeneralTmp.hh"
#include "G4AtomicShells.hh"
#include "G4NuclearLevel.hh"
#include "G4NuclearLevelStore.hh"
#include "G4Pow.hh"

G4DiscreteGammaTransition::G4DiscreteGammaTransition(const G4NuclearLevel& level, G4int Z, G4int A): 
  _nucleusZ(Z), _orbitE(-1), _bondE(0.), _aGamma(true), _icm(false), _gammaEnergy(0.), 
  _level(level), _excitation(0.),  _gammaCreationTime(0.),_A(A),_Z(Z)
{
  _levelManager = 0;
  _verbose = 0;
  //JMQ: added tolerence in the mismatch
  //VI:  increased tolerence 
  _tolerance = 10*CLHEP::keV;
}

G4DiscreteGammaTransition::~G4DiscreteGammaTransition() 
{}

void G4DiscreteGammaTransition::SelectGamma()
{
  // default gamma 
  _aGamma = true;    
  _gammaEnergy = 0.;
  
  G4int nGammas = _level.NumberOfGammas();
  if (nGammas > 0)
    {
      G4int iGamma = 0;
      if(1 < nGammas) {
	G4double random = G4UniformRand();
      
	//G4cout << "G4DiscreteGammaTransition::SelectGamma  N= " 
	//       << nGammas << " rand= " << random << G4endl;
	for(iGamma=0; iGamma<nGammas; ++iGamma)
	  {
	    //G4cout << iGamma << "  prob= " 
	    //	   << (_level.GammaCumulativeProbabilities())[iGamma] << G4endl;
	    if(random <= (_level.GammaCumulativeProbabilities())[iGamma])
	      { break; }
	  }
      }
      /*     
      G4cout << "Elevel(MeV)= " << _level.Energy()/MeV
	     << " Etran(MeV)= " << (_level.GammaEnergies())[iGamma]/MeV
	     << " Eexc(MeV)= " << _excitation/MeV << G4endl;
      */

      // VI: do not apply correction here in order do not make 
      //     double correction
      //G4double eCorrection = _level.Energy() - _excitation;      
      //_gammaEnergy = (_level.GammaEnergies())[iGamma] - eCorrection;
      _gammaEnergy = (_level.GammaEnergies())[iGamma];
            
      //JMQ: 
      //1)If chosen gamma energy is close enough to excitation energy, 
      //  the later is used instead for gamma dacey to gs (it guarantees 
      //  energy conservation)
      //2)For energy conservation, level energy differences instead of  
      //  tabulated gamma energies must be used (origin of final fake photons)
      
      // VI: remove fake photons - applied only for the last transition
      //     do not applied on each transition
      if(std::fabs(_excitation - _gammaEnergy) < _tolerance) { 
	_gammaEnergy =_excitation;
      }

      //  JMQ: Warning: the following check is needed to avoid loops:
      //  Due essentially to missing nuclear levels in data files, it is
      //  possible that _gammaEnergy is so low as the nucleus doesn't change
      //  its level after the transition.
      //  When such case is found, force the full deexcitation of the nucleus.
      //
      //    NOTE: you should force the transition to the next lower level,
      //          but this change needs a more complex revision of actual 
      //          design.
      //          I leave this for a later revision.

      // VI: the check has no sence and we make this very simple
      if (_gammaEnergy < _tolerance) { 
	_gammaEnergy = _excitation; 
      }

      //G4cout << "G4DiscreteGammaTransition::SelectGamma: " << _gammaEnergy 
      //	     << " _icm: " << _icm << G4endl;

      // now decide whether Internal Coversion electron should be emitted instead
      if (_icm) {
	G4double random = G4UniformRand();
	if ( random <= (_level.TotalConvertionProbabilities())[iGamma]
	     *(_level.GammaWeights())[iGamma]
	     /((_level.TotalConvertionProbabilities())[iGamma]*(_level.GammaWeights())[iGamma]
	       +(_level.GammaWeights())[iGamma])) 
	  {
	    G4int iShell = 9;
	    random = G4UniformRand() ;
	    if ( random <= (_level.KConvertionProbabilities())[iGamma]) 
	      { iShell = 0;}
	    else if ( random <= (_level.L1ConvertionProbabilities())[iGamma]) 
	      { iShell = 1;}
	    else if ( random <= (_level.L2ConvertionProbabilities())[iGamma]) 
	      { iShell = 2;}
	    else if ( random <= (_level.L3ConvertionProbabilities())[iGamma]) 
	      { iShell = 3;}	
	    else if ( random <= (_level.M1ConvertionProbabilities())[iGamma]) 
	      { iShell = 4;}
	    else if ( random <= (_level.M2ConvertionProbabilities())[iGamma]) 
	      { iShell = 5;}
	    else if ( random <= (_level.M3ConvertionProbabilities())[iGamma]) 
	      { iShell = 6;}
	    else if ( random <= (_level.M4ConvertionProbabilities())[iGamma]) 
	      { iShell = 7;}
	    else if ( random <= (_level.M5ConvertionProbabilities())[iGamma]) 
	      { iShell = 8;}
	    // the following is needed to match the ishell to that used in  G4AtomicShells
	    if ( iShell == 9) {
	      if ( (_nucleusZ < 28) && (_nucleusZ > 20)) {
		iShell--;
	      } else if ( _nucleusZ == 20 || _nucleusZ == 19 ) {
		iShell = iShell -2;
	      }
	    }
	    //L.Desorgher 02/11/2011
	    //Atomic shell information is available in Geant4 only up top Z=100
	    //To extend the photo evaporation code to Z>100  the call 
	    // to G4AtomicShells::GetBindingEnergy should be forbidden for Z>100
	    _bondE = 0.;
	    if (_nucleusZ <=100)
	    _bondE = G4AtomicShells::GetBindingEnergy(_nucleusZ, iShell);
	    if (_verbose > 0) {
	      G4cout << "G4DiscreteGammaTransition: _nucleusZ = " <<_nucleusZ 
		     << " , iShell = " << iShell  
		     << " , Shell binding energy = " << _bondE/keV
		     << " keV " << G4endl;
	    }

	    // 09.05.2010 VI : it is an error - cannot subtract bond energy from 
	    //                 transition energy here
	    //_gammaEnergy = _gammaEnergy - _bondE; 
	    //G4cout << "_gammaEnergy = " << _gammaEnergy << G4endl;

	    _orbitE = iShell;	  
	    _aGamma = false ;   // emitted is not a gamma now 
	  }
      }
      
      G4double tau = _level.HalfLife() / G4Pow::GetInstance()->logZ(2);

      //09.05.2010 VI rewrite samling of decay time 
      //              assuming ordinary exponential low
      _gammaCreationTime = 0.;      
      if(tau > 0.0) {  _gammaCreationTime = -tau*std::log(G4UniformRand()); }

    }
  return;
}

G4double G4DiscreteGammaTransition::GetGammaEnergy()
{
  return _gammaEnergy;
}

G4double G4DiscreteGammaTransition::GetGammaCreationTime()
{
  return _gammaCreationTime;
}

void G4DiscreteGammaTransition::SetEnergyFrom(G4double energy)
{
  _excitation = energy;
}






