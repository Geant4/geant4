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
//        09 Sep. 2002, Fan Lei  (flei@space.qinetiq.com)
//              Added renormalization to determine whether transition leads to
//              electron or gamma in SelectGamma()
//
//        21 Nov. 2001, Fan Lei (flei@space.qinetiq.com)
//              i) added G4int _nucleusZ initialise it through the constructor
//              ii) modified SelectGamma() to allow the generation of conversion electrons      
//              iii) added #include G4AtomicShells.hh
//      
//        15 April 1999, Alessandro Brunengo (Alessandro.Brunengo@ge.infn.it)
//              Added creation time evaluation for products of evaporation
//      
// -------------------------------------------------------------------

#include "G4DiscreteGammaTransition.hh"
#include "Randomize.hh"
#include "G4RandGeneralTmp.hh"
#include "G4AtomicShells.hh"

G4DiscreteGammaTransition::G4DiscreteGammaTransition(const G4NuclearLevel& level): 
  _gammaEnergy(0.), _level(level), _excitation(0.), _gammaCreationTime(0.)
{ }

G4DiscreteGammaTransition::G4DiscreteGammaTransition(const G4NuclearLevel& level, G4int Z): 
  _nucleusZ(Z), _orbitE(-1), _bondE(0.), _aGamma(true), _icm(false), _gammaEnergy(0.), 
  _level(level), _excitation(0.),  _gammaCreationTime(0.)
{
  _verbose = 0;
}


G4DiscreteGammaTransition::~G4DiscreteGammaTransition() 
{ }


void G4DiscreteGammaTransition::SelectGamma()
{

  _gammaEnergy = 0.;

  G4int nGammas = _level.NumberOfGammas();
  if (nGammas > 0)
    {
      G4double random = G4UniformRand();

      G4int iGamma = 0;
      for(iGamma=0;iGamma < nGammas;iGamma++)
	{
	  if(random <= (_level.GammaCumulativeProbabilities())[iGamma])
	    break;
	}


      // Small correction due to the fact that there are mismatches between 
      // nominal level energies and emitted gamma energies

      G4double eCorrection = _level.Energy() - _excitation;

      _gammaEnergy = (_level.GammaEnergies())[iGamma] - eCorrection;

      //  Warning: the following check is needed to avoid loops:
      //  Due essentially to missing nuclear levels in data files, it is
      //  possible that _gammaEnergy is so low as the nucleus doesn't change
      //  its level after the transition.
      //  When such case is found, force the full deexcitation of the nucleus.
      //
      //    NOTE: you should force the transition to the next lower level,
      //          but this change needs a more complex revision of actual design.
      //          I leave this for a later revision.

      if (_gammaEnergy < _level.Energy()*10e-5) _gammaEnergy = _excitation;
      // now decide whether Internal Coversion electron should be emitted instead
      if (_icm) {
	random = G4UniformRand() ;
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
	    if (_verbose > 0)
	      G4cout << "G4DiscreteGammaTransition: _nucleusZ = " <<_nucleusZ 
		     << " , iShell = " << iShell  
		     << " , Shell binding energy = " << G4AtomicShells::GetBindingEnergy(_nucleusZ, iShell) / keV
		     << " keV " << G4endl;
	    _bondE = G4AtomicShells::GetBindingEnergy(_nucleusZ, iShell);
	    _gammaEnergy = _gammaEnergy - _bondE;
	    _orbitE = iShell;	  
	    _aGamma = false ;   // emitted is not a gamma now 
	  }
      }
    
      G4double tau = _level.HalfLife() / std::log(2.0);

      G4double tMin = 0;
      G4double tMax = 10.0 * tau;
      //  Original code, not very efficent
      //  G4int nBins = 200;
      //G4double sampleArray[200];

      //  for(G4int i = 0;i<nBins;i++)
      //{
      //  G4double t = tMin + ((tMax-tMin)/nBins)*i;
      //  sampleArray[i] = (std::exp(-t/tau))/tau;
      // }

      //  G4RandGeneralTmp randGeneral(sampleArray, nBins);
      //G4double random = randGeneral.shoot();
  
      //_gammaCreationTime = tMin + (tMax - tMin) * random;

      // new code by Fan Lei
      //
      if (tau != 0 ) 
      {
	  random = G4UniformRand() ;
	  _gammaCreationTime = -(std::log(random*(std::exp(-tMax/tau) - std::exp(-tMin/tau)) + 
					std::exp(-tMin/tau)));
	  //  if(_verbose > 10)
	  //    G4cout << "*---*---* G4DiscreteTransition: _gammaCreationTime = "
	  //	   << _gammaCreationTime/second << G4endl;
       } else { _gammaCreationTime=0.; }
    }
  return;
}


//G4bool G4DiscreteGammaTransition::IsAGamma()
//{
//  return _aGamma;
//}


G4double G4DiscreteGammaTransition::GetGammaEnergy()
{
  return _gammaEnergy;
}

G4double G4DiscreteGammaTransition::GetGammaCreationTime()
{
  return _gammaCreationTime;
}

void G4DiscreteGammaTransition::SetEnergyFrom(const G4double energy)
{
  _excitation = energy;
  return;
}






