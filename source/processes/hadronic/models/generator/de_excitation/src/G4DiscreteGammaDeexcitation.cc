// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// -------------------------------------------------------------------
//      GEANT 4 class file 
//
//      For information related to this code contact:
//      CERN, IT Division, ASD group
//      CERN, Geneva, Switzerland
//
//      File name:     G4DiscreteGammaDeexcitation
//
//      Author:        Maria Grazia Pia (pia@genova.infn.it)
// 
//      Creation date: 23 October 1998
//
//      Modifications: 
//      
// -------------------------------------------------------------------

#include "G4DiscreteGammaDeexcitation.hh"
#include "G4DiscreteGammaTransition.hh"
#include "G4NuclearLevelManager.hh"


G4DiscreteGammaDeexcitation::G4DiscreteGammaDeexcitation(): _Z(0),_A(0)
{
  _tolerance = 0.1 * MeV;
}


G4DiscreteGammaDeexcitation::~G4DiscreteGammaDeexcitation() {}


G4VGammaTransition* G4DiscreteGammaDeexcitation::CreateTransition()
{
  G4Fragment nucleus = GetNucleus();
  G4int A = nucleus.GetA();
  G4int Z = nucleus.GetZ();

  if (_levelManager.IsValid(Z,A)) 
    { 
      if (_verbose > 1)
	G4cout 
	  << "G4DiscreteGammaDeexcitation::CreateTransition - (A,Z) is valid " 
	  << endl;
      
      if (_A != A || _Z != Z) 
	{
	  _levelManager.SetNucleus(Z,A);
	  _A = A;
	  _Z = Z;
	}
      
      G4double excitation = nucleus.GetExcitationEnergy();
      //      const G4NuclearLevel* level =_levelManager.NearestLevel(excitation, _tolerance);
      const G4NuclearLevel* level =_levelManager.NearestLevel(excitation);

      if (level != 0)  
	{ 
	  if (_verbose > 0)
	    G4cout 
	      << "G4DiscreteGammaDeexcitation::CreateTransition - Created from level energy " 
	      << level->Energy() << ", excitation is " 
	      << excitation << endl;
	  return new G4DiscreteGammaTransition(*level);
	}
      else 
	{ 
	  if (_verbose > 0)
	    G4cout 
	      << "G4DiscreteGammaDeexcitation::CreateTransition - No transition created from "
	      << excitation << " within tolerance " << _tolerance << endl;
	  
	  return 0; 
	}
    }
  else return 0;
}


G4bool G4DiscreteGammaDeexcitation::CanDoTransition() const
{

  G4bool canDo = true;

  if (_transition == 0) 
    {
      canDo = false;

      if (_verbose > 0)
	G4cout 
	  << "G4DiscreteGammaDeexcitation::CanDoTransition - Null transition " 
	  << endl;
    }

  G4Fragment nucleus = GetNucleus();


  G4double A = nucleus.GetA();
  G4double Z = nucleus.GetZ();
  if (A <2 || Z<3 || Z>92)
    {
      canDo = false;
      if (_verbose > 0) 
	G4cout 
	  << "G4DiscreteGammaDeexcitation::CanDoTransition - n/p/H/>U"
	  << endl;
    }

  G4double excitation = nucleus.GetExcitationEnergy();
  if (excitation <= 0.) 
    {
      canDo = false;
      if (_verbose > 0) 
	G4cout 
	  << "G4DiscreteGammaDeexcitation::CanDoTransition -  Excitation <= 0" 
	  << endl;
    }

  if (excitation > _levelManager.MaxLevelEnergy() + _tolerance) canDo = false;
  if (excitation < _levelManager.MinLevelEnergy() - _tolerance) canDo = false;  
  // The following is a protection to avoid looping in case of elements with very low
  // ensdf levels
  if (excitation < _levelManager.MinLevelEnergy() * 0.9) canDo = false;  
  
  if (_verbose > 0)
    {
      G4cout << "G4DiscreteGammaDeexcitation::CanDoTransition -  Excitation " 
	     << excitation << ", Min-Max are " 
	     << _levelManager.MinLevelEnergy() << " "
	     << _levelManager.MaxLevelEnergy() << endl;
    }
  
  if (canDo)
    { if (_verbose > 0) 
      G4cout <<"G4DiscreteGammaDeexcitation::CanDoTransition - CanDo" << endl; }
  //  else
  //    { 
  //      delete _transition;
  //      _transition = 0;
  //    }

  return canDo;
      
}

