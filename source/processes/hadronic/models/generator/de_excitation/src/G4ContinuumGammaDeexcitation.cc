// This code implementation is the intellectual property of
// the GEANT4 collaboration.
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
//      File name:     G4ContinuumGammaDeexcitation
//
//      Authors:       Carlo Dallapiccola (dallapiccola@umdhep.umd.edu)
//                     Maria Grazia Pia (pia@genova.infn.it)
// 
//      Creation date: 23 October 1998
//
//      Modifications: 
//      
// -------------------------------------------------------------------
//
// Class G4ContinuumGammaDeexcitation.cc
//
//       Concrete class derived from G4VGammaDeexcitation
//
//
#include "G4ContinuumGammaDeexcitation.hh"

#include "G4Gamma.hh"
#include "G4ContinuumGammaTransition.hh"
#include "G4NuclearLevelManager.hh"
#include "G4Fragment.hh"      
#include "G4ConstantLevelDensityParameter.hh"

//
// Constructor
//

G4ContinuumGammaDeexcitation::G4ContinuumGammaDeexcitation(): _nucleusZ(0), _nucleusA(0)
{ }


G4ContinuumGammaDeexcitation::~G4ContinuumGammaDeexcitation() 
{ }


G4VGammaTransition* G4ContinuumGammaDeexcitation::CreateTransition()
{
  G4Fragment nucleus = GetNucleus();
  G4int Z = nucleus.GetZ();
  G4int A = nucleus.GetA();
  G4double excitation = nucleus.GetExcitationEnergy();

  if (_nucleusA != A || _nucleusZ != Z)
    {
      _levelManager.SetNucleus(Z,A);
      _nucleusA = A;
      _nucleusZ = Z;
    }

  if (_verbose > 1)
    G4cout << "G4ContinuumGammaDeexcitation::CreateTransition - Created" << G4endl;

  return new G4ContinuumGammaTransition(_levelManager,Z,A,excitation,_verbose );
}
    

G4bool G4ContinuumGammaDeexcitation::CanDoTransition() const
{
 G4bool canDo = true;

  if (_transition == 0) 
    {
      canDo = false;

      if (_verbose > 0)
	G4cout 
	  << "G4ContinuumGammaDeexcitation::CanDoTransition - Null transition "
	  << G4endl;
    }

  G4Fragment nucleus = GetNucleus();
  G4double excitation = nucleus.GetExcitationEnergy();

  G4double A = nucleus.GetA();
  G4double Z = nucleus.GetZ();
  if (A <2 || Z<3)
    {
      canDo = false;
      if (_verbose > 0) 
	G4cout 
	  << "G4ContinuumGammaDeexcitation::CanDoTransition - n/p/H"
	  << G4endl;
    }



  if (excitation <= 0.) 
    {
      canDo = false;
      if (_verbose > 0) 
	G4cout 
	  << "G4ContinuumGammaDeexcitation::CanDoTransition -  Excitation <= 0"
	  << G4endl;
    }

  if (excitation <= _levelManager.MaxLevelEnergy()) 
    {
      canDo = false;  
      if (_verbose > 0)
      G4cout << "G4ContinuumGammaDeexcitation::CanDoTransition -  Excitation " 
	     << excitation << " below max discrete level " 
	     << _levelManager.MaxLevelEnergy() << G4endl;
    }
  
  if (canDo)
    { if (_verbose > 1) 
      G4cout <<"G4ContinuumGammaDeexcitation::CanDoTransition - CanDo" 
	     << G4endl; 
    }

  return canDo;

}



