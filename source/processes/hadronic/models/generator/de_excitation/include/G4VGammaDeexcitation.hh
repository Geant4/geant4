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
//      File name:     G4VGammaDeexcitation
//
//      Author:        Maria Grazia Pia (pia@genova.infn.it)
// 
//      Creation date: 23 October 1998
//
//      Modifications: 
//      
// -------------------------------------------------------------------

#ifndef G4VGAMMADEEXCITATION_HH
#define G4VGAMMADEEXCITATION_HH

#include "globals.hh"
#include "G4VGammaTransition.hh"
#include "G4Fragment.hh"
#include "G4FragmentVector.hh"

class G4VGammaDeexcitation 
{

public:

  G4VGammaDeexcitation();

  virtual ~G4VGammaDeexcitation();
  
  virtual G4VGammaTransition* CreateTransition() = 0;
  virtual G4bool CanDoTransition() const = 0;

  // Single gamma transition
  virtual G4FragmentVector* DoTransition();

  // Chain of gamma transitions
  virtual G4FragmentVector* DoChain();

  virtual G4Fragment* GenerateGamma();

  virtual const G4Fragment& GetNucleus() const;

  virtual void SetNucleus(const G4Fragment& nucleus);

  virtual void SetVerboseLevel(G4int verbose);
  

protected:

  void Initialize();
  void UpdateNucleus(const G4Fragment* gamma);
  void Update(const G4Fragment* gamma);

  G4VGammaTransition* _transition;  // Owned pointer
  G4int _verbose;

private:  

  G4Fragment _nucleus;

  G4VGammaDeexcitation(const G4VGammaDeexcitation &right);

  const G4VGammaDeexcitation& operator=(const G4VGammaDeexcitation &right);
  G4bool operator==(const G4VGammaDeexcitation &right) const;
  G4bool operator!=(const G4VGammaDeexcitation &right) const;
  
};

#endif


















