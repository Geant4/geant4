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
// * authors in the GEANT4 collaboration.                             *
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
//      File name:     G4VGammaDeexcitation
//
//      Author:        Maria Grazia Pia (pia@genova.infn.it)
// 
//      Creation date: 23 October 1998
//
//      Modifications: 
//      
//        15 April 1999, Alessandro Brunengo (Alessandro.Brunengo@ge.infn.it)
//              Added creation time evaluation for products of evaporation
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
  void Update();

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


















