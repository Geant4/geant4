// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4PiMinusStopAbsorption.hh,v 1.1 1999-01-07 16:13:39 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// -------------------------------------------------------------------
//      GEANT 4 class file --- Copyright CERN 1998
//      CERN Geneva Switzerland
//
//      For information related to this code contact:
//      CERN, IT Division, ASD group
//
//      File name:     G4PiMinusStopAbsorption.hh 
//
//      Author:        Maria Grazia Pia (pia@genova.infn.it)
// 
//      Creation date: 12 May 1998
//
//      Modifications: 
//      13 Sep 1998 - Changed DoAbsorption
//
// -------------------------------------------------------------------

#ifndef G4PIMINUSSTOPABSORPTION_HH
#define G4PIMINUSSTOPABSORPTION_HH

#include <rw/tpordvec.h>

#include "globals.hh"
#include "G4DynamicParticle.hh"
#include "G4DynamicParticleVector.hh"
#include "G4PiMinusStopMaterial.hh"
#include "G4ThreeVector.hh"

class G4PiMinusStopAbsorption
{  

public:

  // Constructor
  G4PiMinusStopAbsorption(G4PiMinusStopMaterial* materialAlgo, const G4double Z, const G4double A);

  // Destructor
  ~G4PiMinusStopAbsorption();

  // Return final absorption products
  G4DynamicParticleVector* DoAbsorption();

  // Energy involved in the absorption process
  G4double Energy();

  // Return nucleus recoil momentum
  G4ThreeVector RecoilMomentum();

  // Number of protons in the absorption products
  G4int NProtons();

  // Number of neutrons in the absorption products
  G4int NNeutrons();

  void SetVerboseLevel(G4int level);

private:

  // Hide assignment operator as private 
  G4PiMinusStopAbsorption& operator=(const G4PiMinusStopAbsorption &right);

  // Copy constructor
  G4PiMinusStopAbsorption(const G4PiMinusStopAbsorption& );

  G4PiMinusStopMaterial* _materialAlgo; // owned pointer
  G4DynamicParticleVector* _absorptionProducts;

  G4double _A;
  G4double _Z;

  G4int _level;
};
 
#endif




















