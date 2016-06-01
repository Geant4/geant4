// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4DistributionGenerator.hh,v 2.1 1998/07/12 03:07:08 urbi Exp $
// GEANT4 tag $Name: geant4-00 $
//
// -------------------------------------------------------------------
//      GEANT 4 class file --- Copyright CERN 1998
//      CERN Geneva Switzerland
//
//      For information related to this code contact:
//      CERN, IT Division, ASD group
//
//      File name:     G4DistributionGenerator.hh 
//
//      Author:        Maria Grazia Pia (pia@genova.infn.it)
// 
//      Creation date: 12 May 1998
//
//      Modifications: 
// -------------------------------------------------------------------

#ifndef G4DISTRIBUTIONGENERATOR_HH
#define G4DISTRIBUTIONGENERATOR_HH

#include <rw/tvordvec.h>

#include "globals.hh"
#include "Randomize.hh"

class G4DistributionGenerator
{  

public:

  // Constructor
  G4DistributionGenerator(RWTValOrderedVector<G4double>& x,
  			  RWTValOrderedVector<G4double>& values);
  G4DistributionGenerator();

  // Destructor
  ~G4DistributionGenerator();

  G4double Generate(G4double ranflat);

private:

  // Hide assignment operator as private 
  //  G4DistributionGenerator& operator=(const G4DistributionGenerator &right);

  // Copy constructor
  //  G4DistributionGenerator(const G4DistributionGenerator& );

  RWTValOrderedVector<G4double> _x;
  RWTValOrderedVector<G4double> _cumProb;

};
 
#endif
