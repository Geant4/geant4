// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4PiMinusStopMaterial.hh,v 1.1 1999-01-07 16:13:40 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// -------------------------------------------------------------------
//      GEANT 4 class file --- Copyright CERN 1998
//      CERN Geneva Switzerland
//
//      For information related to this code contact:
//      CERN, IT Division, ASD group
//
//      File name:     G4PiMinusStopMaterial.hh 
//
//      Author:        Maria Grazia Pia (pia@genova.infn.it)
// 
//      Creation date: 18 May 1998
//
//      Modifications: 
//      13 Sep 1998 - MGP   Modified P4Vector
//      
// -------------------------------------------------------------------

#ifndef G4PIMINUSSTOPMATERIAL_HH
#define G4PIMINUSSTOPMATERIAL_HH

#include <rw/tvordvec.h>
#include <rw/tpordvec.h>
#include <rw/cstring.h>
#include "globals.hh"
#include "G4LorentzVector.hh"
//#include "G4String.hh"
#include "G4DistributionGenerator.hh"
#include "G4ParticleDefinition.hh"


class G4PiMinusStopMaterial 
{  

private:

  // Hide assignment operator as private 
  G4PiMinusStopMaterial& operator=(const G4PiMinusStopMaterial &right);

  // Copy constructor
  G4PiMinusStopMaterial(const G4PiMinusStopMaterial& );

public:

  // Constructor
  G4PiMinusStopMaterial();

  // Destructor
  virtual ~G4PiMinusStopMaterial();

  // Definitions of absorption products
  virtual RWTPtrOrderedVector<G4ParticleDefinition>* DefinitionVector();

  // 4-vectors of absorption products
  virtual RWTPtrOrderedVector<G4LorentzVector>* P4Vector(const G4double binding, 
							 const G4double mass);

  // Number of final nucleons, out of generated absorption products
  virtual G4double FinalNucleons()=0;

protected:
   
  RWTPtrOrderedVector<G4ParticleDefinition>* _definitions;
  RWTPtrOrderedVector<G4LorentzVector>* _momenta; 
  G4DistributionGenerator* _distributionE;
  G4DistributionGenerator* _distributionAngle;
  G4double _R; 

  G4double GenerateAngle(G4double range);
  G4LorentzVector MakeP4(G4double p, G4double theta, 
			 G4double phi, G4double e);
  G4double RecoilEnergy(const G4double mass);

private:

};
 

#endif
