// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4PiMinusStopCo.cc,v 1.2 1999-11-11 15:37:46 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// -------------------------------------------------------------------
//      GEANT 4 class file --- Copyright CERN 1998
//      CERN Geneva Switzerland
//
//      For information related to this code contact:
//      CERN, IT Division, ASD group
//
//      File name:     G4PiMinusStopCo
//
//      Author:        Maria Grazia Pia (pia@genova.infn.it)
// 
//      Creation date: 8 May 1998
//
//      Modifications: 
// -------------------------------------------------------------------

#include "G4ios.hh"

#include "G4PiMinusStopCo.hh"

#include "g4rw/tpordvec.h"
#include "g4rw/tvordvec.h"
#include "g4rw/cstring.h"

#include "globals.hh"
#include "Randomize.hh"
#include "G4Proton.hh"
#include "G4Neutron.hh"
#include "G4ParticleTypes.hh"
#include "G4ReactionKinematics.hh"
#include "G4DynamicParticleVector.hh"
#include "G4LorentzVector.hh"
#include "G4NucleiPropertiesTable.hh"
#include "G4PiMinusStopMaterial.hh"
#include "G4DistributionGenerator.hh"

// np/pp production ratio
// Experimental values: 
// R(np/pp) (E. Gadioli et al., Phys Rev C36 (1987) 741
G4double G4PiMinusStopCo::npRatio = 2.5;


 
// Average numbers of final nucleons detected, for N-pair absorption
// (Hartmann et al., Nucl. Phys. A300 (1978) 345
G4double G4PiMinusStopCo::nFinalNucleons = 1.38;

// Kinetic energy (MeV) distributions measured for coincident nucleon 
// emission
// P. Heusi et al., Nucl. Phys. A407(1983) 429

G4int G4PiMinusStopCo::eKinEntries = 11;

G4double G4PiMinusStopCo::eKinData[11] = {0.085, 0.09, 0.09, 0.15, 
				      0.1, 0.09, 0.08, 
				      0.04,  0.03,  0.02, 0.01};

G4double G4PiMinusStopCo::eKin[12] = { 15., 17.5, 25.,  33.,  
				   42.,  52.,  62., 
				   75., 85., 95., 105. };


// Opening angle distributions measured for coincident nucleon emission
// (P.Heusi et al., Nucl. Phys. A407 (1983) 429

G4int G4PiMinusStopCo::angleEntries = 7;

G4double G4PiMinusStopCo::angleData[7] = 
{6., 8., 9., 10., 25., 40., 45. };

G4double G4PiMinusStopCo::angle[8] = { 0.5, 0.7, 0.87, 1.4, 2.1, 
				   2.44, 2.8, 3.1415927 };



// Constructor

G4PiMinusStopCo::G4PiMinusStopCo()
  
{
  // Cluster size: nucleon pair, alpha, triton etc.
  // First implementation: interaction with nucleon pair only
  _clusterSize = 2;

  // R ratio
  _R = 1. / (1. + npRatio);

  _definitions = new G4RWTPtrOrderedVector<G4ParticleDefinition>();
  _momenta = new G4RWTPtrOrderedVector<G4LorentzVector>();

  G4RWTValOrderedVector<double> eKinVector;
  G4RWTValOrderedVector<double> eKinDataVector;
  int i;
  for (i=0; i<eKinEntries; i++)
    {
      eKinVector.insert(eKin[i]);
      eKinDataVector.insert(eKinData[i]);
    }
  eKinVector.insert(eKin[eKinEntries]);
  _distributionE = new G4DistributionGenerator(eKinVector,eKinDataVector);

  G4RWTValOrderedVector<double> angleVector;
  G4RWTValOrderedVector<double> angleDataVector;
  for (i=0; i<angleEntries; i++)
    {
      angleVector.insert(angle[i]);
      angleDataVector.insert(angleData[i]);
    }
  angleVector.insert(angle[angleEntries]);
  _distributionAngle = new G4DistributionGenerator(angleVector,angleDataVector);
}


// Destructor

G4PiMinusStopCo::~G4PiMinusStopCo()
{}

G4double G4PiMinusStopCo::FinalNucleons()
{
  return nFinalNucleons;
}

