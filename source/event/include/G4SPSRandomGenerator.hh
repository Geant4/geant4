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
///////////////////////////////////////////////////////////////////////////////
//
// MODULE:        G4SPSRandomGenerator.hh
//
// Version:      1.0
// Date:         5/02/04
// Author:       Fan Lei 
// Organisation: QinetiQ ltd.
// Customer:     ESA/ESTEC
//
///////////////////////////////////////////////////////////////////////////////
//
// CHANGE HISTORY
// --------------
//
//
// Version 1.0, 05/02/2004, Fan Lei, Created.
//    Based on the G4GeneralParticleSource class in Geant4 v6.0
//
///////////////////////////////////////////////////////////////////////////////
//
// Class Description:
//
// Special random number generator used by G4GeneralParticleSource to allow 
// biasing applied at the lowest level for all distributions.
//
///////////////////////////////////////////////////////////////////////////////
//
// MEMBER FUNCTIONS
// ----------------
//
// G4SPSRandomGenerator ()
//    Constructor: Initializes variables
//
// ~G4SPSRandomGenerator ()
//    Destructor: 
//
// void SetXBias(G4ThreeVector)
//    Allows the user to re-distribute the random
//    numbers used to generate x co-ordinates.
//
// void SetYBias(G4ThreeVector)
//    Allows the user to re-distribute the random
//    numbers used to generate y co-ordinates.
//
// void SetZBias(G4ThreeVector)
//    Allows the user to re-distribute the random
//    numbers used to generate z co-ordinates.
//
// void SetThetaBias(G4ThreeVector)
//    Allows the user to re-distribute the random
//    numbers used to generate values of theta.
//
// void SetPhiBias(G4ThreeVector)
//    Allows the user to re-distribute the random
//    numbers used to generate values of phi.
//
// void SetEnergyBias(G4ThreeVector)
//    Allows the user to re-distribute the random
//    numbers used to generate the energies.
//
// G4double GenRandX()
//    Generates the random number for x, with or without biasing.
//
// G4double GenRandY()
//    Generates the random number for y, with or without biasing.
//
// G4double GenRandZ()
//    Generates the random number for z, with or without biasing.
//
// G4double GenRandTheta()
//    Generates the random number for theta, with or without biasing.
//
// G4double GenRandPhi()
//    Generates the random number for phi, with or without biasing.
//
// G4double GenRandEnergy()
//    Generates the random number for energy, with or without biasing.
//
//  inline G4double GetBiasWeight()
//    Returns the weight change after biasing
// 
//  void ReSetHist(G4String);
//    Re-sets the histogram for user defined distribution
//
// void SetVerbosity(G4int)
//    Sets the verbosity level.
//
///////////////////////////////////////////////////////////////////////////////
//
#ifndef G4SPSRandomGenerator_h
#define G4SPSRandomGenerator_h 1

#include "G4PhysicsOrderedFreeVector.hh"
#include "G4DataInterpolation.hh"

class G4SPSRandomGenerator
{
public:
  G4SPSRandomGenerator (); 
  ~G4SPSRandomGenerator ();

  //  static G4SPSRandomGenerator* getInstance ();

  // Biasing Methods 
  void SetXBias(G4ThreeVector);
  void SetYBias(G4ThreeVector);
  void SetZBias(G4ThreeVector);
  void SetThetaBias(G4ThreeVector);
  void SetPhiBias(G4ThreeVector);
  void SetEnergyBias(G4ThreeVector);
  G4double GenRandX();
  G4double GenRandY();
  G4double GenRandZ();
  G4double GenRandTheta();
  G4double GenRandPhi();
  G4double GenRandEnergy();

  inline G4double GetBiasWeight() 
  { return bweights[0]*bweights[1]*bweights[2]*bweights[3]*bweights[4]*bweights[5];};
 
 // method to re-set the histograms
  void ReSetHist(G4String);

  // Set the verbosity level.
   void SetVerbosity(G4int a) {verbosityLevel = a; } ;

private:

  //  static G4SPSRandomGenerator  *instance;

  G4bool XBias, IPDFXBias;
  G4PhysicsOrderedFreeVector XBiasH;
  G4PhysicsOrderedFreeVector IPDFXBiasH;
  G4bool YBias, IPDFYBias;
  G4PhysicsOrderedFreeVector YBiasH;
  G4PhysicsOrderedFreeVector IPDFYBiasH;
  G4bool ZBias, IPDFZBias;
  G4PhysicsOrderedFreeVector ZBiasH;
  G4PhysicsOrderedFreeVector IPDFZBiasH;
  G4bool ThetaBias, IPDFThetaBias;
  G4PhysicsOrderedFreeVector ThetaBiasH;
  G4PhysicsOrderedFreeVector IPDFThetaBiasH;
  G4bool PhiBias, IPDFPhiBias;
  G4PhysicsOrderedFreeVector PhiBiasH;
  G4PhysicsOrderedFreeVector IPDFPhiBiasH;
  G4bool EnergyBias, IPDFEnergyBias;
  G4PhysicsOrderedFreeVector EnergyBiasH;
  G4PhysicsOrderedFreeVector IPDFEnergyBiasH;

  G4double bweights[6]; //record x,y,z,theta,phi,energy weights

  // Verbosity
  G4int verbosityLevel;

  G4PhysicsOrderedFreeVector ZeroPhysVector ; // for re-set only 
  
};


#endif




