//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
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
// 26/10/2004 F Lei
//    Created separated the theta, phi generators for position distributions.
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
// void SetPosThetaBias(G4ThreeVector)
//    Allows the user to re-distribute the random
//    numbers used to generate values of theta for position distribution.
//
// void SetPosPhiBias(G4ThreeVector)
//    Allows the user to re-distribute the random
//    numbers used to generate values of phi for position distribution.
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
// G4double GenRandPosTheta()
//    Generates the random number for theta, with or without biasing for position distribution.
//
// G4double GenRandPosPhi()
//    Generates the random number for phi, with or without biasing for position distribution.
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

class G4SPSRandomGenerator {
public:
	G4SPSRandomGenerator();
	~G4SPSRandomGenerator();

	//  static G4SPSRandomGenerator* getInstance ();

	// Biasing Methods
	void SetXBias(G4ThreeVector);
	void SetYBias(G4ThreeVector);
	void SetZBias(G4ThreeVector);
	void SetThetaBias(G4ThreeVector);
	void SetPhiBias(G4ThreeVector);
	void SetEnergyBias(G4ThreeVector);
	void SetPosThetaBias(G4ThreeVector);
	void SetPosPhiBias(G4ThreeVector);
	G4double GenRandX();
	G4double GenRandY();
	G4double GenRandZ();
	G4double GenRandTheta();
	G4double GenRandPhi();
	G4double GenRandEnergy();
	G4double GenRandPosTheta();
	G4double GenRandPosPhi();

	inline void SetIntensityWeight(G4double weight) {
		bweights[8] = weight;
	}
	;

	inline G4double GetBiasWeight() {
		return bweights[0] * bweights[1] * bweights[2] * bweights[3]
				* bweights[4] * bweights[5] * bweights[6] * bweights[7]
				* bweights[8];
	}
	;

	// method to re-set the histograms
	void ReSetHist(G4String);

	// Set the verbosity level.
	void SetVerbosity(G4int a) {
		verbosityLevel = a;
	}
	;

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
	G4bool PosThetaBias, IPDFPosThetaBias;
	G4PhysicsOrderedFreeVector PosThetaBiasH;
	G4PhysicsOrderedFreeVector IPDFPosThetaBiasH;
	G4bool PosPhiBias, IPDFPosPhiBias;
	G4PhysicsOrderedFreeVector PosPhiBiasH;
	G4PhysicsOrderedFreeVector IPDFPosPhiBiasH;

	G4double alpha;   // for biasing energy

	G4double bweights[9]; //record x,y,z,theta,phi,energy,posThet,posPhi,intensity weights

	// Verbosity
	G4int verbosityLevel;

	G4PhysicsOrderedFreeVector ZeroPhysVector; // for re-set only

};

#endif

