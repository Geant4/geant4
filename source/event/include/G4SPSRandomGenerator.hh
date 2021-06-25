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
// G4SPSRandomGenerator
//
// Class Description:
//
// Special random number generator used by G4GeneralParticleSource to allow 
// biasing applied at the lowest level for all distributions.
// This is a shared class between threads.
// Only one thread should use the set-methods here.
// Note that this is exactly what is achieved using UI commands.
// If you use the set methods to set defaults in your
// application take care that only one thread is executing them.
// In addition take care of calling these methods before the run is started
// Do not use the setters during the event loop

// Author: Fan Lei, QinetiQ ltd.
// Customer: ESA/ESTEC
// History:
// - 05/02/2004, Fan Lei - Created.
//     Based on the G4GeneralParticleSource class
// - 06/06/2014, Andrea Dotti
//     Added a mutex to protect access to shared resources (data members).
//     Getters and Setters are mutex'd but not the GetRand* methods,
//     because it is assumed these are called only during the event loop
//     during which the status of this class is invariant
// --------------------------------------------------------------------
#ifndef G4SPSRandomGenerator_hh
#define G4SPSRandomGenerator_hh 1

#include "G4PhysicsFreeVector.hh"
#include "G4DataInterpolation.hh"
#include "G4ThreeVector.hh"
#include "G4Threading.hh"
#include "G4Cache.hh"

class G4SPSRandomGenerator
{
  public:

    G4SPSRandomGenerator();
      // Constructor: initializes variables

   ~G4SPSRandomGenerator();
      // Destructor

    // Biasing Methods

    void SetXBias(const G4ThreeVector&);
      // Allows the user to re-distribute the random
      // numbers used to generate x co-ordinates

    void SetYBias(const G4ThreeVector&);
      // Allows the user to re-distribute the random
      // numbers used to generate y co-ordinates

    void SetZBias(const G4ThreeVector&);
      // Allows the user to re-distribute the random
      // numbers used to generate z co-ordinates

    void SetThetaBias(const G4ThreeVector&);
      // Allows the user to re-distribute the random
      // numbers used to generate values of theta

    void SetPhiBias(const G4ThreeVector&);
      // Allows the user to re-distribute the random
      // numbers used to generate values of phi

    void SetEnergyBias(const G4ThreeVector&);
    // Allows the user to re-distribute the random
    // numbers used to generate the energies

    void SetPosThetaBias(const G4ThreeVector&);
      // Allows the user to re-distribute the random
      // numbers used to generate values of theta for position distribution

    void SetPosPhiBias(const G4ThreeVector&);
      // Allows the user to re-distribute the random
      // numbers used to generate values of phi for position distribution

    G4double GenRandX();
      // Generates the random number for x, with or without biasing

    G4double GenRandY();
      // Generates the random number for y, with or without biasing

    G4double GenRandZ();
      // Generates the random number for z, with or without biasing

    G4double GenRandTheta();
      // Generates the random number for theta, with or without biasing

    G4double GenRandPhi();
      // Generates the random number for phi, with or without biasing

    G4double GenRandEnergy();
      // Generates the random number for energy, with or without biasing

    G4double GenRandPosTheta();
      // Generates the random number for theta, with or without biasing
      // for position distribution

    G4double GenRandPosPhi();
      // Generates the random number for phi, with or without biasing
      // for position distribution

    void SetIntensityWeight(G4double weight);

    G4double GetBiasWeight() const ;
      // Returns the weight change after biasing

        // method to re-set the histograms
    void ReSetHist(const G4String&);
      // Resets the histogram for user defined distribution

    void SetVerbosity(G4int a);
      // Sets the verbosity level

  private:

    // Encapsulate in a struct to guarantee that correct
    // initial state is set via constructor
    //
    struct a_check
    {
      G4bool val;
      a_check() { val = false; }
    };
    
    // See .cc for an explanation of this in method GenRandX()
    //
    G4Cache<a_check> local_IPDFXBias;
    G4bool XBias, IPDFXBias;
    G4PhysicsFreeVector XBiasH;
    G4PhysicsFreeVector IPDFXBiasH;
    G4Cache<a_check> local_IPDFYBias;
    G4bool YBias, IPDFYBias;
    G4PhysicsFreeVector YBiasH;
    G4PhysicsFreeVector IPDFYBiasH;
    G4Cache<a_check> local_IPDFZBias;
    G4bool ZBias, IPDFZBias;
    G4PhysicsFreeVector ZBiasH;
    G4PhysicsFreeVector IPDFZBiasH;
    G4Cache<a_check> local_IPDFThetaBias;
    G4bool ThetaBias, IPDFThetaBias;
    G4PhysicsFreeVector ThetaBiasH;
    G4PhysicsFreeVector IPDFThetaBiasH;
    G4Cache<a_check> local_IPDFPhiBias;
    G4bool PhiBias, IPDFPhiBias;
    G4PhysicsFreeVector PhiBiasH;
    G4PhysicsFreeVector IPDFPhiBiasH;
    G4Cache<a_check> local_IPDFEnergyBias;
    G4bool EnergyBias, IPDFEnergyBias;
    G4PhysicsFreeVector EnergyBiasH;
    G4PhysicsFreeVector IPDFEnergyBiasH;
    G4Cache<a_check> local_IPDFPosThetaBias;
    G4bool PosThetaBias, IPDFPosThetaBias;
    G4PhysicsFreeVector PosThetaBiasH;
    G4PhysicsFreeVector IPDFPosThetaBiasH;
    G4Cache<a_check> local_IPDFPosPhiBias;
    G4bool PosPhiBias, IPDFPosPhiBias;
    G4PhysicsFreeVector PosPhiBiasH;
    G4PhysicsFreeVector IPDFPosPhiBiasH;

    struct bweights_t
    {
      G4double w[9];
      bweights_t();
      G4double& operator[] (const int i);
    };
    G4Cache<bweights_t> bweights;
      // record x,y,z,theta,phi,energy,posThet,posPhi,intensity weights

    G4int verbosityLevel;
      // Verbosity
 
    G4Mutex mutex;
      // Protect shared resources
};

#endif
