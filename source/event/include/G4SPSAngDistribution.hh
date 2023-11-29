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
// G4SPSAngDistribution
//
// Class description:
//
// To generate the direction of a primary vertex according to the
// defined distribution. This is a shared class between threads.
// Only one thread should use the set-methods here.
// Note that this is exactly what is achieved using UI commands.
// If you use the set methods to set defaults in your
// application take care that only one thread is executing them.
// In addition take care of calling these methods before the run is started
// Do not use these setters during the event loop.

// Author: Fan Lei, QinetiQ ltd.
// Customer: ESA/ESTEC
// History:
// - 05/02/2004, Fan Lei - Created.
//    Based on the G4GeneralParticleSource class.
// - 26/10/2004, Fan Lei
//    Added a "focused" option to allow all primary particles pointing to 
//    a user specified focusing point. 
// - 06/06/2014, Andrea Dotti, SLAC
//    For thread safety: this is a shared object,
//    Added mutex to control access to shared resources (data members).
//    in Getters and Setters, mutex is NOT used in GenerateOne() because it is
//    assumed that properties are not changed during event loop.
// --------------------------------------------------------------------
#ifndef G4SPSAngDistribution_hh
#define G4SPSAngDistribution_hh 1

#include "G4PhysicsFreeVector.hh"
#include "G4DataInterpolation.hh"
#include "G4ParticleMomentum.hh"

#include "G4SPSPosDistribution.hh"
#include "G4SPSRandomGenerator.hh"

#include "G4Threading.hh"
#include "G4AutoLock.hh"

class G4SPSAngDistribution 
{

  public:

    G4SPSAngDistribution (); 
      // Constructor: Initializes variables
   ~G4SPSAngDistribution ();
      // Destructor

    // Angular Distribution Methods

    void SetAngDistType(const G4String&);
      // Used to set the type of angular distribution wanted. Arguments
      // are iso, cos, beam and user for isotropic, cosine-law, beam and
      // user-defined respectively.

    void DefineAngRefAxes(const G4String&, const G4ThreeVector&);
      // Used in a similar way as SetPosRot() to define vectors, one x'
      // and one in the plane x'y', to create a rotated set of axes for
      // the angular distribution.

    void SetMinTheta(G4double);
      // Sets the minimum value for the angle theta.

    void SetMinPhi(G4double);
      // Sets the minimum value for phi.  

    void SetMaxTheta(G4double);
      // Sets the maximum value for theta.

    void SetMaxPhi(G4double);
      // Sets the maximum value for phi.

    void SetBeamSigmaInAngR(G4double);
      // Sets the sigma for 1D beam.

    void SetBeamSigmaInAngX(G4double);
      // Sets the first sigma for 2D beam.

    void SetBeamSigmaInAngY(G4double);
      // Sets the second sigma for 2D beam.

    void UserDefAngTheta(const G4ThreeVector&);
      // This method allows the user to define a histogram in Theta.

    void UserDefAngPhi(const G4ThreeVector&);
      // This method allows the user to define a histogram in phi.

    void SetFocusPoint(const G4ThreeVector&);
    void SetParticleMomentumDirection(const G4ParticleMomentum& aMomDirection);
    void SetUseUserAngAxis(G4bool);
    void SetUserWRTSurface(G4bool);
      // Allows user to have user-defined spectra either with respect to the
      // coordinate system (default) or with respect to the surface normal.

    void SetPosDistribution(G4SPSPosDistribution* a);
      // Sets the required position generator, required for determining
      // the cosine-law distribution.

    void SetBiasRndm(G4SPSRandomGenerator* a);
      // Sets the biased random number generator.

    G4ParticleMomentum GenerateOne();
      // Generates one random direction.

    void ReSetHist(const G4String&);
      // Re-sets the histogram for user defined distribution..

    void SetVerbosity(G4int a);
      // Sets the verbosity level.

    // Some accessors

    G4String GetDistType(); 
    G4double GetMinTheta();
    G4double GetMaxTheta();
    G4double GetMinPhi();
    G4double GetMaxPhi();
    G4ThreeVector GetDirection();

  private:

  // These methods generate the momentum vectors for the particles.

    void GenerateFocusedFlux(G4ParticleMomentum& outputMom);
      // This method generates momentum vectors for particles pointing to
      // an user specified focusing point.

    void GenerateIsotropicFlux(G4ParticleMomentum& outputMom);
      // This method generates momentum vectors for particles according
      // to an isotropic distribution.

    void GenerateCosineLawFlux(G4ParticleMomentum& outputMom);
      // This method generates momentum vectors for particles according
      // to a cosine-law distribution.

    void GenerateBeamFlux(G4ParticleMomentum& outputMom);
    void GeneratePlanarFlux(G4ParticleMomentum& outputMom);

    void GenerateUserDefFlux(G4ParticleMomentum& outputMom);
      // Controls generation of momentum vectors according to user-defined
      // distributions.

    G4double GenerateUserDefTheta();
      // Generates the theta angle according to a user-defined distribution.

    G4double GenerateUserDefPhi();
      // Generates phi according to a user-defined distribution.

  private:

    // Angular distribution variables.

    G4String AngDistType; // String to hold Ang dist type iso, cos, user
    G4ThreeVector AngRef1, AngRef2, AngRef3; // Reference axes for ang dist
    G4double MinTheta, MaxTheta, MinPhi, MaxPhi; // min/max theta/phi
    G4double DR,DX,DY ; // Standard deviations for beam divergence 
    G4double Theta{0.}, Phi{0.}; // Store these for use with DEBUG
    G4ThreeVector FocusPoint ; // the focusing point in mother coordinates
    G4bool IPDFThetaExist, IPDFPhiExist; // tell whether IPDF histos exist
    G4PhysicsFreeVector UDefThetaH; // Theta histo data
    G4PhysicsFreeVector IPDFThetaH; //Cumulative Theta histogram.
    G4PhysicsFreeVector UDefPhiH; // Phi histo bins
    G4PhysicsFreeVector IPDFPhiH; // Cumulative phi histogram.
    G4String UserDistType; //String to hold user distributions
    G4bool UserWRTSurface; // G4bool to tell whether user wants distribution wrt
                           // surface normals or co-ordinate system
    G4bool UserAngRef; // Set to true when user defines a new coordinates

    G4ParticleMomentum particle_momentum_direction;

    G4SPSPosDistribution* posDist = nullptr; // need for the cosine-law distrib.
    G4SPSRandomGenerator* angRndm = nullptr; // biased random generator

    G4int verbosityLevel; // Verbosity

    G4PhysicsFreeVector ZeroPhysVector; // for re-set only

    G4Mutex mutex; // Protect access to shared resources
};

#endif
