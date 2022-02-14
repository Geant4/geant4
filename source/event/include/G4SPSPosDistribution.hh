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
// G4SPSPosDistribution
//
// Class Description:
//
// To generate the position of a primary vertex according to
// the defined distribution. This is a shared class between threads.
// Only one thread should use the set-methods here.
// Note that this is exactly what is achieved using UI commands.
// If you use the set methods to set defaults in your application take care
// that only one thread is executing them.
// In addition take care of calling these methods before the run is started
// Do not use these setters during the event loop

// Author: Fan Lei, QinetiQ ltd.
// Customer: ESA/ESTEC
// History:
// - 05/02/2004, Fan Lei, Created.
//     Based on the G4GeneralParticleSource class.
// - 06/06/2014, Andrea Dotti
//     Added mutex to control access to shared resources (data members).
//     In Getters and Setters, mutex is NOT used in GenerateOne because
//     it is assumed that properties are not changed during event loop.
// - 13/02/2017, Maxime Chauvin
//     Added surface and volume shape "EllipticCylinder"
// --------------------------------------------------------------------
#ifndef G4SPSPosDistribution_hh
#define G4SPSPosDistribution_hh 1

#include "G4Navigator.hh"
#include "G4SPSRandomGenerator.hh"
#include "G4Threading.hh"
#include "G4Cache.hh"

class G4SPSPosDistribution
{
  public:

    G4SPSPosDistribution();
      // Constructor: initializes data and instantiates the Navigator class

   ~G4SPSPosDistribution();
      // Destructor

    // Methods to create source position dist

    void SetPosDisType(const G4String&);
      // Allows user to choose Point, Plane, Surface or Volume source
      // position distributions

    void SetPosDisShape(const G4String&);
      // Allows the user to choose the particular shape they wish for the
      // position distribution. Choices are: Square, Circle, Ellipse,
      // Rectangle, Sphere, Ellipsoid, Cylinder, Parallelepiped

    void SetCentreCoords(const G4ThreeVector&);
      // Sets the coordinates of the centre of the position distribution

    void SetPosRot1(const G4ThreeVector&); 
      // Used to specify the coordinate system for the position distribution
      // along with SetPosRot2. Sets the vector x' and need not be a unit vector

    void SetPosRot2(const G4ThreeVector&); 
      // Used in connection with SetPosRot1. This sets a vector in the plane
      // x'y'. By a series of cross products x', y', z' are generated. Again
      // need not be a unit vector

    void SetHalfX(G4double);
      // Sets the half length in x

    void SetHalfY(G4double);
      // Sets the half length in y

    void SetHalfZ(G4double);
      // Sets the half length in z

    void SetRadius(G4double);
      // Sets the radius where appropriate for source distribution shapes

    void SetRadius0(G4double);
      // Sets the inner radius where appropriate for source distribution shapes

    void SetBeamSigmaInR(G4double);
      // Sets the sigma for 1D beam

    void SetBeamSigmaInX(G4double);
      // Sets the first sigma for 2D beam

    void SetBeamSigmaInY(G4double);
      // Sets the second sigma for 2D beam

    void SetParAlpha(G4double);
      // Sets the angle Alpha in the Parallelepiped shapes

    void SetParTheta(G4double);
      // Sets the angle Theta in the Parallelepiped shapes

    void SetParPhi(G4double);
      // Sets the angle Phi in the Parallelepiped shapes

    void ConfineSourceToVolume(const G4String&);
      // Used to confine the start positions to a particular volume

    void SetBiasRndm (G4SPSRandomGenerator* a);
      // Sets the biased random number generator

    void SetVerbosity(G4int a);
      // Sets the verbosity level

    G4ThreeVector GenerateOne();
      // Generate one random position

    const G4String& GetPosDisType() const;
    const G4String& GetPosDisShape() const;
    const G4ThreeVector& GetCentreCoords() const;
    G4double GetHalfX() const;
    G4double GetHalfY() const;
    G4double GetHalfZ() const;
    G4double GetRadius() const;
    inline G4double GetRadius0() const { return Radius0; }
    inline G4double GetParAlpha() const { return ParAlpha; }
    inline G4double GetParTheta() const { return ParTheta; }
    inline G4double GetParPhi()   const { return ParPhi; }
    inline const G4ThreeVector& GetRotx() const { return Rotx; }
    inline const G4ThreeVector& GetRoty() const { return Roty; }
    inline const G4ThreeVector& GetRotz() const { return Rotz; }
    inline G4bool GetConfined() const { return Confine; }
    inline const G4String& GetConfineVolume() const { return VolName; }

    const G4ThreeVector& GetSideRefVec1() const;
    const G4ThreeVector& GetSideRefVec2() const;
    const G4ThreeVector& GetSideRefVec3() const;
    const G4String& GetSourcePosType() const;
    const G4ThreeVector& GetParticlePos() const;

  private:

    void GenerateRotationMatrices();

    // The following functions generate the source position
    //
    void GeneratePointSource(G4ThreeVector& outoutPos);
    void GeneratePointsInBeam(G4ThreeVector& outoutPos);
    void GeneratePointsInPlane(G4ThreeVector& outoutPos);
    void GeneratePointsOnSurface(G4ThreeVector& outputPos);
    void GeneratePointsInVolume(G4ThreeVector& outputPos);

    G4bool IsSourceConfined(G4ThreeVector& outputPos);

  private:

    // NOTE:
    // This is a shared resource, however setters that
    // changes the parameters via UI commands are by design
    // thread-safe because only one thread will call these methods
    // See G4GeneralParticleSourceMessenger constructor for an explanation
    //
    struct thread_data_t  // Caching of some data
    {
      G4ThreeVector CSideRefVec1;
      G4ThreeVector CSideRefVec2;
      G4ThreeVector CSideRefVec3;
      G4ThreeVector CParticlePos;
      thread_data_t();
    };

    G4String SourcePosType;
      // Point, Plane, Surface, Volume
    G4String Shape;
      // Circle, Square, Rectangle, etc...
    G4ThreeVector CentreCoords;
      // Coordinates of centre of input shape
    G4ThreeVector Rotx, Roty, Rotz;
      // Unit vectors defining rotation matrix
    G4double halfx, halfy, halfz;
      // Half lengths
    G4double Radius;
      // Radius for circles or spheres
    G4double Radius0;
      // The inner radius of an annulus
    G4double SR, SX, SY;
      // Standard deviation in radial, x, y for beam type source
    G4double ParAlpha,  ParTheta, ParPhi;
      // Angle for Right Parallellepipeds
    G4bool Confine = false;
      // If true confines source distribution to VolName
    G4String VolName;
      // Volume name
    G4int verbosityLevel;
      // Verbosity
    G4SPSRandomGenerator* PosRndm = nullptr;
     // Biased random generator

    G4Cache<thread_data_t> ThreadData;
    G4Mutex a_mutex;
};

#endif
