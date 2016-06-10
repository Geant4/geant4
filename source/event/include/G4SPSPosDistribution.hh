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
// MODULE:        G4SPSPosDistribution.hh
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
// 06/06/2014  A Dotti
//    For thread safety: this is a shared object,
//    mutex has been added to control access to shared resources (data members).
//    in Getters and Setters, mutex is NOT used in GenerateOne because it is
//    assumed that properties are not changed during event loop.
//
// Version 1.0, 05/02/2004, Fan Lei, Created.
//    Based on the G4GeneralParticleSource class in Geant4 v6.0
//
///////////////////////////////////////////////////////////////////////////////
//
//
// Class Description:
//
// To generate the position of a primary vertex according to the defined distribution 
//
///////////////////////////////////////////////////////////////////////////////
//
// MEMBER FUNCTIONS
// ----------------
//
// G4SPSPosDistribution ()
//    Constructor: Initializes variables and instantiates the Navigator class
//
// ~G4SPSPosDistribution ()
//    Destructor: 
//
// void SetPosDisType(G4String)
//    Allows user to choose Point, Plane, Surface or Volume source
//    position distributions.
//
// void SetPosDisShape(G4String)
//    Allows the user to choose the particular shape they wish for the
//    position distribution. Choices are Square, Circle, Ellipse, Rectangle,
//    Sphere, Ellipsoid, Cylinder, Parallelepiped.
//
// void SetCentreCoords(G4ThreeVector)
//    Sets the co-ordinates of the centre of the position distribution.
//
// void SetPosRot1(G4ThreeVector)
//    Used to specify the co-ordinate system for the position distribution
//    along with SetPosRot2. SetPosRot1 sets the vector x' and need not be
//    a unit vector.
//
// void SetPosRot2(G4ThreeVector)
//    Used in connection with SetPosRot1. This sets a vector in the plane
//    x'y'. By a series of cross products x', y', z' are generated. Again
//    need not be a unit vector.
// 
// void SetHalfX(G4double)
//    Sets the half length in x.
//
// void SetHalfY(G4double)
//    Sets the half length in y.
//
// void SetHalfZ(G4double)
//    Sets the half length in z.
//
// void SetRadius(G4double)
//    Sets the radius where appropriate for source distribution shapes.
//
// void SetRadius0(G4double)
//    Sets the inner radius where appropriate for source distribution shapes.
//
//  void SetBeamSigmaInR(G4double);
//    Sets the sigma for 1D beam
//
//  void SetBeamSigmaInX(G4double);
//    Sets the first sigma for 2D beam
// 
//  void SetBeamSigmaInY(G4double);
//    Sets the second sigma for 2D beam
//
// void SetParAlpha(G4double)
//    Sets the angle Alpha in the Parallelepiped shapes.
//
// void SetParTheta(G4double)
//    Sets the angle Theta in the Parallelepiped shapes.
//
// void SetParPhi(G4double)
//    Sets the angle Phi in the Parallelepiped shapes.
//
// void ConfineSourceToVolume(G4String)
//    Used to confine the start positions to a particular volume.
//
//  void SetBiasRndm (G4SPSRandomGenerator* a) { posRndm = a ; };
//    Sets the biased random number generator
//
//  G4ThreeVector GenerateOne();
//    Generate one random position
//
// void SetVerbosity(G4int)
//    Sets the verbosity level.
//
///////////////////////////////////////////////////////////////////////////////
//
#ifndef G4SPSPosDistribution_h
#define G4SPSPosDistribution_h 1

#include "G4Navigator.hh"
#include "G4SPSRandomGenerator.hh"
#include "G4Threading.hh"
#include "G4Cache.hh"

class G4SPSPosDistribution
{
public:
  G4SPSPosDistribution (); 
  ~G4SPSPosDistribution ();

  // methods to create source position dist.
  void SetPosDisType(G4String); // Point, Plane, Surface, Volume
  G4String GetPosDisType();
  void SetPosDisShape(G4String);
  G4String GetPosDisShape();
  // SetPosDisShape - Square, Circle, Annulus, Ellipse, Rectangle, Sphere,
  // Ellipsoid, Cylinder, Right (parallelepiped).
  void SetCentreCoords(G4ThreeVector);
  G4ThreeVector GetCentreCoords();
  void SetPosRot1(G4ThreeVector); 
  void SetPosRot2(G4ThreeVector); 
  void SetHalfX(G4double);
  G4double GetHalfX();
  void SetHalfY(G4double);
  G4double GetHalfY();
  void SetHalfZ(G4double);
  G4double GetHalfZ();
  void SetRadius(G4double);
  G4double GetRadius();
  void SetRadius0(G4double);
  void SetBeamSigmaInR(G4double);
  void SetBeamSigmaInX(G4double);
  void SetBeamSigmaInY(G4double);
  void SetParAlpha(G4double);
  void SetParTheta(G4double);
  void SetParPhi(G4double);
  void ConfineSourceToVolume(G4String);
  //
  void SetBiasRndm (G4SPSRandomGenerator* a);
  // Set the verbosity level.
  void SetVerbosity(G4int a);
  //
  G4ThreeVector GenerateOne();
    
    G4ThreeVector GetSideRefVec1() const;
    G4ThreeVector GetSideRefVec2() const;
    G4ThreeVector GetSideRefVec3() const;
    G4String GetSourcePosType();
    G4ThreeVector GetParticlePos() const;

private:

  void GenerateRotationMatrices();
  // the following routines generate the source position
  void GeneratePointSource(G4ThreeVector& outoutPos);
  void GeneratePointsInBeam(G4ThreeVector& outoutPos);
  void GeneratePointsInPlane(G4ThreeVector& outoutPos);
  void GeneratePointsOnSurface(G4ThreeVector& outputPos);
  void GeneratePointsInVolume(G4ThreeVector& outputPos);

  G4bool IsSourceConfined(G4ThreeVector& outputPos);

    void InitThreadLocalCache() const;
    void HandleThreadLocalCache( const G4ThreeVector& v1 , const G4ThreeVector& v2, const G4ThreeVector& v3 );
private:

  // Position distribution Variables
  G4String SourcePosType; //Point,Plane,Surface,Volume
  G4String Shape; //Circle,Square,Rectangle etc..
  G4double halfx, halfy, halfz; //half lengths
  G4double Radius; //Radius for circles or spheres
  G4double Radius0; // The inner radius of an annulus
  G4double SR,SX,SY; // Standard deviation in raduial, x, y for beam type source
  G4ThreeVector CentreCoords; // Coords of centre of input shape
  G4ThreeVector Rotx, Roty, Rotz; // Unit vectors defining rotation matrix
  G4double ParAlpha, ParTheta, ParPhi; //Angle for Right Parallellepipeds
  G4bool Confine; //If true confines source distribution to VolName
  G4String VolName;
  //Thread local caches
    G4Cache<G4ThreeVector*> CSideRefVec1, CSideRefVec2, CSideRefVec3;
    G4Cache<G4ThreeVector> CParticlePos;

  G4ThreeVector particle_position; // the final particle position
  //
  //G4Navigator *gNavigator;
  //
  G4SPSRandomGenerator* posRndm; // biased random generator
  // Verbosity
  G4int verbosityLevel;
    
  G4Mutex mutex; // protect access to shared resources

};

#endif




