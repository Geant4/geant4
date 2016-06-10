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
//
// $Id: G4ErrorSurfaceTrajParam.hh 69014 2013-04-15 09:42:51Z gcosmo $
//
//
// Class Description:
//
// Holds the 5 independent variables of the trajectory for a
// G4ErrorSurfaceTrajState object. It is not used for anything but for
// printing, but anyhow it is updated everytime the position and momentum
// are updated 

// History:
// - Created:  P. Arce
// --------------------------------------------------------------------

#ifndef G4ErrorSurfaceTrajParam_hh
#define G4ErrorSurfaceTrajParam_hh

#include "G4Point3D.hh"
#include "G4Vector3D.hh"
#include "G4Plane3D.hh"
#include "G4ThreeVector.hh"

#include "globals.hh"
#include "G4Track.hh"

class G4ErrorSurfaceTrajParam
{
 public:  // with description

  G4ErrorSurfaceTrajParam()
   : fInvP(0.), fPV(0.), fPW(0.), fV(0.), fW(0.) {}
  G4ErrorSurfaceTrajParam( const G4Point3D& pos, const G4Vector3D& mom,
                           const G4Vector3D& vecV, const G4Vector3D& vecW );
  G4ErrorSurfaceTrajParam( const G4Point3D& pos, const G4Vector3D& mom,
                           const G4Plane3D& plane );
  virtual ~G4ErrorSurfaceTrajParam(){}

  friend
    std::ostream& operator<<(std::ostream&, const G4ErrorSurfaceTrajParam& ts);
  
  // Get and Set methods 

  void SetParameters( const G4Point3D& pos, const G4Vector3D& mom,
                      const G4Vector3D& vecV, const G4Vector3D& vecW );
  void SetParameters( const G4Point3D& pos, const G4Vector3D& mom,
                      const G4Plane3D& plane );

  G4Vector3D GetDirection() const { return fDir; }
  G4Vector3D GetPlaneNormal() const { return fVectorV.cross(fVectorW); }
  G4Vector3D GetVectorV() const { return fVectorV; }
  G4Vector3D GetVectorW() const { return fVectorW; }
  G4double GetPV() const{ return fPV; }
  G4double GetPW() const{ return fPW; }
  G4double GetV() const{ return fV; }
  G4double GetW() const{ return fW; }
  G4double GetInvP() const{ return fInvP; }

 private:

  G4ThreeVector fDir;
  G4Vector3D fVectorV; //one of the vectors defining the plane
  G4Vector3D fVectorW; //one of the vectors defining the plane
  G4double fInvP; // inverse of momentum
  G4double fPV; // projection of momentum in one direction
  G4double fPW; // projection of momentum in one direction
  G4double fV; // projection of position in one direction
  G4double fW; // projection of position in one direction

};

#endif
