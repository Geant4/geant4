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
// $Id: G4ErrorSurfaceTrajState.hh 69014 2013-04-15 09:42:51Z gcosmo $
//
//
// Class description:
//
// Represents a trajectory state on a surface.
// It can be represented by the 5 variables:
//      1/p, v', w', v, w
// where v'=dv/du and w'=dw/du in an orthonormal coordinate system
// with axis u, v and w.

// History:
//
// - Created:   P. Arce 
// --------------------------------------------------------------------

#ifndef G4ErrorSurfaceTrajState_hh
#define G4ErrorSurfaceTrajState_hh

#include "globals.hh"

#include "G4ErrorTrajState.hh"
#include "G4ErrorSurfaceTrajParam.hh"
#include "G4ErrorFreeTrajState.hh"

#include "G4Point3D.hh"
#include "G4Vector3D.hh"
#include "G4Plane3D.hh"

class G4ErrorSurfaceTrajState : public G4ErrorTrajState
{

 public:  // with description

  G4ErrorSurfaceTrajState( const G4String& partType, const G4Point3D& pos,
                           const G4Vector3D& mom, const G4Plane3D& plane,
                           const G4ErrorTrajErr& errmat = G4ErrorTrajErr(5,0) );
    // Constructor by providing particle, position, momentum and
    // G4Plane3D surface

  G4ErrorSurfaceTrajState( const G4String& partType, const G4Point3D& pos,
                           const G4Vector3D& mom, const G4Vector3D& vecV,
                           const G4Vector3D& vecW,
                           const G4ErrorTrajErr& errmat = G4ErrorTrajErr(5,0) );
    // Constructor by providing particle, position, momentum and
    // two vectors on surface

  G4ErrorSurfaceTrajState( G4ErrorFreeTrajState& tpSC, const G4Plane3D& plane );
    // Constructor by providing G4ErrorFreeTrajState and G4Plane3D surface

  G4ErrorSurfaceTrajState( G4ErrorFreeTrajState& tpSC, const G4Vector3D& vecV,
                           const G4Vector3D& vecW  , G4ErrorMatrix &transfM);
    // Constructor by providing G4ErrorFreeTrajState and two vectors on surface

  ~G4ErrorSurfaceTrajState(){}
  G4ErrorMatrix   BuildErrorMatrix( G4ErrorFreeTrajState& tpSC, const G4Vector3D& vecV,
				    const G4Vector3D& vecW );
    // Build the error matrix from a free state plus the vectors of the surface

  virtual void Dump( std::ostream& out = G4cout ) const;

  friend
    std::ostream& operator<<(std::ostream&, const G4ErrorSurfaceTrajState& ts);

  // Set and Get methods 

  G4ErrorSurfaceTrajParam GetParameters() const
    { return fTrajParam; }
  void SetParameters( const G4Point3D& pos, const G4Vector3D& mom,
                      const G4Vector3D& vecV, const G4Vector3D& vecW )
    {
      fPosition = pos;
      fMomentum = mom;
      fTrajParam.SetParameters( pos, mom, vecV, vecW );
    }

  void SetParameters( const G4Point3D& pos, const G4Vector3D& mom,
                      const G4Plane3D& plane )
    {
      fPosition = pos;
      fMomentum = mom;
      fTrajParam.SetParameters( pos, mom, plane );
    }

  G4Vector3D GetVectorV() const
    { return fTrajParam.GetVectorV(); }

  G4Vector3D GetVectorW() const
    { return fTrajParam.GetVectorW(); }

  virtual void SetPosition( const G4Point3D pos )
    { SetParameters( pos, fMomentum, GetVectorV(), GetVectorW() ); }

  virtual void SetMomentum( const G4Vector3D& mom )
    { SetParameters( fPosition, mom, GetVectorV(), GetVectorW() ); }
    
 private:

  void Init();
    // Define Trajectory State type and build charge

 private:

  G4ErrorSurfaceTrajParam fTrajParam;

};

#endif
