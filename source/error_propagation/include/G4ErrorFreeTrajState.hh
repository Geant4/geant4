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
// $Id: G4ErrorFreeTrajState.hh 66892 2013-01-17 10:57:59Z gunter $
//
// Class Description:
//
// Represents a free G4ErrorTrajState
// It can be represented by the 5 variables
//     1/p, lambda, phi, y_perp, z_perp
// where lambda and phi are the dip and azimuthal angles related
// to the momentum components in the following way:
//            p_x = p cos(lambda) cos(phi)  ! lambda = 90 - theta
//            p_y = p cos(lambda) sin(phi)
//            p_z = p sin(lambda)
// y_perp and z_perp are the coordinates of the trajectory in a
// local orthonormal reference frame with the x_perp axis along the
// particle direction, the y_perp being parallel to the x-y plane.
//
// This class also takes care of propagating the error associated to 
// the trajectory 

// History:
// - Created:   P. Arce
// --------------------------------------------------------------------

#ifndef G4ErrorFreeTrajState_hh
#define G4ErrorFreeTrajState_hh

#include "globals.hh"

#include "G4ErrorMatrix.hh"

#include "G4ErrorTrajState.hh"
#include "G4ErrorFreeTrajParam.hh"

#include "G4Point3D.hh"
#include "G4Vector3D.hh"

class G4ErrorSurfaceTrajState;

class G4ErrorFreeTrajState : public G4ErrorTrajState
{
 public:  // with description

  G4ErrorFreeTrajState() : theFirstStep(true) {}
  G4ErrorFreeTrajState( const G4String& partName,
                        const G4Point3D& pos,
                        const G4Vector3D& mom,
                        const G4ErrorTrajErr& errmat = G4ErrorTrajErr(5,0) );
    // Constructor by providing particle, position and momentum

  G4ErrorFreeTrajState( const G4ErrorSurfaceTrajState& tpOS );
    // Constructor by providing G4ErrorSurfaceTrajState

  ~G4ErrorFreeTrajState(){}

  virtual G4int Update( const G4Track* aTrack );
    // update parameters from G4Track

  virtual G4int PropagateError( const G4Track* aTrack );
    // propagate the error along the step

  virtual void Dump( std::ostream& out = G4cout ) const;
    // dump TrajState parameters

  friend
    std::ostream& operator<<(std::ostream&, const G4ErrorFreeTrajState& ts);

  // Set and Get methods 

  virtual void SetPosition( const G4Point3D pos )
    { SetParameters( pos, fMomentum ); }

  virtual void SetMomentum( const G4Vector3D& mom )
    { SetParameters( fPosition, mom ); }

  void SetParameters( const G4Point3D& pos, const G4Vector3D& mom )
    {
      fPosition = pos;
      fMomentum = mom;
      fTrajParam.SetParameters( pos, mom );
    }

  G4ErrorFreeTrajParam GetParameters() const
    { return fTrajParam; }

  G4ErrorMatrix GetTransfMat() const
    { return theTransfMat; }

 private:  

  void Init();
    // define TrajState type and build charge

  G4int PropagateErrorMSC( const G4Track* aTrack );
    // add the error associated to multiple scattering

  void CalculateEffectiveZandA( const G4Material* mate, double& effZ, double& effA );
    // calculate effective Z and A (needed by PropagateErrorMSC)
  
  G4int PropagateErrorIoni( const G4Track* aTrack );
    // add the error associated to ionization energy loss


 private:

  G4ErrorFreeTrajParam fTrajParam;

  G4ErrorMatrix theTransfMat;

  G4bool theFirstStep; // to count if transf mat is updated or initialized
};

#endif
