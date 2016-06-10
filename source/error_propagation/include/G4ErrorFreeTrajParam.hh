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
// $Id: G4ErrorFreeTrajParam.hh 66892 2013-01-17 10:57:59Z gunter $
//
// Class Description:
//
// Holds the 5 independent variables of the trajectory for a
// G4ErrorFreeTrajState object. It is not used for anything but for
// printing, but anyhow it is updated everytime the position and
// momentum are updated.

// History:
// - Created: Pedro Arce, September 2004
// --------------------------------------------------------------------

#ifndef G4ErrorFreeTrajParam_hh
#define G4ErrorFreeTrajParam_hh

#include "G4Point3D.hh"
#include "G4Vector3D.hh"

#include "globals.hh"
#include "G4Track.hh"

class G4ErrorFreeTrajParam
{
 public:  // with description

  G4ErrorFreeTrajParam()
   : fInvP(0.), fLambda(0.), fPhi(0.), fYPerp(0.), fZPerp(0.){}
  G4ErrorFreeTrajParam( const G4Point3D& pos, const G4Vector3D& mom );
    // build parameters from position and momentum

  virtual ~G4ErrorFreeTrajParam(){}

  void Update( const G4Track* aTrack );
    // update parameters from G4Track

  friend
    std::ostream& operator<<(std::ostream&, const G4ErrorFreeTrajParam& ts);
  
  // Set and Get methods 

  void SetParameters( const G4Point3D& pos, const G4Vector3D& mom );

  G4Vector3D GetDirection() const { return fDir;}

  G4double GetInvP() const { return fInvP; }
  G4double GetLambda() const { return fLambda; }
  G4double GetPhi() const { return fPhi; }
  G4double GetYPerp() const { return fYPerp; }
  G4double GetZPerp() const { return fZPerp; }

 private:

  G4Vector3D fDir; //direction to which YPerp, ZPerp refer
  G4double fInvP; // inverse of momentum
  G4double fLambda; // 90 - theta angle of direction
  G4double fPhi; // phi angle of direction
  G4double fYPerp; // Y coordinate
  G4double fZPerp; // Z coordinate
};

#endif
