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
// $Id: G4ErrorFreeTrajParam.cc,v 1.3 2007-09-24 16:24:45 arce Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ------------------------------------------------------------
//      GEANT 4 class implementation file 
// ------------------------------------------------------------
//

#include "G4ErrorFreeTrajParam.hh"
#include "G4ThreeVector.hh"

#include <iomanip>

//------------------------------------------------------------------------
G4ErrorFreeTrajParam::G4ErrorFreeTrajParam( const G4Point3D& pos,
                                            const G4Vector3D& mom )
{
  SetParameters( pos, mom );
}


//------------------------------------------------------------------------
void G4ErrorFreeTrajParam::SetParameters( const G4Point3D& pos,
                                          const G4Vector3D& mom )
{
  fDir = mom;
  fInvP = 1./mom.mag();
  fLambda = 90.*deg - mom.theta();
  fPhi = mom.phi();
  G4Vector3D vxPerp(0.,0.,0.);
  if( mom.mag() > 0.) {
    vxPerp = mom/mom.mag();
  }
  G4Vector3D vyPerp = G4Vector3D( -vxPerp.y(), vxPerp.x(), 0.);
  vyPerp /= vyPerp.mag();
  G4Vector3D vzPerp = vxPerp.cross( vyPerp );
  vzPerp /= vzPerp.mag();
  // check if right handed
  //  fXPerp = pos.proj( mom );
  G4ThreeVector posv(pos);
  if( vyPerp.mag() != 0. ) {
    fYPerp = posv.project( vyPerp ).mag();
    fZPerp = posv.project( vzPerp ).mag();
  } else {
    fYPerp = 0.;
    fZPerp = 0.;
  }
}

//------------------------------------------------------------------------
void G4ErrorFreeTrajParam::Update( const G4Track* aTrack )
{
  SetParameters( aTrack->GetPosition(), aTrack->GetMomentum() );

}


//------------------------------------------------------------------------
std::ostream& operator<<(std::ostream& out, const G4ErrorFreeTrajParam& tp)
{
  //  long mode = out.setf(std::ios::fixed,std::ios::floatfield);
  
  //  out << tp.theType;
  //  out << std::setprecision(5) << std::setw(10);
  out << std::setprecision(8) << " InvP= " << tp.fInvP << " Theta= "
      << tp.fLambda << " Phi= " << tp.fPhi << " YPerp= " << tp.fYPerp
      << " ZPerp= " << tp.fZPerp << G4endl;
  out << " momentum direction= " << tp.fDir << G4endl;
    
  return out;
}
