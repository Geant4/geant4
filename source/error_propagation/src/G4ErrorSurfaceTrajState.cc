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
// $Id: G4ErrorSurfaceTrajState.cc 69014 2013-04-15 09:42:51Z gcosmo $
//
// ------------------------------------------------------------
//      GEANT 4 class implementation file 
// ------------------------------------------------------------
//

#include "G4ErrorSurfaceTrajState.hh"
#include "G4ErrorPropagatorData.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4Field.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"

#include "G4ErrorMatrix.hh"

#include <iomanip>

//------------------------------------------------------------------------
G4ErrorSurfaceTrajState::
G4ErrorSurfaceTrajState( const G4String& partType, const G4Point3D& pos,
                         const G4Vector3D& mom, const G4Vector3D& vecU,
                         const G4Vector3D& vecV, const G4ErrorTrajErr& errmat)
  : G4ErrorTrajState( partType, pos, mom, errmat )
{
  Init();
  fTrajParam = G4ErrorSurfaceTrajParam( pos, mom, vecU, vecV );
}


//------------------------------------------------------------------------
G4ErrorSurfaceTrajState::
G4ErrorSurfaceTrajState( const G4String& partType, const G4Point3D& pos,
                         const G4Vector3D& mom, const G4Plane3D& plane,
                         const G4ErrorTrajErr& errmat )
  : G4ErrorTrajState( partType, pos, mom, errmat )
{
  Init();
  fTrajParam = G4ErrorSurfaceTrajParam( pos, mom, plane );

}


//------------------------------------------------------------------------
G4ErrorSurfaceTrajState::
G4ErrorSurfaceTrajState( G4ErrorFreeTrajState& tpSC, const G4Plane3D& plane )
  : G4ErrorTrajState( tpSC.GetParticleType(), tpSC.GetPosition(),
                      tpSC.GetMomentum() )
{
  //  fParticleType = tpSC.GetParticleType();
  //  fPosition = tpSC.GetPosition();
  //  fMomentum = tpSC.GetMomentum();
  fTrajParam = G4ErrorSurfaceTrajParam( fPosition, fMomentum, plane );
  Init();

  //----- Get the error matrix in SC coordinates
  BuildErrorMatrix( tpSC, GetVectorV(), GetVectorW() );
}
 

//------------------------------------------------------------------------
G4ErrorSurfaceTrajState::
G4ErrorSurfaceTrajState( G4ErrorFreeTrajState& tpSC, const G4Vector3D& vecU,
                         const G4Vector3D& vecV , G4ErrorMatrix &transfM)
  : G4ErrorTrajState( tpSC.GetParticleType(), tpSC.GetPosition(),
                      tpSC.GetMomentum() )
{
  Init(); // needed to define charge sign 
  fTrajParam = G4ErrorSurfaceTrajParam( fPosition, fMomentum, vecU, vecV );
  //----- Get the error matrix in SC coordinates
  transfM= BuildErrorMatrix( tpSC, vecU, vecV );
}


//------------------------------------------------------------------------
G4ErrorMatrix  G4ErrorSurfaceTrajState::
BuildErrorMatrix( G4ErrorFreeTrajState& tpSC, const G4Vector3D&,
                  const G4Vector3D& )
{
  G4double sclambda = tpSC.GetParameters().GetLambda();
  G4double scphi = tpSC.GetParameters().GetPhi();
  if( G4ErrorPropagatorData::GetErrorPropagatorData()->GetMode() == G4ErrorMode_PropBackwards ){
    sclambda *= -1;
    scphi += CLHEP::pi;
  }
  G4double cosLambda = std::cos( sclambda );
  G4double sinLambda = std::sin( sclambda );
  G4double sinPhi = std::sin( scphi );
  G4double cosPhi = std::cos( scphi );

#ifdef G4EVERBOSE
  if( iverbose >= 4) G4cout << " PM " << fMomentum.mag() << " pLambda " << sclambda << " pPhi " << scphi << G4endl;
#endif

  G4ThreeVector vTN( cosLambda*cosPhi, cosLambda*sinPhi,sinLambda );
  G4ThreeVector vUN( -sinPhi, cosPhi, 0. );
  G4ThreeVector vVN( -vTN.z()*vUN.y(),vTN.z()*vUN.x(), cosLambda );
  
#ifdef G4EVERBOSE
  if( iverbose >= 4) G4cout << " SC2SD: vTN " << vTN << " vUN " << vUN << " vVN " << vVN << G4endl;
#endif
  G4double UJ = vUN*GetVectorV();
  G4double UK = vUN*GetVectorW();
  G4double VJ = vVN*GetVectorV();
  G4double VK = vVN*GetVectorW();


  //--- Get transformation first
  G4ErrorMatrix transfM(5, 5, 0 );
  //--- Get magnetic field
  const G4Field* field = G4TransportationManager::GetTransportationManager()->GetFieldManager()->GetDetectorField();

  G4Vector3D vectorU = GetVectorV().cross( GetVectorW() );
  G4double T1R = 1. / ( vTN * vectorU );

#ifdef G4EVERBOSE
  if( iverbose >= 4) G4cout << "surf vectors:U,V,W " << vectorU << " " <<  GetVectorV() << " " << GetVectorW() << "  T1R:"<<T1R<<G4endl;
#endif


  if( fCharge != 0 && field ) {
    G4double pos[3]; pos[0] = fPosition.x()*cm; pos[1] = fPosition.y()*cm; pos[2] = fPosition.z()*cm;
    G4double Hd[3];
    field->GetFieldValue( pos, Hd );
    G4ThreeVector H = G4ThreeVector( Hd[0], Hd[1], Hd[2] ) / tesla *10.;  //in kilogauus
    G4double magH = H.mag();
    G4double invP = 1./(fMomentum.mag()/GeV);
    G4double magHM = magH * invP;
    if( magH != 0. ) {
      G4double magHM2 = fCharge / magH;
      G4double Q = -magHM * c_light/ (km/ns);
#ifdef G4EVERBOSE
      if( iverbose >= 4) G4cout << GeV <<  " Q " << Q << " magHM " << magHM << " c_light/(km/ns) " << c_light/(km/ns) << G4endl;      
#endif

      G4double sinz = -H*vUN * magHM2;
      G4double cosz =  H*vVN * magHM2;
      G4double T3R = Q * std::pow(T1R,3);
      G4double UI = vUN * vectorU;
      G4double VI = vVN * vectorU;
#ifdef G4EVERBOSE
      if( iverbose >= 4) {
        G4cout << " T1R " << T1R << " T3R " << T3R << G4endl;
        G4cout << " UI " << UI << " VI " << VI << " vectorU " << vectorU << G4endl;
        G4cout << " UJ " << UJ << " VJ " << VJ << G4endl;
        G4cout << " UK " << UK << " VK " << VK << G4endl;
      }
#endif

      transfM[1][3] = -UI*( VK*cosz-UK*sinz)*T3R;
      transfM[1][4] = -VI*( VK*cosz-UK*sinz)*T3R;
      transfM[2][3] = UI*( VJ*cosz-UJ*sinz)*T3R;
      transfM[2][4] = VI*( VJ*cosz-UJ*sinz)*T3R;
    }
  }

  G4double T2R = T1R * T1R;
  transfM[0][0] = 1.;
  transfM[1][1] = -UK*T2R;
  transfM[1][2] = VK*cosLambda*T2R;
  transfM[2][1] = UJ*T2R;
  transfM[2][2] = -VJ*cosLambda*T2R;
  transfM[3][3] = VK*T1R;
  transfM[3][4] = -UK*T1R;
  transfM[4][3] = -VJ*T1R;
  transfM[4][4] = UJ*T1R;

#ifdef G4EVERBOSE
  if( iverbose >= 4) G4cout << " SC2SD transf matrix A= " << transfM << G4endl;
#endif
  fError = G4ErrorTrajErr( tpSC.GetError().similarity( transfM ) );

#ifdef G4EVERBOSE
  if( iverbose >= 1) G4cout << "G4EP: error matrix SC2SD " << fError << G4endl;
  if( iverbose >= 4) G4cout << "G4ErrorSurfaceTrajState from SC " << *this << G4endl;
#endif

  return transfM; // want to use trasnfM for the reverse transform
}


//------------------------------------------------------------------------
void G4ErrorSurfaceTrajState::Init()
{
 theTSType = G4eTS_OS;
 BuildCharge();
}


//------------------------------------------------------------------------
void G4ErrorSurfaceTrajState::Dump( std::ostream& out ) const
{
  out << *this;
}


//------------------------------------------------------------------------
std::ostream& operator<<(std::ostream& out, const G4ErrorSurfaceTrajState& ts)
{
  std::ios::fmtflags oldFlags = out.flags();
  out.setf(std::ios::fixed,std::ios::floatfield);
  
  ts.DumpPosMomError( out );
 
  out << " G4ErrorSurfaceTrajState: Params: " << ts.fTrajParam << G4endl;
  out.flags(oldFlags);
  return out;
}
