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
// ------------------------------------------------------------
//      GEANT 4 class implementation file
// ------------------------------------------------------------
//

#include <iomanip>

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4Field.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4GeometryTolerance.hh"
#include "G4Material.hh"
#include "G4ErrorPropagatorData.hh"
#include "G4ErrorFreeTrajState.hh"
#include "G4ErrorFreeTrajParam.hh"
#include "G4ErrorSurfaceTrajState.hh"
#include "G4ErrorMatrix.hh"

//------------------------------------------------------------------------
G4ErrorFreeTrajState::G4ErrorFreeTrajState(const G4String& partName,
                                           const G4Point3D& pos,
                                           const G4Vector3D& mom,
                                           const G4ErrorTrajErr& errmat)
  : G4ErrorTrajState(partName, pos, mom, errmat)
{
  fTrajParam = G4ErrorFreeTrajParam(pos, mom);
  Init();
}

//------------------------------------------------------------------------
G4ErrorFreeTrajState::G4ErrorFreeTrajState(const G4ErrorSurfaceTrajState& tpSD)
  : G4ErrorTrajState(tpSD.GetParticleType(), tpSD.GetPosition(),
                     tpSD.GetMomentum())
{
  // G4ThreeVector planeNormal = tpSD.GetPlaneNormal();
  // G4double fPt = tpSD.GetMomentum()*planeNormal;//mom projected on normal to
  // plane G4ErrorSurfaceTrajParam tpSDparam = tpSD.GetParameters();
  // G4ThreeVector Psc = fPt * planeNormal +
  // tpSDparam.GetPU()*tpSDparam.GetVectorU() + tpSD.GetPV()*tpSD.GetVectorW();

  fTrajParam = G4ErrorFreeTrajParam(fPosition, fMomentum);
  Init();

  //----- Get the error matrix in SC coordinates
  G4ErrorSurfaceTrajParam tpSDparam = tpSD.GetParameters();
  G4double mom                      = fMomentum.mag();
  G4double mom2                     = fMomentum.mag2();
  G4double TVW1 =
    std::sqrt(mom2 / (mom2 + tpSDparam.GetPV() * tpSDparam.GetPV() +
                      tpSDparam.GetPW() * tpSDparam.GetPW()));
  G4ThreeVector vTVW(TVW1, tpSDparam.GetPV() / mom * TVW1,
                     tpSDparam.GetPW() / mom * TVW1);
  G4Vector3D vectorU = tpSDparam.GetVectorV().cross(tpSDparam.GetVectorW());
  G4Vector3D vTN     = vTVW.x() * vectorU + vTVW.y() * tpSDparam.GetVectorV() +
                   vTVW.z() * tpSDparam.GetVectorW();

#ifdef G4EVERBOSE
  if(iverbose >= 5)
  {
    G4double pc2 = std::asin(vTN.z());
    G4double pc3 = std::atan(vTN.y() / vTN.x());

    G4cout << " CHECK: pc2 " << pc2 << " = " << GetParameters().GetLambda()
           << " diff " << pc2 - GetParameters().GetLambda() << G4endl;
    G4cout << " CHECK: pc3 " << pc3 << " = " << GetParameters().GetPhi()
           << " diff " << pc3 - GetParameters().GetPhi() << G4endl;
  }
#endif

  //--- Get the unit vectors perp to P
  G4double cosl = std::cos(GetParameters().GetLambda());
  if(cosl < 1.E-30)
    cosl = 1.E-30;
  G4double cosl1 = 1. / cosl;
  G4Vector3D vUN(-vTN.y() * cosl1, vTN.x() * cosl1, 0.);
  G4Vector3D vVN(-vTN.z() * vUN.y(), vTN.z() * vUN.x(), cosl);

  G4Vector3D vUperp = G4Vector3D(-fMomentum.y(), fMomentum.x(), 0.);
  G4Vector3D vVperp = vUperp.cross(fMomentum);
  vUperp *= 1. / vUperp.mag();
  vVperp *= 1. / vVperp.mag();

#ifdef G4EVERBOSE
  if(iverbose >= 5)
  {
    G4cout << " CHECK: vUN " << vUN << " = " << vUperp << " diff "
           << (vUN - vUperp).mag() << G4endl;
    G4cout << " CHECK: vVN " << vVN << " = " << vVperp << " diff "
           << (vVN - vVperp).mag() << G4endl;
  }
#endif

  // get the dot products of vectors perpendicular to direction and vector
  // defining SD plane
  G4double dUU = vUperp * tpSD.GetVectorV();
  G4double dUV = vUperp * tpSD.GetVectorW();
  G4double dVU = vVperp * tpSD.GetVectorV();
  G4double dVV = vVperp * tpSD.GetVectorW();

  //--- Get transformation first
  G4ErrorMatrix transfM(5, 5, 1);
  //--- Get magnetic field
  const G4Field* field = G4TransportationManager::GetTransportationManager()
                           ->GetFieldManager()
                           ->GetDetectorField();
  G4ThreeVector dir    = fTrajParam.GetDirection();
  G4double invCosTheta = 1. / std::cos(dir.theta());
  G4cout << " dir=" << dir << " invCosTheta " << invCosTheta << G4endl;

  if(fCharge != 0 && field)
  {
    G4double pos1[3];
    pos1[0] = fPosition.x() * cm;
    pos1[1] = fPosition.y() * cm;
    pos1[2] = fPosition.z() * cm;
    G4double h1[3];
    field->GetFieldValue(pos1, h1);
    G4ThreeVector HPre = G4ThreeVector(h1[0], h1[1], h1[2]) / tesla * 10.;
    G4double magHPre   = HPre.mag();
    G4double invP      = 1. / fMomentum.mag();
    G4double magHPreM  = magHPre * invP;
    if(magHPre != 0.)
    {
      G4double magHPreM2 = fCharge / magHPre;

      G4double Q    = -magHPreM * c_light;
      G4double sinz = -HPre * vUperp * magHPreM2;
      G4double cosz = HPre * vVperp * magHPreM2;

      transfM[1][3] = -Q * dir.y() * sinz;
      transfM[1][4] = -Q * dir.z() * sinz;
      transfM[2][3] = -Q * dir.y() * cosz * invCosTheta;
      transfM[2][4] = -Q * dir.z() * cosz * invCosTheta;
    }
  }

  transfM[0][0] = 1.;
  transfM[1][1] = dir.x() * dVU;
  transfM[1][2] = dir.x() * dVV;
  transfM[2][1] = dir.x() * dUU * invCosTheta;
  transfM[2][2] = dir.x() * dUV * invCosTheta;
  transfM[3][3] = dUU;
  transfM[3][4] = dUV;
  transfM[4][3] = dVU;
  transfM[4][4] = dVV;

  fError = G4ErrorTrajErr(tpSD.GetError().similarity(transfM));

#ifdef G4EVERBOSE
  if(iverbose >= 1)
    G4cout << "error matrix SD2SC " << fError << G4endl;
  if(iverbose >= 4)
    G4cout << "G4ErrorFreeTrajState from SD " << *this << G4endl;
#endif
}

//------------------------------------------------------------------------
void G4ErrorFreeTrajState::Init()
{
  theTSType = G4eTS_FREE;
  BuildCharge();
  theTransfMat = G4ErrorMatrix(5, 5, 0);
  theFirstStep = true;
}

//------------------------------------------------------------------------
void G4ErrorFreeTrajState::Dump(std::ostream& out) const { out << *this; }

//------------------------------------------------------------------------
G4int G4ErrorFreeTrajState::Update(const G4Track* aTrack)
{
  G4int ierr = 0;
  fTrajParam.Update(aTrack);
  UpdatePosMom(aTrack->GetPosition(), aTrack->GetMomentum());
  return ierr;
}

//------------------------------------------------------------------------
std::ostream& operator<<(std::ostream& out, const G4ErrorFreeTrajState& ts)
{
  std::ios::fmtflags orig_flags = out.flags();

  out.setf(std::ios::fixed, std::ios::floatfield);

  ts.DumpPosMomError(out);

  out << " G4ErrorFreeTrajState: Params: " << ts.fTrajParam << G4endl;

  out.flags(orig_flags);

  return out;
}

//------------------------------------------------------------------------
G4int G4ErrorFreeTrajState::PropagateError(const G4Track* aTrack)
{
  G4double stepLengthCm = aTrack->GetStep()->GetStepLength() / cm;
  if(G4ErrorPropagatorData::GetErrorPropagatorData()->GetStage() ==
     G4ErrorStage_Deflation)
    stepLengthCm *= -1.;

  G4double kCarTolerance =
    G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();

  if(std::fabs(stepLengthCm) <= kCarTolerance / cm)
    return 0;

#ifdef G4EVERBOSE
  if(iverbose >= 2)
    G4cout << "  G4ErrorFreeTrajState::PropagateError " << G4endl;
  G4cout << "G4EP: iverbose=" << iverbose << G4endl;
#endif

  // * *** ERROR PROPAGATION ON A HELIX ASSUMING SC VARIABLES
  G4Point3D vposPost = aTrack->GetPosition() / cm;
  G4Vector3D vpPost  = aTrack->GetMomentum() / GeV;
  //  G4Point3D vposPre = fPosition/cm;
  //  G4Vector3D vpPre = fMomentum/GeV;
  G4Point3D vposPre = aTrack->GetStep()->GetPreStepPoint()->GetPosition() / cm;
  G4Vector3D vpPre  = aTrack->GetStep()->GetPreStepPoint()->GetMomentum() / GeV;
  // correct to avoid propagation along Z
  if(vpPre.mag() == vpPre.z())
    vpPre.setX(1.E-6 * MeV);
  if(vpPost.mag() == vpPost.z())
    vpPost.setX(1.E-6 * MeV);

  G4double pPre  = vpPre.mag();
  G4double pPost = vpPost.mag();
#ifdef G4EVERBOSE
  if(iverbose >= 2)
  {
    G4cout << "G4EP: vposPre " << vposPre << G4endl << "G4EP: vposPost "
           << vposPost << G4endl;
    G4cout << "G4EP: vpPre " << vpPre << G4endl << "G4EP: vpPost " << vpPost
           << G4endl;
    G4cout << " err start step " << fError << G4endl;
    G4cout << "G4EP: stepLengthCm " << stepLengthCm << G4endl;
  }
#endif

  if(pPre == 0. || pPost == 0)
    return 2;
  G4double pInvPre   = 1. / pPre;
  G4double pInvPost  = 1. / pPost;
  G4double deltaPInv = pInvPost - pInvPre;
  if(iverbose >= 2)
    G4cout << "G4EP:  pInvPre" << pInvPre << "  pInvPost:" << pInvPost
           << "  deltaPInv:" << deltaPInv << G4endl;

  G4Vector3D vpPreNorm  = vpPre * pInvPre;
  G4Vector3D vpPostNorm = vpPost * pInvPost;
  if(iverbose >= 2)
    G4cout << "G4EP: vpPreNorm " << vpPreNorm << " vpPostNorm " << vpPostNorm
           << G4endl;
  // return if propagation along Z??
  if(1. - std::fabs(vpPreNorm.z()) < kCarTolerance)
    return 4;
  if(1. - std::fabs(vpPostNorm.z()) < kCarTolerance)
    return 4;
  G4double sinpPre =
    std::sin(vpPreNorm.theta());  // cosine perpendicular to pPre = sine pPre
  G4double sinpPost =
    std::sin(vpPostNorm.theta());  // cosine perpendicular to pPost = sine pPost
  G4double sinpPostInv = 1. / std::sin(vpPostNorm.theta());

#ifdef G4EVERBOSE
  if(iverbose >= 2)
    G4cout << "G4EP: cosl " << sinpPre << " cosl0 " << sinpPost << G4endl;
#endif
  //* *** DEFINE TRANSFORMATION MATRIX BETWEEN X1 AND X2 FOR
  //* *** NEUTRAL PARTICLE OR FIELDFREE REGION
  G4ErrorMatrix transf(5, 5, 0);

  transf[3][2] = stepLengthCm * sinpPost;
  transf[4][1] = stepLengthCm;
  for(auto ii = 0; ii < 5; ++ii)
  {
    transf[ii][ii] = 1.;
  }
#ifdef G4EVERBOSE
  if(iverbose >= 2)
  {
    G4cout << "G4EP: transf matrix neutral " << transf;
  }
#endif

  //  charge X propagation direction
  G4double charge = aTrack->GetDynamicParticle()->GetCharge();
  if(G4ErrorPropagatorData::GetErrorPropagatorData()->GetMode() ==
     G4ErrorMode_PropBackwards)
  {
    charge *= -1.;
  }
  //  G4cout << " charge " << charge << G4endl;
  // t check if particle has charge
  // t  if( charge == 0 ) goto 45;
  // check if the magnetic field is = 0.

  // position is from geant4, it is assumed to be in mm (for debugging,
  // eventually it will not be transformed) it is assumed vposPre[] is in cm and
  // pos1[] is in mm.
  G4double pos1[3];
  pos1[0] = vposPre.x() * cm;
  pos1[1] = vposPre.y() * cm;
  pos1[2] = vposPre.z() * cm;
  G4double pos2[3];
  pos2[0] = vposPost.x() * cm;
  pos2[1] = vposPost.y() * cm;
  pos2[2] = vposPost.z() * cm;
  G4double h1[3], h2[3];

  const G4Field* field = G4TransportationManager::GetTransportationManager()
                           ->GetFieldManager()
                           ->GetDetectorField();
  if(!field)
    return 0;  // goto 45

  // calculate transformation except it NEUTRAL PARTICLE OR FIELDFREE REGION
  if(charge != 0. && field)
  {
    field->GetFieldValue(pos1,
                         h1);  // here pos1[], pos2[] are in mm, not changed
    field->GetFieldValue(pos2, h2);
    G4ThreeVector HPre =
      G4ThreeVector(h1[0], h1[1], h1[2]) / tesla *
      10.;  // 10. is to get same dimensions as GEANT3 (kilogauss)
    G4ThreeVector HPost = G4ThreeVector(h2[0], h2[1], h2[2]) / tesla * 10.;
    G4double magHPre    = HPre.mag();
    G4double magHPost   = HPost.mag();
#ifdef G4EVERBOSE
    if(iverbose >= 2)
    {
      G4cout << "G4EP: h1 = " << h1[0] << ", " << h1[1] << ", " << h1[2]
             << G4endl;
      G4cout << "G4EP: pos1/mm = " << pos1[0] << ", " << pos1[1] << ", "
             << pos1[2] << G4endl;
      G4cout << "G4EP: pos2/mm = " << pos2[0] << ", " << pos2[1] << ", "
             << pos2[2] << G4endl;
      G4cout << "G4EP: B-filed in KGauss HPre " << HPre << G4endl
             << "G4EP: in KGauss HPost " << HPost << G4endl;
    }
#endif

    if(magHPre + magHPost != 0.)
    {
      //* *** CHECK WHETHER H*ALFA/P IS TOO DIFFERENT AT X1 AND X2
      G4double gam;
      if(magHPost != 0.)
      {
        gam = HPost * vpPostNorm / magHPost;
      }
      else
      {
        gam = HPre * vpPreNorm / magHPre;
      }

      // G4eMagneticLimitsProcess will limit the step, but based on an straight
      // line trajectory
      G4double alphaSqr  = 1. - gam * gam;
      G4double diffHSqr  = (HPre * pInvPre - HPost * pInvPost).mag2();
      G4double delhp6Sqr = 300. * 300.;
#ifdef G4EVERBOSE
      if(iverbose >= 2)
      {
        G4cout << " G4EP: gam " << gam << " alphaSqr " << alphaSqr
               << " diffHSqr " << diffHSqr << G4endl;
        G4cout << " alpha= " << std::sqrt(alphaSqr) << G4endl;
      }
#endif
      if(diffHSqr * alphaSqr > delhp6Sqr)
        return 3;

      //* *** DEFINE AVERAGE MAGNETIC FIELD AND GRADIENT
      G4double pInvAver = 1. / (pInvPre + pInvPost);
      G4double CFACT8   = 2.997925E-4;
      // G4double HAver
      G4ThreeVector vHAverNorm((HPre * pInvPre + HPost * pInvPost) * pInvAver *
                               charge * CFACT8);
      G4double HAver    = vHAverNorm.mag();
      G4double invHAver = 1. / HAver;
      vHAverNorm *= invHAver;
#ifdef G4EVERBOSE
      if(iverbose >= 2)
        G4cout << " G4EP: HaverNorm " << vHAverNorm << " magHAver " << HAver
               << " charge " << charge << G4endl;
#endif

      G4double pAver        = (pPre + pPost) * 0.5;
      G4double QAver        = -HAver / pAver;
      G4double thetaAver    = QAver * stepLengthCm;
      G4double sinThetaAver = std::sin(thetaAver);
      G4double cosThetaAver = std::cos(thetaAver);
      G4double gamma        = vHAverNorm * vpPostNorm;
      G4ThreeVector AN2     = vHAverNorm.cross(vpPostNorm);

#ifdef G4EVERBOSE
      if(iverbose >= 2)
        G4cout << " G4EP: AN2 " << AN2 << "  gamma:" << gamma
               << "  theta=" << thetaAver << G4endl;
#endif
      G4double AU = 1. / vpPreNorm.perp();
      // t  G4ThreeVector vU( vpPreNorm.cross( G4ThreeVector(0.,0.,1.) ) * AU );
      G4ThreeVector vUPre(-AU * vpPreNorm.y(), AU * vpPreNorm.x(), 0.);
      G4ThreeVector vVPre(-vpPreNorm.z() * vUPre.y(), vpPreNorm.z() * vUPre.x(),
                          vpPreNorm.x() * vUPre.y() -
                            vpPreNorm.y() * vUPre.x());

      //
      AU = 1. / vpPostNorm.perp();
      // t  G4ThreeVector vU( vpPostNorm.cross( G4ThreeVector(0.,0.,1.) ) * AU
      // );
      G4ThreeVector vUPost(-AU * vpPostNorm.y(), AU * vpPostNorm.x(), 0.);
      G4ThreeVector vVPost(
        -vpPostNorm.z() * vUPost.y(), vpPostNorm.z() * vUPost.x(),
        vpPostNorm.x() * vUPost.y() - vpPostNorm.y() * vUPost.x());
#ifdef G4EVERBOSE
      G4cout << " vpPostNorm " << vpPostNorm << G4endl;
      if(iverbose >= 2)
        G4cout << " G4EP: AU " << AU << " vUPre " << vUPre << " vVPre " << vVPre
               << " vUPost " << vUPost << " vVPost " << vVPost << G4endl;
#endif
      G4Point3D deltaPos(vposPre - vposPost);

      // * *** COMPLETE TRANSFORMATION MATRIX BETWEEN ERRORS AT X1 AND X2
      // * *** FIELD GRADIENT PERPENDICULAR TO TRACK IS PRESENTLY NOT
      // * *** TAKEN INTO ACCOUNT

      G4double QP = QAver * pAver;  // = -HAver
#ifdef G4EVERBOSE
      if(iverbose >= 2)
        G4cout << " G4EP: QP " << QP << " QAver " << QAver << " pAver " << pAver
               << G4endl;
#endif
      G4double ANV =
        -(vHAverNorm.x() * vUPost.x() + vHAverNorm.y() * vUPost.y());
      G4double ANU =
        (vHAverNorm.x() * vVPost.x() + vHAverNorm.y() * vVPost.y() +
         vHAverNorm.z() * vVPost.z());
      G4double OMcosThetaAver = 1. - cosThetaAver;
#ifdef G4EVERBOSE
      if(iverbose >= 2)
        G4cout << "G4EP: OMcosThetaAver " << OMcosThetaAver << " cosThetaAver "
               << cosThetaAver << " thetaAver " << thetaAver << " QAver "
               << QAver << " stepLengthCm " << stepLengthCm << G4endl;
#endif
      G4double TMSINT = thetaAver - sinThetaAver;
#ifdef G4EVERBOSE
      if(iverbose >= 2)
        G4cout << " G4EP: ANV " << ANV << " ANU " << ANU << G4endl;
#endif

      G4ThreeVector vHUPre(
        -vHAverNorm.z() * vUPre.y(), vHAverNorm.z() * vUPre.x(),
        vHAverNorm.x() * vUPre.y() - vHAverNorm.y() * vUPre.x());
#ifdef G4EVERBOSE
      //    if( iverbose >= 2 ) G4cout << "G4EP: HUPre(1) " << vHUPre.x() << " "
      //    << vHAverNorm.z() << " " << vUPre.y() << G4endl;
#endif
      G4ThreeVector vHVPre(
        vHAverNorm.y() * vVPre.z() - vHAverNorm.z() * vVPre.y(),
        vHAverNorm.z() * vVPre.x() - vHAverNorm.x() * vVPre.z(),
        vHAverNorm.x() * vVPre.y() - vHAverNorm.y() * vVPre.x());
#ifdef G4EVERBOSE
      if(iverbose >= 2)
        G4cout << " G4EP: HUPre " << vHUPre << " HVPre " << vHVPre << G4endl;
#endif

      //------------------- COMPUTE MATRIX
      //---------- 1/P

      transf[0][0] =
        1. -
        deltaPInv * pAver *
          (1. + (vpPostNorm.x() * deltaPos.x() + vpPostNorm.y() * deltaPos.y() +
                 vpPostNorm.z() * deltaPos.z()) /
                  stepLengthCm) +
        2. * deltaPInv * pAver;

      transf[0][1] =
        -deltaPInv / thetaAver *
        (TMSINT * gamma *
           (vHAverNorm.x() * vVPre.x() + vHAverNorm.y() * vVPre.y() +
            vHAverNorm.z() * vVPre.z()) +
         sinThetaAver *
           (vVPre.x() * vpPostNorm.x() + vVPre.y() * vpPostNorm.y() +
            vVPre.z() * vpPostNorm.z()) +
         OMcosThetaAver *
           (vHVPre.x() * vpPostNorm.x() + vHVPre.y() * vpPostNorm.y() +
            vHVPre.z() * vpPostNorm.z()));

      transf[0][2] =
        -sinpPre * deltaPInv / thetaAver *
        (TMSINT * gamma *
           (vHAverNorm.x() * vUPre.x() + vHAverNorm.y() * vUPre.y()) +
         sinThetaAver *
           (vUPre.x() * vpPostNorm.x() + vUPre.y() * vpPostNorm.y()) +
         OMcosThetaAver *
           (vHUPre.x() * vpPostNorm.x() + vHUPre.y() * vpPostNorm.y() +
            vHUPre.z() * vpPostNorm.z()));

      transf[0][3] = -deltaPInv / stepLengthCm *
                     (vUPre.x() * vpPostNorm.x() + vUPre.y() * vpPostNorm.y());

      transf[0][4] = -deltaPInv / stepLengthCm *
                     (vVPre.x() * vpPostNorm.x() + vVPre.y() * vpPostNorm.y() +
                      vVPre.z() * vpPostNorm.z());

      // ***   Lambda
      transf[1][0] =
        -QP * ANV *
        (vpPostNorm.x() * deltaPos.x() + vpPostNorm.y() * deltaPos.y() +
         vpPostNorm.z() * deltaPos.z()) *
        (1. + deltaPInv * pAver);
#ifdef G4EVERBOSE
      if(iverbose >= 3)
        G4cout << "ctransf10= " << transf[1][0] << " " << -QP << " " << ANV
               << " " << vpPostNorm.x() << " " << deltaPos.x() << " "
               << vpPostNorm.y() << " " << deltaPos.y() << " " << vpPostNorm.z()
               << " " << deltaPos.z() << " " << deltaPInv << " " << pAver
               << G4endl;
#endif

      transf[1][1] =
        cosThetaAver * (vVPre.x() * vVPost.x() + vVPre.y() * vVPost.y() +
                        vVPre.z() * vVPost.z()) +
        sinThetaAver * (vHVPre.x() * vVPost.x() + vHVPre.y() * vVPost.y() +
                        vHVPre.z() * vVPost.z()) +
        OMcosThetaAver *
          (vHAverNorm.x() * vVPre.x() + vHAverNorm.y() * vVPre.y() +
           vHAverNorm.z() * vVPre.z()) *
          (vHAverNorm.x() * vVPost.x() + vHAverNorm.y() * vVPost.y() +
           vHAverNorm.z() * vVPost.z()) +
        ANV * (-sinThetaAver *
                 (vVPre.x() * vpPostNorm.x() + vVPre.y() * vpPostNorm.y() +
                  vVPre.z() * vpPostNorm.z()) +
               OMcosThetaAver * (vVPre.x() * AN2.x() + vVPre.y() * AN2.y() +
                                 vVPre.z() * AN2.z()) -
               TMSINT * gamma *
                 (vHAverNorm.x() * vVPre.x() + vHAverNorm.y() * vVPre.y() +
                  vHAverNorm.z() * vVPre.z()));

      transf[1][2] =
        cosThetaAver * (vUPre.x() * vVPost.x() + vUPre.y() * vVPost.y()) +
        sinThetaAver * (vHUPre.x() * vVPost.x() + vHUPre.y() * vVPost.y() +
                        vHUPre.z() * vVPost.z()) +
        OMcosThetaAver *
          (vHAverNorm.x() * vUPre.x() + vHAverNorm.y() * vUPre.y()) *
          (vHAverNorm.x() * vVPost.x() + vHAverNorm.y() * vVPost.y() +
           vHAverNorm.z() * vVPost.z()) +
        ANV * (-sinThetaAver *
                 (vUPre.x() * vpPostNorm.x() + vUPre.y() * vpPostNorm.y()) +
               OMcosThetaAver * (vUPre.x() * AN2.x() + vUPre.y() * AN2.y()) -
               TMSINT * gamma *
                 (vHAverNorm.x() * vUPre.x() + vHAverNorm.y() * vUPre.y()));
      transf[1][2] = sinpPre * transf[1][2];

      transf[1][3] = -QAver * ANV *
                     (vUPre.x() * vpPostNorm.x() + vUPre.y() * vpPostNorm.y());

      transf[1][4] = -QAver * ANV *
                     (vVPre.x() * vpPostNorm.x() + vVPre.y() * vpPostNorm.y() +
                      vVPre.z() * vpPostNorm.z());

      // ***   Phi

      transf[2][0] =
        -QP * ANU *
        (vpPostNorm.x() * deltaPos.x() + vpPostNorm.y() * deltaPos.y() +
         vpPostNorm.z() * deltaPos.z()) *
        sinpPostInv * (1. + deltaPInv * pAver);
#ifdef G4EVERBOSE
      if(iverbose >= 3)
        G4cout << "ctransf20= " << transf[2][0] << " " << -QP << " " << ANU
               << " " << vpPostNorm.x() << " " << deltaPos.x() << " "
               << vpPostNorm.y() << " " << deltaPos.y() << " " << vpPostNorm.z()
               << " " << deltaPos.z() << " " << sinpPostInv << " " << deltaPInv
               << " " << pAver << G4endl;
#endif
      transf[2][1] =
        cosThetaAver * (vVPre.x() * vUPost.x() + vVPre.y() * vUPost.y()) +
        sinThetaAver * (vHVPre.x() * vUPost.x() + vHVPre.y() * vUPost.y()) +
        OMcosThetaAver *
          (vHAverNorm.x() * vVPre.x() + vHAverNorm.y() * vVPre.y() +
           vHAverNorm.z() * vVPre.z()) *
          (vHAverNorm.x() * vUPost.x() + vHAverNorm.y() * vUPost.y()) +
        ANU * (-sinThetaAver *
                 (vVPre.x() * vpPostNorm.x() + vVPre.y() * vpPostNorm.y() +
                  vVPre.z() * vpPostNorm.z()) +
               OMcosThetaAver * (vVPre.x() * AN2.x() + vVPre.y() * AN2.y() +
                                 vVPre.z() * AN2.z()) -
               TMSINT * gamma *
                 (vHAverNorm.x() * vVPre.x() + vHAverNorm.y() * vVPre.y() +
                  vHAverNorm.z() * vVPre.z()));
      transf[2][1] = sinpPostInv * transf[2][1];

      transf[2][2] =
        cosThetaAver * (vUPre.x() * vUPost.x() + vUPre.y() * vUPost.y()) +
        sinThetaAver * (vHUPre.x() * vUPost.x() + vHUPre.y() * vUPost.y()) +
        OMcosThetaAver *
          (vHAverNorm.x() * vUPre.x() + vHAverNorm.y() * vUPre.y()) *
          (vHAverNorm.x() * vUPost.x() + vHAverNorm.y() * vUPost.y()) +
        ANU * (-sinThetaAver *
                 (vUPre.x() * vpPostNorm.x() + vUPre.y() * vpPostNorm.y()) +
               OMcosThetaAver * (vUPre.x() * AN2.x() + vUPre.y() * AN2.y()) -
               TMSINT * gamma *
                 (vHAverNorm.x() * vUPre.x() + vHAverNorm.y() * vUPre.y()));
      transf[2][2] = sinpPostInv * sinpPre * transf[2][2];

      transf[2][3] = -QAver * ANU *
                     (vUPre.x() * vpPostNorm.x() + vUPre.y() * vpPostNorm.y()) *
                     sinpPostInv;
#ifdef G4EVERBOSE
      if(iverbose >= 3)
        G4cout << "ctransf23= " << transf[2][3] << " " << -QAver << " " << ANU
               << " " << vUPre.x() << " " << vpPostNorm.x() << " " << vUPre.y()
               << " " << vpPostNorm.y() << " " << sinpPostInv << G4endl;
#endif

      transf[2][4] = -QAver * ANU *
                     (vVPre.x() * vpPostNorm.x() + vVPre.y() * vpPostNorm.y() +
                      vVPre.z() * vpPostNorm.z()) *
                     sinpPostInv;

      // ***   Yt

      transf[3][0] = pAver *
                     (vUPost.x() * deltaPos.x() + vUPost.y() * deltaPos.y()) *
                     (1. + deltaPInv * pAver);
#ifdef G4EVERBOSE
      if(iverbose >= 3)
        G4cout << "ctransf30= " << transf[3][0] << " " << pAver << " "
               << vUPost.x() << " " << deltaPos.x() << " " << vUPost.y() << " "
               << deltaPos.y() << " " << deltaPInv << " " << pAver << G4endl;
#endif

      transf[3][1] =
        (sinThetaAver * (vVPre.x() * vUPost.x() + vVPre.y() * vUPost.y()) +
         OMcosThetaAver * (vHVPre.x() * vUPost.x() + vHVPre.y() * vUPost.y()) +
         TMSINT * (vHAverNorm.x() * vUPost.x() + vHAverNorm.y() * vUPost.y()) *
           (vHAverNorm.x() * vVPre.x() + vHAverNorm.y() * vVPre.y() +
            vHAverNorm.z() * vVPre.z())) /
        QAver;

      transf[3][2] =
        (sinThetaAver * (vUPre.x() * vUPost.x() + vUPre.y() * vUPost.y()) +
         OMcosThetaAver * (vHUPre.x() * vUPost.x() + vHUPre.y() * vUPost.y()) +
         TMSINT * (vHAverNorm.x() * vUPost.x() + vHAverNorm.y() * vUPost.y()) *
           (vHAverNorm.x() * vUPre.x() + vHAverNorm.y() * vUPre.y())) *
        sinpPre / QAver;
#ifdef G4EVERBOSE
      if(iverbose >= 3)
        G4cout << "ctransf32= " << transf[3][2] << " " << sinThetaAver << " "
               << vUPre.x() << " " << vUPost.x() << " " << vUPre.y() << " "
               << vUPost.y() << " " << OMcosThetaAver << " " << vHUPre.x()
               << " " << vUPost.x() << " " << vHUPre.y() << " " << vUPost.y()
               << " " << TMSINT << " " << vHAverNorm.x() << " " << vUPost.x()
               << " " << vHAverNorm.y() << " " << vUPost.y() << " "
               << vHAverNorm.x() << " " << vUPre.x() << " " << vHAverNorm.y()
               << " " << vUPre.y() << " " << sinpPre << " " << QAver << G4endl;
#endif

      transf[3][3] = (vUPre.x() * vUPost.x() + vUPre.y() * vUPost.y());

      transf[3][4] = (vVPre.x() * vUPost.x() + vVPre.y() * vUPost.y());

      // ***   Zt
      transf[4][0] = pAver *
                     (vVPost.x() * deltaPos.x() + vVPost.y() * deltaPos.y() +
                      vVPost.z() * deltaPos.z()) *
                     (1. + deltaPInv * pAver);

      transf[4][1] =
        (sinThetaAver * (vVPre.x() * vVPost.x() + vVPre.y() * vVPost.y() +
                         vVPre.z() * vVPost.z()) +
         OMcosThetaAver * (vHVPre.x() * vVPost.x() + vHVPre.y() * vVPost.y() +
                           vHVPre.z() * vVPost.z()) +
         TMSINT *
           (vHAverNorm.x() * vVPost.x() + vHAverNorm.y() * vVPost.y() +
            vHAverNorm.z() * vVPost.z()) *
           (vHAverNorm.x() * vVPre.x() + vHAverNorm.y() * vVPre.y() +
            vHAverNorm.z() * vVPre.z())) /
        QAver;
#ifdef G4EVERBOSE
      if(iverbose >= 3)
        G4cout << "ctransf41= " << transf[4][1] << " " << sinThetaAver << " "
               << OMcosThetaAver << " " << TMSINT << " " << vVPre << " "
               << vVPost << " " << vHVPre << " " << vHAverNorm << " " << QAver
               << G4endl;
#endif

      transf[4][2] =
        (sinThetaAver * (vUPre.x() * vVPost.x() + vUPre.y() * vVPost.y()) +
         OMcosThetaAver * (vHUPre.x() * vVPost.x() + vHUPre.y() * vVPost.y() +
                           vHUPre.z() * vVPost.z()) +
         TMSINT *
           (vHAverNorm.x() * vVPost.x() + vHAverNorm.y() * vVPost.y() +
            vHAverNorm.z() * vVPost.z()) *
           (vHAverNorm.x() * vUPre.x() + vHAverNorm.y() * vUPre.y())) *
        sinpPre / QAver;

      transf[4][3] = (vUPre.x() * vVPost.x() + vUPre.y() * vVPost.y());

      transf[4][4] = (vVPre.x() * vVPost.x() + vVPre.y() * vVPost.y() +
                      vVPre.z() * vVPost.z());
      //   if(iverbose >= 3) G4cout <<"ctransf44= " << transf[4][4] <<" "<<
      //   vVPre.x()  <<" "<<vVPost.x() <<" "<< vVPre.y() <<" "<< vVPost.y() <<"
      //   "<< vVPre.z() <<" "<< vVPost.z() << G4endl;

#ifdef G4EVERBOSE
      if(iverbose >= 1)
        G4cout << "G4EP: transf matrix computed " << transf << G4endl;
#endif
      /*    for( G4int ii=0;ii<5;ii++){
        for( G4int jj=0;jj<5;jj++){
          G4cout << transf[ii][jj] << " ";
        }
        G4cout << G4endl;
        } */
    }
  }
  // end of calculate transformation except it NEUTRAL PARTICLE OR FIELDFREE
  // REGION
  /*  if( iverbose >= 1 ) G4cout << "G4EP: transf not updated but initialized "
  << theFirstStep << G4endl; if( theFirstStep ) { theTransfMat = transf;
    theFirstStep = false;
  }else{
    theTransfMat = theTransfMat * transf;
    if( iverbose >= 1 ) G4cout << "G4EP: transf matrix accumulated" <<
  theTransfMat << G4endl;
  }
  */
  theTransfMat = transf;
#ifdef G4EVERBOSE
  if(iverbose >= 1)
    G4cout << "G4EP: error matrix before transformation " << fError << G4endl;
  if(iverbose >= 2)
    G4cout << " tf * err " << theTransfMat * fError << G4endl
           << " transf matrix " << theTransfMat.T() << G4endl;
#endif

  fError = fError.similarity(theTransfMat).T();
  //-    fError = transf * fError * transf.T();
#ifdef G4EVERBOSE
  if(iverbose >= 1)
    G4cout << "G4EP: error matrix propagated " << fError << G4endl;
#endif

  //? S = B*S*BT S.similarity(B)
  //? R = S
  // not needed * *** TRANSFORM ERROR MATRIX FROM INTERNAL TO EXTERNAL
  // VARIABLES;

  PropagateErrorMSC(aTrack);

  PropagateErrorIoni(aTrack);

  return 0;
}

//------------------------------------------------------------------------
G4int G4ErrorFreeTrajState::PropagateErrorMSC(const G4Track* aTrack)
{
  G4ThreeVector vpPre   = aTrack->GetMomentum() / GeV;
  G4double pPre         = vpPre.mag();
  G4double pBeta        = pPre * pPre / (aTrack->GetTotalEnergy() / GeV);
  G4double stepLengthCm = aTrack->GetStep()->GetStepLength() / cm;

  G4Material* mate = aTrack->GetVolume()->GetLogicalVolume()->GetMaterial();
  G4double effZ, effA;
  CalculateEffectiveZandA(mate, effZ, effA);

#ifdef G4EVERBOSE
  if(iverbose >= 4)
    G4cout << "material "
           << mate->GetName()
           //<< " " << mate->GetZ() << " "  << mate->GetA()
           << " effZ:" << effZ << " effA:" << effA
           << " dens(g/mole):" << mate->GetDensity() / g * mole
           << " Radlen/cm:" << mate->GetRadlen() / cm << " nuclLen/cm"
           << mate->GetNuclearInterLength() / cm << G4endl;
#endif

  G4double RI = stepLengthCm / (mate->GetRadlen() / cm);
#ifdef G4EVERBOSE
  if(iverbose >= 4)
    G4cout << std::setprecision(6) << std::setw(6) << "G4EP:MSC: RI=X/X0 " << RI
           << " stepLengthCm " << stepLengthCm << " radlen/cm "
           << (mate->GetRadlen() / cm) << " RI*1.e10:" << RI * 1.e10 << G4endl;
#endif
  G4double charge = aTrack->GetDynamicParticle()->GetCharge();
  G4double DD     = 1.8496E-4 * RI * (charge / pBeta * charge / pBeta);
#ifdef G4EVERBOSE
  if(iverbose >= 3)
    G4cout << "G4EP:MSC: D*1E6= " << DD * 1.E6 << " pBeta " << pBeta << G4endl;
#endif
  G4double S1 = DD * stepLengthCm * stepLengthCm / 3.;
  G4double S2 = DD;
  G4double S3 = DD * stepLengthCm / 2.;

  G4double CLA =
    std::sqrt(vpPre.x() * vpPre.x() + vpPre.y() * vpPre.y()) / pPre;
#ifdef G4EVERBOSE
  if(iverbose >= 2)
    G4cout << std::setw(6) << "G4EP:MSC: RI " << RI << " S1 " << S1 << " S2 "
           << S2 << " S3 " << S3 << " CLA " << CLA << G4endl;
#endif
  fError[1][1] += S2;
  fError[1][4] -= S3;
  fError[2][2] += S2 / CLA / CLA;
  fError[2][3] += S3 / CLA;
  fError[3][3] += S1;
  fError[4][4] += S1;

#ifdef G4EVERBOSE
  if(iverbose >= 2)
    G4cout << "G4EP:MSC: error matrix propagated msc " << fError << G4endl;
#endif

  return 0;
}

//------------------------------------------------------------------------
void G4ErrorFreeTrajState::CalculateEffectiveZandA(const G4Material* mate,
                                                   G4double& effZ,
                                                   G4double& effA)
{
  effZ = 0.;
  effA = 0.;
  auto nelem = mate->GetNumberOfElements();
  const G4double* fracVec = mate->GetFractionVector();
  for(G4int ii = 0; ii < (G4int)nelem; ++ii)
  {
    effZ += mate->GetElement(ii)->GetZ() * fracVec[ii];
    effA += mate->GetElement(ii)->GetA() * fracVec[ii] / g * mole;
  }
}

//------------------------------------------------------------------------
G4int G4ErrorFreeTrajState::PropagateErrorIoni(const G4Track* aTrack)
{
  G4double stepLengthCm = aTrack->GetStep()->GetStepLength() / cm;
#ifdef G4EVERBOSE
  G4double DEDX2;
  if(stepLengthCm < 1.E-7)
  {
    DEDX2 = 0.;
  }
#endif
  //  *     Calculate xi factor (KeV).
  G4Material* mate = aTrack->GetVolume()->GetLogicalVolume()->GetMaterial();
  G4double effZ, effA;
  CalculateEffectiveZandA(mate, effZ, effA);

  G4double Etot  = aTrack->GetTotalEnergy() / GeV;
  G4double beta  = aTrack->GetMomentum().mag() / GeV / Etot;
  G4double mass  = aTrack->GetDynamicParticle()->GetMass() / GeV;
  G4double gamma = Etot / mass;

  // *     Calculate xi factor (keV).
  G4double XI = 153.5 * effZ * stepLengthCm * (mate->GetDensity() / mg * mole) /
                (effA * beta * beta);

#ifdef G4EVERBOSE
  if(iverbose >= 2)
  {
    G4cout << "G4EP:IONI: XI/keV " << XI << " beta " << beta << " gamma "
           << gamma << G4endl;
    G4cout << " density " << (mate->GetDensity() / mg * mole) << " effA "
           << effA << " step " << stepLengthCm << G4endl;
  }
#endif
  // *     Maximum energy transfer to atomic electron (KeV).
  G4double eta       = beta * gamma;
  G4double etasq     = eta * eta;
  G4double eMass     = 0.51099906 / GeV;
  G4double massRatio = eMass / mass;
  G4double F1        = 2 * eMass * etasq;
  G4double F2        = 1. + 2. * massRatio * gamma + massRatio * massRatio;
  G4double Emax      = 1.E+6 * F1 / F2;  // now in keV

  //  * *** and now sigma**2  in GeV
  G4double dedxSq =
    XI * Emax * (1. - (beta * beta / 2.)) * 1.E-12;  // now in GeV^2
  /*The above  formula for var(1/p) good for dens scatterers. However, for MIPS
    passing through a gas it leads to overestimation. Further more for incident
    electrons the Emax is almost equal to incident energy. This leads  to
    k=Xi/Emax  as small as e-6  and gradually the cov matrix  explodes.

    http://www2.pv.infn.it/~rotondi/kalman_1.pdf

    Since I do not have enough info at the moment to implement Landau &
    sub-Landau models for k=Xi/Emax <0.01 I'll saturate k at this value for now
  */

  if(XI / Emax < 0.01)
    dedxSq *=
      XI / Emax * 100;  // Quench for low Elos, see above: newVar=odVar *k/0.01

#ifdef G4EVERBOSE
  if(iverbose >= 2)
    G4cout << "G4EP:IONI: DEDX^2(GeV^2) " << dedxSq << " emass/GeV: " << eMass
           << " Emax/keV: " << Emax << "  k=Xi/Emax=" << XI / Emax << G4endl;

#endif

  G4double pPre6 =
    (aTrack->GetStep()->GetPreStepPoint()->GetMomentum() / GeV).mag();
  pPre6 = std::pow(pPre6, 6);
  // Apply it to error
  fError[0][0] += Etot * Etot * dedxSq / pPre6;
#ifdef G4EVERBOSE
  if(iverbose >= 2)
    G4cout << "G4:IONI Etot/GeV: " << Etot << " err_dedx^2/GeV^2: " << dedxSq
           << " p^6: " << pPre6 << G4endl;
  if(iverbose >= 2)
    G4cout << "G4EP:IONI: error2_from_ionisation "
           << (Etot * Etot * dedxSq) / pPre6 << G4endl;
#endif

  return 0;
}
