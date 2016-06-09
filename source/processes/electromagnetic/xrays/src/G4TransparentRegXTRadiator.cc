//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4TransparentRegXTRadiator.cc,v 1.6 2005/05/29 15:46:41 grichine Exp $
// GEANT4 tag $Name: geant4-07-01 $
//

#include <complex>

#include "G4TransparentRegXTRadiator.hh"
#include "Randomize.hh"
#include "G4Integrator.hh"
#include "G4Gamma.hh"

using namespace std;

////////////////////////////////////////////////////////////////////////////
//
// Constructor, destructor

G4TransparentRegXTRadiator::G4TransparentRegXTRadiator(G4LogicalVolume *anEnvelope,
					 G4Material* foilMat,G4Material* gasMat, 
                                         G4double a, G4double b, G4int n,
                                         const G4String& processName) :
  G4VXTRenergyLoss(anEnvelope,foilMat,gasMat,a,b,n,processName)
{
  G4cout<<"Regular transparent X-ray TR  radiator EM process is called"<<G4endl;

  // Build energy and angular integral spectra of X-ray TR photons from
  // a radiator

  fAlphaPlate = 10000;
  fAlphaGas   = 1000;

  BuildTable();
}

///////////////////////////////////////////////////////////////////////////

G4TransparentRegXTRadiator::~G4TransparentRegXTRadiator()
{
  ;
}

///////////////////////////////////////////////////////////////////////////
//
//

G4double G4TransparentRegXTRadiator::SpectralXTRdEdx(G4double energy)
{
  G4double result, sum = 0., tmp, cof1, cof2, cofMin, cofPHC,aMa, bMb, sigma;
  G4int k, kMax, kMin;

  aMa = fPlateThick*GetPlateLinearPhotoAbs(energy);
  bMb = fGasThick*GetGasLinearPhotoAbs(energy);
  sigma = aMa + bMb;
   
  cofPHC  = 4*pi*hbarc;
  tmp     = (fSigma1 - fSigma2)/cofPHC/energy;  
  cof1    = fPlateThick*tmp;
  cof2    = fGasThick*tmp;

  cofMin  =  energy*(fPlateThick + fGasThick)/fGamma/fGamma;
  cofMin += (fPlateThick*fSigma1 + fGasThick*fSigma2)/energy;
  cofMin /= cofPHC;

  //  if (fGamma < 1200) kMin = G4int(cofMin);  // 1200 ?
  // else               kMin = 1;


  kMin = G4int(cofMin);
  if (cofMin > kMin) kMin++;

  // tmp  = (fPlateThick + fGasThick)*energy*fMaxThetaTR;
  // tmp /= cofPHC;
  // kMax = G4int(tmp);
  // if(kMax < 0) kMax = 0;
  // kMax += kMin;
  

  kMax = kMin + 19; //  9; // kMin + G4int(tmp);

  // tmp /= fGamma;
  // if( G4int(tmp) < kMin ) kMin = G4int(tmp);
  // G4cout<<"kMin = "<<kMin<<";    kMax = "<<kMax<<G4endl;

  for( k = kMin; k <= kMax; k++ )
  {
    tmp    = pi*fPlateThick*(k + cof2)/(fPlateThick + fGasThick);
    result = (k - cof1)*(k - cof1)*(k + cof2)*(k + cof2);
    sum   += sin(tmp)*sin(tmp)*abs(k-cofMin)/result;
    //  G4cout<<"k = "<<k<<";    sum = "<<sum<<G4endl;    
  }
  result = 4*( cof1 + cof2 )*( cof1 + cof2 )*sum/energy;
  // result *= ( 1 - exp(-0.5*fPlateNumber*sigma) )/( 1 - exp(-0.5*sigma) );  // fPlateNumber;
  result *= fPlateNumber; // *exp(-0.5*fPlateNumber*sigma); // +1-exp(-0.5*fPlateNumber*sigma); 
  /*  
  fEnergy = energy;
  //  G4Integrator<G4VXTRenergyLoss,G4double(G4VXTRenergyLoss::*)(G4double)> integral;
  G4Integrator<G4TransparentRegXTRadiator,G4double(G4VXTRenergyLoss::*)(G4double)> integral;
 
  tmp = integral.Legendre96(this,&G4VXTRenergyLoss::SpectralAngleXTRdEdx,
                             0.0,0.3*fMaxThetaTR) +
      integral.Legendre96(this,&G4VXTRenergyLoss::SpectralAngleXTRdEdx,
                             0.3*fMaxThetaTR,0.6*fMaxThetaTR) +         
      integral.Legendre96(this,&G4VXTRenergyLoss::SpectralAngleXTRdEdx,
                             0.6*fMaxThetaTR,fMaxThetaTR) ;
  result += tmp;
  */
  return result;
}


///////////////////////////////////////////////////////////////////////////
//
// Approximation for radiator interference factor for the case of
// fully Regular radiator. The plate and gas gap thicknesses are fixed .
// The mean values of the plate and gas gap thicknesses 
// are supposed to be about XTR formation zones but much less than 
// mean absorption length of XTR photons in coresponding material.

G4double 
G4TransparentRegXTRadiator::GetStackFactor( G4double energy, 
                                         G4double gamma, G4double varAngle )
{
  /*
  G4double result, Za, Zb, Ma, Mb, sigma;
  
  Za = GetPlateFormationZone(energy,gamma,varAngle);
  Zb = GetGasFormationZone(energy,gamma,varAngle);
  Ma = GetPlateLinearPhotoAbs(energy);
  Mb = GetGasLinearPhotoAbs(energy);
  sigma = Ma*fPlateThick + Mb*fGasThick;

  G4complex Ca(1.0+0.5*fPlateThick*Ma/fAlphaPlate,fPlateThick/Za/fAlphaPlate); 
  G4complex Cb(1.0+0.5*fGasThick*Mb/fAlphaGas,fGasThick/Zb/fAlphaGas); 

  G4complex Ha = pow(Ca,-fAlphaPlate);  
  G4complex Hb = pow(Cb,-fAlphaGas);
  G4complex H  = Ha*Hb;
  G4complex F1 =   (1.0 - Ha)*(1.0 - Hb )/(1.0 - H)
                 * G4double(fPlateNumber) ;
  G4complex F2 =   (1.0-Ha)*(1.0-Ha)*Hb/(1.0-H)/(1.0-H)
                 * (1.0 - exp(-0.5*fPlateNumber*sigma)) ;
  //    *(1.0 - pow(H,fPlateNumber)) ;
    G4complex R  = (F1 + F2)*OneInterfaceXTRdEdx(energy,gamma,varAngle);
  // G4complex R  = F2*OneInterfaceXTRdEdx(energy,gamma,varAngle);
  result       = 2.0*real(R);  
  return      result;
  */
   // numerically unstable result

  G4double result, Qa, Qb, Q, aZa, bZb, aMa, bMb, D, sigma; 
 
  aZa   = fPlateThick/GetPlateFormationZone(energy,gamma,varAngle);
  bZb   = fGasThick/GetGasFormationZone(energy,gamma,varAngle);
  aMa   = fPlateThick*GetPlateLinearPhotoAbs(energy);
  bMb   = fGasThick*GetGasLinearPhotoAbs(energy);
  sigma = aMa*fPlateThick + bMb*fGasThick;
  Qa    = exp(-0.5*aMa);
  Qb    = exp(-0.5*bMb);
  Q     = Qa*Qb;

  G4complex Ha( Qa*cos(aZa), -Qa*sin(aZa)   );  
  G4complex Hb( Qb*cos(bZb), -Qb*sin(bZb)    );
  G4complex H  = Ha*Hb;
  G4complex Hs = conj(H);
  D            = 1.0 /( (1 - Q)*(1 - Q) + 
                  4*Q*sin(0.5*(aZa + bZb))*sin(0.5*(aZa + bZb)) );
  G4complex F1 = (1.0 - Ha)*(1.0 - Hb)*(1.0 - Hs)
                 * G4double(fPlateNumber)*D;
  G4complex F2 = (1.0 - Ha)*(1.0 - Ha)*Hb*(1.0 - Hs)*(1.0 - Hs)
                   // * (1.0 - pow(H,fPlateNumber)) * D*D;
                 * (1.0 - exp(-0.5*fPlateNumber*sigma)) * D*D;
  G4complex R  = (F1 + F2)*OneInterfaceXTRdEdx(energy,gamma,varAngle);
  result       = 2.0*real(R); 
  return      result;
  
}


//
//
////////////////////////////////////////////////////////////////////////////








