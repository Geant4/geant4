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
/// \file electromagnetic/TestEm10/src/XTRTransparentRegRadModel.cc
/// \brief Implementation of the XTRTransparentRegRadModel class
//
//

#include <complex>

#include "XTRTransparentRegRadModel.hh"
#include "Randomize.hh"
#include "G4Integrator.hh"
#include "G4Gamma.hh"
#include "G4PhysicalConstants.hh"

using namespace std;

////////////////////////////////////////////////////////////////////////////
//
// Constructor, destructor

XTRTransparentRegRadModel::XTRTransparentRegRadModel(G4LogicalVolume *anEnvelope,
                                         G4Material* foilMat,G4Material* gasMat, 
                                         G4double a, G4double b, G4int n,
                                         const G4String& processName) :
  G4VXTRenergyLoss(anEnvelope,foilMat,gasMat,a,b,n,processName)
{
  G4cout<<"Regular transparent X-ray TR  radiator EM process is called"<<G4endl;

  // Build energy and angular integral spectra of X-ray TR photons from
  // a radiator
  fExitFlux   = true;
  fAlphaPlate = 10000;
  fAlphaGas   = 1000;

  //  BuildTable();
}

///////////////////////////////////////////////////////////////////////////

XTRTransparentRegRadModel::~XTRTransparentRegRadModel()
{
  ;
}

///////////////////////////////////////////////////////////////////////////
//
//

G4double XTRTransparentRegRadModel::SpectralXTRdEdx(G4double energy)
{
  G4double result, sum = 0., tmp, cof1, cof2, cofMin, cofPHC,aMa, bMb, sigma;
  G4int k, kMax, kMin;

  aMa = GetPlateLinearPhotoAbs(energy);
  bMb = GetGasLinearPhotoAbs(energy);

  // if(fCompton)
  {
    aMa += GetPlateCompton(energy);
    bMb += GetGasCompton(energy);
  }
  aMa *= fPlateThick;
  bMb *= fGasThick;

  sigma = aMa + bMb;
   
  cofPHC  = 4*pi*hbarc;
  cofPHC *= 200./197.;
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
  

  kMax = kMin + 9; // 5; // 9; //   kMin + G4int(tmp);

  // tmp /= fGamma;
  // if( G4int(tmp) < kMin ) kMin = G4int(tmp);
  // G4cout<<"kMin = "<<kMin<<";    kMax = "<<kMax<<G4endl;

  for( k = kMin; k <= kMax; k++ )
  {
    tmp    = pi*fPlateThick*(k + cof2)/(fPlateThick + fGasThick);
    result = (k - cof1)*(k - cof1)*(k + cof2)*(k + cof2);

    if( k == kMin && kMin == G4int(cofMin) )
    {
      sum   += 0.5*sin(tmp)*sin(tmp)*std::abs(k-cofMin)/result;
    }
    else
    {
      sum   += sin(tmp)*sin(tmp)*std::abs(k-cofMin)/result;
    }
    //  G4cout<<"k = "<<k<<";    sum = "<<sum<<G4endl;    
  }
  result = 4.*( cof1 + cof2 )*( cof1 + cof2 )*sum/energy;
  result *= ( 1. - exp(-fPlateNumber*sigma) )/( 1. - exp(-sigma) );  
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
XTRTransparentRegRadModel::GetStackFactor( G4double energy, 
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








