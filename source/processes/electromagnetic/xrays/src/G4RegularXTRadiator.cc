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
// $Id: G4RegularXTRadiator.cc 68037 2013-03-13 14:15:08Z gcosmo $
//

#include <complex>

#include "G4RegularXTRadiator.hh"
#include "G4PhysicalConstants.hh"
#include "Randomize.hh"

#include "G4Gamma.hh"

////////////////////////////////////////////////////////////////////////////
//
// Constructor, destructor

G4RegularXTRadiator::G4RegularXTRadiator(G4LogicalVolume *anEnvelope,
					 G4Material* foilMat,G4Material* gasMat, 
                                         G4double a, G4double b, G4int n,
                                         const G4String& processName) :
  G4VXTRenergyLoss(anEnvelope,foilMat,gasMat,a,b,n,processName)
{
  G4cout<<"Regular X-ray TR  radiator EM process is called"<<G4endl ;

  // Build energy and angular integral spectra of X-ray TR photons from
  // a radiator

  fAlphaPlate = 10000;
  fAlphaGas   = 1000;
  G4cout<<"fAlphaPlate = "<<fAlphaPlate<<" ; fAlphaGas = "<<fAlphaGas<<G4endl ;

  // BuildTable() ;
}

///////////////////////////////////////////////////////////////////////////

G4RegularXTRadiator::~G4RegularXTRadiator()
{
  ;
}

///////////////////////////////////////////////////////////////////////////
//
//

G4double G4RegularXTRadiator::SpectralXTRdEdx(G4double energy)
{
  G4double result, sum = 0., tmp, cof1, cof2, cofMin, cofPHC, theta2, theta2k; 
    G4double aMa, bMb ,sigma, dump;
  G4int k, kMax, kMin;

  aMa = fPlateThick*GetPlateLinearPhotoAbs(energy);
  bMb = fGasThick*GetGasLinearPhotoAbs(energy);
  sigma = 0.5*(aMa + bMb);
  dump = std::exp(-fPlateNumber*sigma);
  if(verboseLevel > 2)  G4cout<<" dump = "<<dump<<G4endl;  
  cofPHC  = 4*pi*hbarc;
  tmp     = (fSigma1 - fSigma2)/cofPHC/energy;  
  cof1    = fPlateThick*tmp;
  cof2    = fGasThick*tmp;

  cofMin  =  energy*(fPlateThick + fGasThick)/fGamma/fGamma;
  cofMin += (fPlateThick*fSigma1 + fGasThick*fSigma2)/energy;
  cofMin /= cofPHC;

  theta2 = cofPHC/(energy*(fPlateThick + fGasThick));

  //  if (fGamma < 1200) kMin = G4int(cofMin);  // 1200 ?
  // else               kMin = 1;


  kMin = G4int(cofMin);
  if (cofMin > kMin) kMin++;

  // tmp  = (fPlateThick + fGasThick)*energy*fMaxThetaTR;
  // tmp /= cofPHC;
  // kMax = G4int(tmp);
  // if(kMax < 0) kMax = 0;
  // kMax += kMin;
  

  kMax = kMin + 49; //  19; // kMin + G4int(tmp);

  // tmp /= fGamma;
  // if( G4int(tmp) < kMin ) kMin = G4int(tmp);

  if(verboseLevel > 2)
  {    
    G4cout<<cof1<<"     "<<cof2<<"        "<<cofMin<<G4endl;
    G4cout<<"kMin = "<<kMin<<";    kMax = "<<kMax<<G4endl;
  }
  for( k = kMin; k <= kMax; k++ )
  {
    tmp    = pi*fPlateThick*(k + cof2)/(fPlateThick + fGasThick);
    result = (k - cof1)*(k - cof1)*(k + cof2)*(k + cof2);
    // tmp = std::sin(tmp)*std::sin(tmp)*std::abs(k-cofMin)/result;
    if( k == kMin && kMin == G4int(cofMin) )
    {
      sum   += 0.5*std::sin(tmp)*std::sin(tmp)*std::abs(k-cofMin)/result;
    }
    else
    {
      sum   += std::sin(tmp)*std::sin(tmp)*std::abs(k-cofMin)/result;
    }
    theta2k = std::sqrt(theta2*std::abs(k-cofMin));

    if(verboseLevel > 2)
    {    
      // G4cout<<"k = "<<k<<"; sqrt(theta2k) = "<<theta2k<<"; tmp = "<<std::sin(tmp)*std::sin(tmp)*std::abs(k-cofMin)/result
      //     <<";    sum = "<<sum<<G4endl;  
      G4cout<<k<<"   "<<theta2k<<"     "<<std::sin(tmp)*std::sin(tmp)*std::abs(k-cofMin)/result
              <<"      "<<sum<<G4endl;  
    }  
  }
  result = 2*( cof1 + cof2 )*( cof1 + cof2 )*sum/energy;
  // result *= ( 1 - std::exp(-0.5*fPlateNumber*sigma) )/( 1 - std::exp(-0.5*sigma) );  
  // fPlateNumber;
  result *= ( 1 - dump + 2*dump*fPlateNumber ); 
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
G4RegularXTRadiator::GetStackFactor( G4double energy, 
                                         G4double gamma, G4double varAngle )
{

  // some gamma (10000/1000) like algorithm

  G4double result, Za, Zb, Ma, Mb;
  
  Za = GetPlateFormationZone(energy,gamma,varAngle);
  Zb = GetGasFormationZone(energy,gamma,varAngle);

  Ma = GetPlateLinearPhotoAbs(energy);
  Mb = GetGasLinearPhotoAbs(energy);


  G4complex Ca(1.0+0.5*fPlateThick*Ma/fAlphaPlate,fPlateThick/Za/fAlphaPlate); 
  G4complex Cb(1.0+0.5*fGasThick*Mb/fAlphaGas,fGasThick/Zb/fAlphaGas); 

  G4complex Ha = std::pow(Ca,-fAlphaPlate);  
  G4complex Hb = std::pow(Cb,-fAlphaGas);
  G4complex H  = Ha*Hb;
  
  G4complex F1 =   (1.0 - Ha)*(1.0 - Hb )/(1.0 - H)
                 * G4double(fPlateNumber);

  G4complex F2 =   (1.0-Ha)*(1.0-Ha)*Hb/(1.0-H)/(1.0-H)
                 * (1.0 - std::pow(H,fPlateNumber));

  G4complex R  = (F1 + F2)*OneInterfaceXTRdEdx(energy,gamma,varAngle);
  
  result       = 2.0*std::real(R);
  
  return      result;
  
  /*
   // numerically stable but slow algorithm

  G4double result, Qa, Qb, Q, aZa, bZb, aMa, bMb;   // , D; 
 
  aZa = fPlateThick/GetPlateFormationZone(energy,gamma,varAngle);
  bZb = fGasThick/GetGasFormationZone(energy,gamma,varAngle);
  aMa = fPlateThick*GetPlateLinearPhotoAbs(energy);
  bMb = fGasThick*GetGasLinearPhotoAbs(energy);
  Qa = std::exp(-aMa);
  Qb = std::exp(-bMb);
  Q  = Qa*Qb;
  G4complex Ha( std::exp(-0.5*aMa)*std::cos(aZa),
               -std::exp(-0.5*aMa)*std::sin(aZa)   );  
  G4complex Hb( std::exp(-0.5*bMb)*std::cos(bZb),
               -std::exp(-0.5*bMb)*std::sin(bZb)    );
  G4complex H  = Ha*Hb;
  
  G4complex Hs = conj(H);
  D            = 1.0 /( (1 - std::sqrt(Q))*(1 - std::sqrt(Q)) + 
                  4*std::sqrt(Q)*std::sin(0.5*(aZa+bZb))*std::sin(0.5*(aZa+bZb)) );
  G4complex F1 = (1.0 - Ha)*(1.0 - Hb)*(1.0 - Hs)
                 * G4double(fPlateNumber)*D;
  G4complex F2 = (1.0-Ha)*(1.0-Ha)*Hb*(1.0-Hs)*(1.0-Hs)
                 * (1.0 - std::pow(H,fPlateNumber)) * D*D;
  G4complex R  = (F1 + F2)*OneInterfaceXTRdEdx(energy,gamma,varAngle);
  

  G4complex S(0.,0.), c(1.,0.);
  G4int k;
  for(k = 1; k < fPlateNumber; k++)
  {
    c *= H;
    S += ( G4double(fPlateNumber) - G4double(k) )*c; 
  }
  G4complex R  = (2.- Ha - 1./Ha)*S + (1. - Ha)*G4double(fPlateNumber);
            R *= OneInterfaceXTRdEdx(energy,gamma,varAngle);
  result       = 2.0*std::real(R); 
  return      result;
  */
}


//
//
////////////////////////////////////////////////////////////////////////////








