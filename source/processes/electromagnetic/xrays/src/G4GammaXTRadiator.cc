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
// $Id: G4GammaXTRadiator.cc,v 1.1 2002-01-22 15:22:53 grichine Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "g4std/complex"

#include "G4GammaXTRadiator.hh"
#include "Randomize.hh"

#include "G4Gamma.hh"

////////////////////////////////////////////////////////////////////////////
//
// Constructor, destructor

G4GammaXTRadiator::G4GammaXTRadiator(G4LogicalVolume* anEnvelope, 
				     G4double alphaPlate,
				     G4double alphaGas,
                                     G4Material* foilMat,G4Material* gasMat, 
                                     G4double a, G4double b, G4int n,
                                     const G4String& processName) :
  G4VXTRenergyLoss(anEnvelope,foilMat,gasMat,a,b,n,processName)
{
  G4cout<<"Gammma distributed X-ray TR radiator model is called"<<G4endl ;

  // Build energy and angular integral spectra of X-ray TR photons from
  // a radiator

  fAlphaPlate = alphaPlate ;
  fAlphaGas   = alphaGas   ;
  G4cout<<"fAlphaPlate = "<<fAlphaPlate<<" ; fAlphaGas = "<<fAlphaGas<<G4endl ;

  BuildTable() ;
}

///////////////////////////////////////////////////////////////////////////

G4GammaXTRadiator::~G4GammaXTRadiator()
{
  ;
}



///////////////////////////////////////////////////////////////////////////
//
// Rough approximation for radiator interference factor for the case of
// fully GamDistr radiator. The plate and gas gap thicknesses are distributed 
// according to exponent. The mean values of the plate and gas gap thicknesses 
// are supposed to be about XTR formation zones but much less than 
// mean absorption length of XTR photons in coresponding material.

G4double 
G4GammaXTRadiator::GetStackFactor( G4double energy, 
                                         G4double gamma, G4double varAngle )
{
  G4double result, Za, Zb, Ma, Mb ;
  
  Za = GetPlateFormationZone(energy,gamma,varAngle) ;
  Zb = GetGasFormationZone(energy,gamma,varAngle) ;

  Ma = GetPlateLinearPhotoAbs(energy) ;
  Mb = GetGasLinearPhotoAbs(energy) ;


  G4complex Ca(1.0+0.5*fPlateThick*Ma/fAlphaPlate,fPlateThick/Za/fAlphaPlate) ; 
  G4complex Cb(1.0+0.5*fGasThick*Mb/fAlphaGas,fGasThick/Zb/fAlphaGas) ; 

  G4complex Ha = G4std::pow(Ca,-fAlphaPlate) ;  
  G4complex Hb = G4std::pow(Cb,-fAlphaGas) ;
  G4complex H  = Ha*Hb ;

  G4complex F1 =   (1.0 - Ha)*(1.0 - Hb )/(1.0 - H)
                 * G4double(fPlateNumber) ;

  G4complex F2 =   (1.0-Ha)*(1.0-Ha)*Hb/(1.0-H)/(1.0-H)
                 * (1.0 - G4std::pow(H,fPlateNumber)) ;

  G4complex R  = (F1 + F2)*OneInterfaceXTRdEdx(energy,gamma,varAngle) ;

  result       = 2.0*G4std::real(R) ;
  
  return      result ;
}


//
//
////////////////////////////////////////////////////////////////////////////








