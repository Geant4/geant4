// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4GamDistrXrayTRmodel.cc,v 1.3 2000-06-15 17:38:39 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "g4std/complex"

#include "G4GamDistrXrayTRmodel.hh"
#include "Randomize.hh"

#include "G4Gamma.hh"

////////////////////////////////////////////////////////////////////////////
//
// Constructor, destructor

G4GamDistrXrayTRmodel::G4GamDistrXrayTRmodel(G4Envelope *anEnvelope, 
					       G4double a, G4double alphaPlate,
                                               G4double b, G4double alphaGas) :
  G4VXrayTRmodel(anEnvelope,a,b)
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

G4GamDistrXrayTRmodel::~G4GamDistrXrayTRmodel()
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
G4GamDistrXrayTRmodel::GetStackFactor( G4double energy, 
                                         G4double gamma, G4double varAngle )
{
  G4double result, Qa, Qb, Q, Za, Zb, Ma, Mb ;
  
  Za = GetPlateFormationZone(energy,gamma,varAngle) ;
  Zb = GetGasFormationZone(energy,gamma,varAngle) ;

  Ma = GetPlateLinearPhotoAbs(energy) ;
  Mb = GetGasLinearPhotoAbs(energy) ;

  Qa = ( 1.0 + fPlateThick*Ma/fAlphaPlate ) ;
  Qa = pow(Qa,-fAlphaPlate) ;
  Qb = ( 1.0 + fGasThick*Mb/fAlphaGas ) ;
  Qb = pow(Qb,-fAlphaGas) ;
  Q  = Qa*Qb ;

  G4complex Ca(1.0+0.5*fPlateThick*Ma/fAlphaPlate,fPlateThick/Za/fAlphaPlate) ; 
  G4complex Cb(1.0+0.5*fGasThick*Mb/fAlphaGas,fGasThick/Zb/fAlphaGas) ; 

  G4complex Ha = G4std::pow(Ca,-fAlphaPlate) ;  
  G4complex Hb = G4std::pow(Cb,-fAlphaGas) ;
  G4complex H  = Ha*Hb ;

  G4complex F1 = ( 0.5*(1+Qa)*(1.0+H) - Ha - Qa*Hb )/(1.0-H) ;

  G4complex F2 = (1.0-Ha)*(Qa-Ha)*Hb/(1.0-H)/(Q-H) ;

  F2          *= pow(Q,G4double(fPlateNumber)) - G4std::pow(H,fPlateNumber) ;

  result      = ( 1 - pow(Q,G4double(fPlateNumber)) )/( 1 - Q ) ;

  result     *= 2.0*G4std::real(F1) ;

  result     += 2.0*G4std::real(F2) ;

  return      result ;
}


//
//
////////////////////////////////////////////////////////////////////////////








