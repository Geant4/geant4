// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4FoamXTRdEdx.cc,v 1.1 2001-02-27 15:23:48 grichine Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "g4std/complex"

#include "G4FoamXTRdEdx.hh"
#include "Randomize.hh"

#include "G4Gamma.hh"

////////////////////////////////////////////////////////////////////////////
//
// Constructor, destructor

G4FoamXTRdEdx::G4FoamXTRdEdx(G4Envelope *anEnvelope, 
                                                  G4double a, G4double b) :
  G4VXTRdEdx(anEnvelope,a,b)
{
  G4cout<<"Foam XTR dE/dx model is called"<<G4endl ;

  // Build energy and angular integral spectra of X-ray TR photons inside
  // a radiator

    BuildTable() ;
}

///////////////////////////////////////////////////////////////////////////

G4FoamXTRdEdx::~G4FoamXTRdEdx()
{
  ;
}



///////////////////////////////////////////////////////////////////////////
//
// Rough approximation for radiator interference factor for the case of
// fully Foam radiator. The plate and gas gap thicknesses are distributed 
// according to exponent. The mean values of the plate and gas gap thicknesses 
// are supposed to be about XTR formation zones but much less than 
// mean absorption length of XTR photons in coresponding material.

G4double 
G4FoamXTRdEdx::GetStackFactor( G4double energy, G4double gamma, G4double varAngle )
{
  G4double result, Za, Zb, Ma, Mb ;
  
  Za = GetPlateFormationZone(energy,gamma,varAngle) ;
  Zb = GetGasFormationZone(energy,gamma,varAngle) ;

  Ma = GetPlateLinearPhotoAbs(energy) ;
  Mb = GetGasLinearPhotoAbs(energy) ;

  G4complex Ca(1.0+0.5*fPlateThick*Ma,fPlateThick/Za) ; 
  G4complex Cb(1.0+0.5*fGasThick*Mb,fGasThick/Zb) ; 

  G4complex Ha = 1.0/Ca ;  
  G4complex Hb = 1.0/Cb ;
  G4complex H  = Ha*Hb ;

  G4complex F1 = (1.0-Ha)*(1.0-Hb)/(1.0-H) ;
  
  F1          *=  G4double(fPlateNumber) ;

  G4complex F2 = (1.0-Ha)*(1.0-Ha)*Hb/(1.0-H)/(1.0-H) ;

  F2          *= 1.0 - G4std::pow(H,fPlateNumber) ;

  G4complex R  = (F1 + F2)*OneInterfaceXTRdEdx(energy,gamma,varAngle) ;

  result       = 2.0*G4std::real(R) ;

  return      result ;
}


//
//
////////////////////////////////////////////////////////////////////////////








