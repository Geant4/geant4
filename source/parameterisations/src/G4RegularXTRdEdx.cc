// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4RegularXTRdEdx.cc,v 1.1 2001-02-27 15:23:48 grichine Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "g4std/complex"

#include "G4RegularXTRdEdx.hh"
#include "Randomize.hh"

#include "G4Gamma.hh"

////////////////////////////////////////////////////////////////////////////
//
// Constructor, destructor

G4RegularXTRdEdx::G4RegularXTRdEdx(G4Envelope *anEnvelope, 
                                                  G4double a, G4double b) :
  G4VXTRdEdx(anEnvelope,a,b)
{
  G4cout<<"Regular X-ray TR dE/dx radiator model is called"<<G4endl ;

  // Build energy and angular integral spectra of X-ray TR photons from
  // a radiator

  BuildTable() ;
}

///////////////////////////////////////////////////////////////////////////

G4RegularXTRdEdx::~G4RegularXTRdEdx()
{
  ;
}



///////////////////////////////////////////////////////////////////////////
//
// Approximation for radiator interference factor for the case of
// fully Regular radiator. The plate and gas gap thicknesses are fixed .
// The mean values of the plate and gas gap thicknesses 
// are supposed to be about XTR formation zones but much less than 
// mean absorption length of XTR photons in coresponding material.

G4double 
G4RegularXTRdEdx::GetStackFactor( G4double energy, 
                                         G4double gamma, G4double varAngle )
{
  G4double result, Qa, Qb, Q, aZa, bZb, aMa, bMb, D, I2 ;
  
  aZa = fPlateThick/GetPlateFormationZone(energy,gamma,varAngle) ;
  bZb = fGasThick/GetGasFormationZone(energy,gamma,varAngle) ;

  aMa = fPlateThick*GetPlateLinearPhotoAbs(energy) ;
  bMb = fGasThick*GetGasLinearPhotoAbs(energy) ;

  Qa = exp(-aMa) ;
  Qb = exp(-bMb) ;
  Q  = Qa*Qb ;

  G4complex Ha( exp(-0.5*aMa)*cos(aZa),
               -exp(-0.5*aMa)*sin(aZa)   ) ; 
 
  G4complex Hb( exp(-0.5*bMb)*cos(bZb),
               -exp(-0.5*bMb)*sin(bZb)    ) ;

  G4complex H  = Ha*Hb ;

  G4complex Hs = G4std::conj(H) ;

  D            = 1.0 /( (1 - sqrt(Q))*(1 - sqrt(Q)) + 
                  4*sqrt(Q)*sin(0.5*(aZa+bZb))*sin(0.5*(aZa+bZb)) ) ;

  G4complex F1 = (1.0 - Ha)*(1.0 - Hb)*(1 - Hs) ;

  F1          *= G4double(fPlateNumber)*D ;

  G4complex F2 = (1.0-Ha)*(1.0-Ha)*Hb*(1.0-Hs)*(1.0-Hs) ;

  F2          *= 1.0 - G4std::pow(H,fPlateNumber) ;

  F2          *= D*D ;

  G4complex R  = (F1 + F2)*OneInterfaceXTRdEdx(energy,gamma,varAngle) ;

  result       = 2.0*G4std::real(R) ;
 
  return      result ;
}


//
//
////////////////////////////////////////////////////////////////////////////








