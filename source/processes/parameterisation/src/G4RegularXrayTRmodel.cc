// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4RegularXrayTRmodel.cc,v 1.2 2000-06-06 08:47:48 mora Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "g4std/complex"

#include "G4RegularXrayTRmodel.hh"
#include "Randomize.hh"

#include "G4Gamma.hh"

typedef G4std::complex<double> G4complex ;


////////////////////////////////////////////////////////////////////////////
//
// Constructor, destructor

G4RegularXrayTRmodel::G4RegularXrayTRmodel(G4Envelope *anEnvelope, 
                                                  G4double a, G4double b) :
  G4VXrayTRmodel(anEnvelope,a,b)
{
  G4cout<<"Regular X-ray TR radiator model is called"<<G4endl ;

  // Build energy and angular integral spectra of X-ray TR photons from
  // a radiator

  BuildTable() ;
}

///////////////////////////////////////////////////////////////////////////

G4RegularXrayTRmodel::~G4RegularXrayTRmodel()
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
G4RegularXrayTRmodel::GetStackFactor( G4double energy, 
                                         G4double gamma, G4double varAngle )
{
  G4double result, Qa, Qb, Q, aZa, bZb, aMa, bMb, I2 ;
  
  aZa = fPlateThick/GetPlateFormationZone(energy,gamma,varAngle) ;
  bZb = fGasThick/GetGasFormationZone(energy,gamma,varAngle) ;

  aMa = fPlateThick*GetPlateLinearPhotoAbs(energy) ;
  bMb = fGasThick*GetGasLinearPhotoAbs(energy) ;

  Qa = exp(-aMa) ;
  Qb = exp(-bMb) ;
  Q  = Qa*Qb ;

  //  G4complex Ca(1.0+0.5*fPlateThick*Ma,fPlateThick/Za) ; 
  //  G4complex Cb(1.0+0.5*fGasThick*Mb,fGasThick/Zb) ; 

  G4complex Ha( exp(-0.5*aMa)*cos(aZa),
               -exp(-0.5*aMa)*sin(aZa)   ) ; 
 
  G4complex Hb( exp(-0.5*bMb)*cos(bZb),
               -exp(-0.5*bMb)*sin(bZb)    ) ;

  G4complex H  = Ha*Hb ;

  G4complex Hs = conj(H) ;

  //  G4complex F1 = ( 0.5*(1+Qa)*(1+H) - Ha - Qa*Hb )/(1-H) ;

  G4complex F2 = (1.0-Ha)*(Qa-Ha)*Hb*(1.0-Hs)*(Q-Hs) ;

  F2          *= pow(Q,fPlateNumber) - pow(H,fPlateNumber) ;

  result       = ( 1 - pow(Q,fPlateNumber) )/( 1 - Q ) ;

  result      *= (1 - Qa)*(1 + Qa - 2*sqrt(Qa)*cos(aZa)) ;

  result      /= (1 - sqrt(Q))*(1 - sqrt(Q)) + 
                  4*sqrt(Q)*sin(0.5*(aZa+bZb))*sin(0.5*(aZa+bZb)) ;

  I2           = 2.0*F2.real() ;

  I2           /= (1 - sqrt(Q))*(1 - sqrt(Q)) + 
                  4*sqrt(Q)*sin(0.5*(aZa+bZb))*sin(0.5*(aZa+bZb)) ;

  I2           /= Q*( (sqrt(Q)-cos(aZa+bZb))*(sqrt(Q)-cos(aZa+bZb)) + 
                      sin(aZa+bZb)*sin(aZa+bZb)   ) ;

  result       += I2 ;

  return      result ;
}


//
//
////////////////////////////////////////////////////////////////////////////








