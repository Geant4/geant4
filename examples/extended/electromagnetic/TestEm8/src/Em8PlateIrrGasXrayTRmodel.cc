// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Em8PlateIrrGasXrayTRmodel.cc,v 1.1 2000-03-03 09:16:15 grichine Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "g4std/complex"

#include "Em8PlateIrrGasXrayTRmodel.hh"
#include "Randomize.hh"

#include "G4Gamma.hh"

typedef G4std::complex<double> G4complex ;


////////////////////////////////////////////////////////////////////////////
//
// Constructor, destructor

Em8PlateIrrGasXrayTRmodel::Em8PlateIrrGasXrayTRmodel(G4Envelope *anEnvelope, 
                                                  G4double a, G4double b) :
  Em8XrayTRmodel(anEnvelope,a,b)
{
  G4cout<<"PlateIrrGas X-ray TR radiator model is called"<<G4endl ;

  // Build energy and angular integral spectra of X-ray TR photons from
  // a radiator

  BuildTable() ;
}

///////////////////////////////////////////////////////////////////////////

Em8PlateIrrGasXrayTRmodel::~Em8PlateIrrGasXrayTRmodel()
{
  ;
}



///////////////////////////////////////////////////////////////////////////
//
// Approximation for radiator interference factor for the case of
// fully PlateIrrGas radiator. The plate and gas gap thicknesses are fixed .
// The mean values of the plate and gas gap thicknesses 
// are supposed to be about XTR formation zones but much less than 
// mean absorption length of XTR photons in coresponding material.

G4double 
Em8PlateIrrGasXrayTRmodel::GetStackFactor( G4double energy, 
                                         G4double gamma, G4double varAngle )
{
  G4double result, Qa, Qb, Q, Za, Zb, Ma, Mb ;
  
  Za = GetPlateFormationZone(energy,gamma,varAngle) ;
  //  Zb = GetGasFormationZone(energy,gamma,varAngle) ;

  Ma = GetPlateLinearPhotoAbs(energy) ;
  Mb = GetGasLinearPhotoAbs(energy) ;

  Qa = exp(-fPlateThick*Ma) ;
  Qb = exp(-fGasThick*Mb) ;
  Q  = Qa*Qb ;

  /* *****************************************************

  //  G4complex Ca(1.0+0.5*fPlateThick*Ma,fPlateThick/Za) ; 
  //  G4complex Cb(1.0+0.5*fGasThick*Mb,fGasThick/Zb) ; 

  G4complex Ha( exp(-0.5*fPlateThick*Ma)*cos(fPlateThick/Za),
              -exp(-0.5*fPlateThick*Ma)*sin(fPlateThick/Za)   ) ; 
 
  G4complex Hb( exp(-0.5*fGasThick*Mb)*cos(fGasThick/Zb),
               -exp(-0.5*fGasThick*Mb)*sin(fGasThick/Za)    ) ;
  G4complex H  = Ha*Hb ;

  G4complex F1 = ( 0.5*(1+Qa)*(1+H) - Ha - Qa*Hb )/(1-H) ;

  G4complex F2 = (1-Ha)*(Qa-Ha)*Hb/(1-H)/(Q-H) ;

  F2          *= pow(Q,fPlateNumber) - pow(H,fPlateNumber) ;

  result      = ( 1 - pow(Q,fPlateNumber) )/( 1 - Q ) ;

  result     *= 2.0*F1.real() ;

  result     += 2.0*F2.real() ;

  ***************************************************************** */

  result      = ( 1 - pow(Q,fPlateNumber) )/( 1 - Q ) ;

  result *= 1 + Qa -2*cos(fPlateThick/Za) ;

  return      result ;
}


//
//
////////////////////////////////////////////////////////////////////////////








