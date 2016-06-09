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
// $Id: G4RegularXrayTRmodel.cc,v 1.2 2004/12/07 09:00:05 gcosmo Exp $
// GEANT4 tag $Name: geant4-07-00-cand-03 $
//

#include <complex>

#include "G4RegularXrayTRmodel.hh"
#include "Randomize.hh"

#include "G4Gamma.hh"

////////////////////////////////////////////////////////////////////////////
//
// Constructor, destructor

G4RegularXrayTRmodel::G4RegularXrayTRmodel(G4Envelope *anEnvelope, 
                                                  G4double a, G4double b) :
  G4VXrayTRadModel(anEnvelope,a,b)
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

  Qa = std::exp(-aMa) ;
  Qb = std::exp(-bMb) ;
  Q  = Qa*Qb ;

  //  G4complex Ca(1.0+0.5*fPlateThick*Ma,fPlateThick/Za) ; 
  //  G4complex Cb(1.0+0.5*fGasThick*Mb,fGasThick/Zb) ; 

  G4complex Ha( std::exp(-0.5*aMa)*std::cos(aZa),
               -std::exp(-0.5*aMa)*std::sin(aZa)   ) ; 
 
  G4complex Hb( std::exp(-0.5*bMb)*std::cos(bZb),
               -std::exp(-0.5*bMb)*std::sin(bZb)    ) ;

  G4complex H  = Ha*Hb ;

  G4complex Hs = std::conj(H) ;

  //  G4complex F1 = ( 0.5*(1+Qa)*(1+H) - Ha - Qa*Hb )/(1-H) ;

  G4complex F2 = (1.0-Ha)*(Qa-Ha)*Hb*(1.0-Hs)*(Q-Hs) ;

  F2          *= std::pow(Q,G4double(fPlateNumber)) - std::pow(H,fPlateNumber) ;

  result       = ( 1 - std::pow(Q,G4double(fPlateNumber)) )/( 1 - Q ) ;

  result      *= (1 - Qa)*(1 + Qa - 2*std::sqrt(Qa)*std::cos(aZa)) ;

  result      /= (1 - std::sqrt(Q))*(1 - std::sqrt(Q)) + 
                  4*std::sqrt(Q)*std::sin(0.5*(aZa+bZb))*std::sin(0.5*(aZa+bZb)) ;

  I2           = 2.0*std::real(F2) ;

  I2           /= (1 - std::sqrt(Q))*(1 - std::sqrt(Q)) + 
                  4*std::sqrt(Q)*std::sin(0.5*(aZa+bZb))*std::sin(0.5*(aZa+bZb)) ;

  I2           /= Q*( (std::sqrt(Q)-std::cos(aZa+bZb))*(std::sqrt(Q)-std::cos(aZa+bZb)) + 
                      std::sin(aZa+bZb)*std::sin(aZa+bZb)   ) ;

  result       += I2 ;

  return      result ;
}


//
//
////////////////////////////////////////////////////////////////////////////








