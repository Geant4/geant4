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
// $Id: G4PlateIrrGasXrayTRmodel.cc,v 1.2 2004/12/07 09:00:05 gcosmo Exp $
// GEANT4 tag $Name: geant4-07-00-cand-03 $
//

#include <complex>

#include "G4PlateIrrGasXrayTRmodel.hh"
#include "Randomize.hh"

#include "G4Gamma.hh"

////////////////////////////////////////////////////////////////////////////
//
// Constructor, destructor

G4PlateIrrGasXrayTRmodel::G4PlateIrrGasXrayTRmodel(G4Envelope *anEnvelope, 
                                                  G4double a, G4double b) :
  G4VXrayTRadModel(anEnvelope,a,b)
{
  G4cout<<"PlateIrrGas X-ray TR radiator model is called"<<G4endl ;

  // Build energy and angular integral spectra of X-ray TR photons from
  // a radiator

  BuildTable() ;
}

///////////////////////////////////////////////////////////////////////////

G4PlateIrrGasXrayTRmodel::~G4PlateIrrGasXrayTRmodel()
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
G4PlateIrrGasXrayTRmodel::GetStackFactor( G4double energy, 
                                         G4double gamma, G4double varAngle )
{
  G4double result, Qa, Qb, Q, Za, Ma, Mb ;
  
  Za = GetPlateFormationZone(energy,gamma,varAngle) ;
  //  Zb = GetGasFormationZone(energy,gamma,varAngle) ;

  Ma = GetPlateLinearPhotoAbs(energy) ;
  Mb = GetGasLinearPhotoAbs(energy) ;

  Qa = std::exp(-fPlateThick*Ma) ;
  Qb = std::exp(-fGasThick*Mb) ;
  Q  = Qa*Qb ;

  /* *****************************************************

  //  G4complex Ca(1.0+0.5*fPlateThick*Ma,fPlateThick/Za) ; 
  //  G4complex Cb(1.0+0.5*fGasThick*Mb,fGasThick/Zb) ; 

  G4complex Ha( std::exp(-0.5*fPlateThick*Ma)*std::cos(fPlateThick/Za),
              -std::exp(-0.5*fPlateThick*Ma)*std::sin(fPlateThick/Za)   ) ; 
 
  G4complex Hb( std::exp(-0.5*fGasThick*Mb)*std::cos(fGasThick/Zb),
               -std::exp(-0.5*fGasThick*Mb)*std::sin(fGasThick/Za)    ) ;
  G4complex H  = Ha*Hb ;

  G4complex F1 = ( 0.5*(1+Qa)*(1+H) - Ha - Qa*Hb )/(1-H) ;

  G4complex F2 = (1-Ha)*(Qa-Ha)*Hb/(1-H)/(Q-H) ;

  F2          *= std::pow(Q,G4double(fPlateNumber)) - std::pow(H,fPlateNumber) ;

  result      = ( 1 - std::pow(Q,G4double(fPlateNumber)) )/( 1 - Q ) ;

  result     *= 2.0*std::real(F1) ;

  result     += 2.0*std::real(F2) ;

  ***************************************************************** */

  result      = ( 1 - std::pow(Q,G4double(fPlateNumber)) )/( 1 - Q ) ;

  result *= 1 + Qa -2*std::sqrt(Qa)*std::cos(fPlateThick/Za) ;

  return      result ;
}


//
//
////////////////////////////////////////////////////////////////////////////








