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
// $Id: G4RegularXTRadiator.cc,v 1.2 2002-01-18 17:26:21 grichine Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "g4std/complex"

#include "G4RegularXTRadiator.hh"
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

  BuildTable() ;
}

///////////////////////////////////////////////////////////////////////////

G4RegularXTRadiator::~G4RegularXTRadiator()
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
G4RegularXTRadiator::GetStackFactor( G4double energy, 
                                         G4double gamma, G4double varAngle )
{
  G4double result, Qa, Qb, Q, aZa, bZb, aMa, bMb, D ;
  
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

  G4complex F1 = (1.0 - Ha)*(1.0 - Hb)*(1.0 - Hs)
                 * G4double(fPlateNumber)*D ;

  G4complex F2 = (1.0-Ha)*(1.0-Ha)*Hb*(1.0-Hs)*(1.0-Hs)
                 * (1.0 - G4std::pow(H,fPlateNumber)) * D*D ;

  G4complex R  = (F1 + F2)*OneInterfaceXTRdEdx(energy,gamma,varAngle) ;

  result       = 2.0*G4std::real(R) ;
 
  return      result ;
}


//
//
////////////////////////////////////////////////////////////////////////////








