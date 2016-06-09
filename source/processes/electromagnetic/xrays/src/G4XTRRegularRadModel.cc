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

#include <complex>

#include "G4XTRRegularRadModel.hh"
#include "Randomize.hh"

#include "G4Gamma.hh"
using namespace std;

////////////////////////////////////////////////////////////////////////////
//
// Constructor, destructor

G4XTRRegularRadModel::G4XTRRegularRadModel(G4LogicalVolume *anEnvelope,
                                         G4Material* foilMat,G4Material* gasMat,
                                         G4double a, G4double b, G4int n,
                                         const G4String& processName) :
  G4VXTRenergyLoss(anEnvelope,foilMat,gasMat,a,b,n,processName)  
{
  G4cout<<" XTR Regular discrete radiator model is called"<<G4endl ;

  fExitFlux = true;

  // Build energy and angular integral spectra of X-ray TR photons from
  // a radiator

  // BuildTable() ;
}

///////////////////////////////////////////////////////////////////////////

G4XTRRegularRadModel::~G4XTRRegularRadModel()
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
G4XTRRegularRadModel::GetStackFactor( G4double energy, 
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

  I2           = 1.; // 2.0*std::real(F2) ;

  I2           /= (1 - std::sqrt(Q))*(1 - std::sqrt(Q)) + 
                  4*std::sqrt(Q)*std::sin(0.5*(aZa+bZb))*std::sin(0.5*(aZa+bZb)) ;

  I2           /= Q*( (std::sqrt(Q)-std::cos(aZa+bZb))*(std::sqrt(Q)-std::cos(aZa+bZb)) + 
                      std::sin(aZa+bZb)*std::sin(aZa+bZb)   ) ;

  G4complex stack  = 2.*I2*F2;
            stack += result;
            stack *= OneInterfaceXTRdEdx(energy,gamma,varAngle);

	    // result       += I2 ;
  result = std::real(stack);

  return      result ;
}


//
//
////////////////////////////////////////////////////////////////////////////








