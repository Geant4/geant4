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
// $Id: G4PlateIrrGasXTRdEdx.cc,v 1.3 2001-09-18 09:02:03 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "g4std/complex"

#include "G4PlateIrrGasXTRdEdx.hh"
#include "Randomize.hh"

#include "G4Gamma.hh"

////////////////////////////////////////////////////////////////////////////
//
// Constructor, destructor

G4PlateIrrGasXTRdEdx::G4PlateIrrGasXTRdEdx(G4Envelope *anEnvelope, 
                                                  G4double a, G4double b) :
  G4VXTRdEdx(anEnvelope,a,b)
{
  G4cout<<"PlateIrrGas X-ray TR radiator model is called"<<G4endl ;

  // Build energy and angular integral spectra of X-ray TR photons from
  // a radiator

  BuildTable() ;
}

///////////////////////////////////////////////////////////////////////////

G4PlateIrrGasXTRdEdx::~G4PlateIrrGasXTRdEdx()
{
  ;
}



///////////////////////////////////////////////////////////////////////////
//
// Approximation for radiator interference factor for the case of
// fully PlateIrrGas radiator. The plate thickness is fixed .
// The mean values of gas gap thicknesses 
// are supposed to be about XTR formation zones but much less than 
// mean absorption length of XTR photons in coresponding material.

G4double 
G4PlateIrrGasXTRdEdx::GetStackFactor( G4double energy, 
                                         G4double gamma, G4double varAngle )
{
  G4double result, Za, Ma ;
  
  Za = GetPlateFormationZone(energy,gamma,varAngle) ;
  Ma = GetPlateLinearPhotoAbs(energy) ;



  

  // G4complex Ca(1.0+0.5*fPlateThick*Ma,fPlateThick/Za) ; 

  G4complex Ha( exp(-0.5*fPlateThick*Ma)*cos(fPlateThick/Za),
               -exp(-0.5*fPlateThick*Ma)*sin(fPlateThick/Za)   ) ; 
 

  G4complex R = (1.0 - Ha) ;

  R          *= G4double(fPlateNumber)*OneInterfaceXTRdEdx(energy,gamma,varAngle) ;

  result       = 2.0*G4std::real(R)  ;


  return      result ;
}


//
//
////////////////////////////////////////////////////////////////////////////
