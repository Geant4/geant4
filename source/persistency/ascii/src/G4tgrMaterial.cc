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
// $Id: G4tgrMaterial.cc 66363 2012-12-18 09:12:54Z gcosmo $
//
//
// class G4tgrMaterial

// History:
// - Created.                                 P.Arce, CIEMAT (November 2007)
// -------------------------------------------------------------------------

#include "G4tgrMaterial.hh"
#include "G4PhysicalConstants.hh"

// -------------------------------------------------------------------------
G4tgrMaterial::G4tgrMaterial()
   : theName("Material"), theDensity(0.), theNoComponents(0),
     theMateType("Material"), theIonisationMeanExcitationEnergy(-1.),
     theState(kStateUndefined), theTemperature(STP_Temperature),
     thePressure(STP_Pressure)
{
}


// -------------------------------------------------------------------------
G4tgrMaterial::~G4tgrMaterial()
{
}


// -------------------------------------------------------------------------
void G4tgrMaterial::SetState( G4String val )
{ 
  if( val == "Undefined" ) {
    theState = kStateUndefined;
  } else if( val == "Solid" ) {
    theState = kStateSolid;
  } else if( val == "Liquid" ) {
    theState = kStateLiquid;
  } else if( val == "Gas" ) {
    theState = kStateGas;
  } else {

    G4Exception("G4tgrMaterial::SetState", "Wrong state", FatalErrorInArgument,
                "Only possible states are Undefined/Solid/Liquid/Gas!");
  }
}
