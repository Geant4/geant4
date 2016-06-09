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
// Class Description
// Name list of Elements for a high precision (based on evaluated data
// libraries) description of themal neutron scattering below 4 eV;
// Based on Thermal neutron scattering files
// from the evaluated nuclear data files ENDF/B-VI, Release2
// To be used in your physics list in case you need this physics.
// In this case you want to register an object of this class with
// the corresponding process.
// Class Description - End

// 15-Nov-06 First implementation is done by T. Koi (SLAC/SCCS)

#include "G4NeutronHPThermalScatteringNames.hh"
#include "G4Neutron.hh"
#include "G4ElementTable.hh"
//#include "G4NeutronHPData.hh"



G4NeutronHPThermalScatteringNames::G4NeutronHPThermalScatteringNames()
{
   names.insert ( std::pair < G4String , G4String > ( "TS_H_of_Water" , "h_water" ) ); 
   names.insert ( std::pair < G4String , G4String > ( "TS_H_of_Polyethylene" , "h_polyethylene" ) ); 
   names.insert ( std::pair < G4String , G4String > ( "TS_C_of_Graphite" , "graphite" ) ); 
}



G4NeutronHPThermalScatteringNames::~G4NeutronHPThermalScatteringNames()
{
;
}



G4bool G4NeutronHPThermalScatteringNames::IsThisThermalElement( G4String aname)
{
   G4bool result = false;
   if ( names.find ( aname ) != names.end() ) result = true; 
   return result;
}
