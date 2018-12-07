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
// This is the *BASIC* version of IORT, a Geant4-based application
//
// Main Authors: G.Russo(a,b), C.Casarino*(c), G.C. Candiano(c), G.A.P. Cirrone(d), F.Romano(d)
// Contributor Authors: S.Guatelli(e)
// Past Authors: G.Arnetta(c), S.E.Mazzaglia(d)
//    
//   (a) Fondazione Istituto San Raffaele G.Giglio, Cefalù, Italy
//   (b) IBFM-CNR , Segrate (Milano), Italy
//   (c) LATO (Laboratorio di Tecnologie Oncologiche), Cefalù, Italy
//   (d) Laboratori Nazionali del Sud of the INFN, Catania, Italy
//   (e) University of Wollongong, Australia
//
//   *Corresponding author, email to carlo.casarino@polooncologicocefalu.it
//////////////////////////////////////////////////////////////////////////////////////////////

#include "Collimator80BeamLineMessenger.hh"
#include "Collimator80BeamLine.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAString.hh"


    Collimator80BeamLineMessenger::Collimator80BeamLineMessenger(Collimator80BeamLine* beamLine)
:collimator80(beamLine)

{
    
    beamLineDir = new G4UIdirectory("/beamLine/");
    beamLineDir -> SetGuidance("set specification of range shifter");  

    FinalCollimatorIORTDir = new G4UIdirectory("/beamLine/FinalCollimatorIORT/");
    FinalCollimatorIORTDir -> SetGuidance("set specification of final collimator");  
    
    innerRadiusFinalCollimatorIORTCmd = new G4UIcmdWithADoubleAndUnit("/beamLine/FinalCollimatorIORT/halfInnerRad",this);
    innerRadiusFinalCollimatorIORTCmd -> SetGuidance("Set size of inner radius ( max 21.5 mm)");
    innerRadiusFinalCollimatorIORTCmd -> SetParameterName("Size",false);
    innerRadiusFinalCollimatorIORTCmd -> SetDefaultUnit("mm");  
    innerRadiusFinalCollimatorIORTCmd -> SetUnitCandidates("mm cm m");  
    innerRadiusFinalCollimatorIORTCmd -> AvailableForStates(G4State_Idle);

    OuterRadiusFinalCollimatorIORTCmd = new G4UIcmdWithADoubleAndUnit("/beamLine/FinalCollimatorIORT/halfOuterRad",this);
    OuterRadiusFinalCollimatorIORTCmd -> SetGuidance("Set size of outer radius ( max 21.5 mm)");
    OuterRadiusFinalCollimatorIORTCmd -> SetParameterName("Size",false);
    OuterRadiusFinalCollimatorIORTCmd -> SetDefaultUnit("mm");  
    OuterRadiusFinalCollimatorIORTCmd -> SetUnitCandidates("mm cm m");  
    OuterRadiusFinalCollimatorIORTCmd -> AvailableForStates(G4State_Idle);

}

Collimator80BeamLineMessenger::~Collimator80BeamLineMessenger()
{ 

    delete OuterRadiusFinalCollimatorIORTCmd; 
    delete innerRadiusFinalCollimatorIORTCmd; 
    delete FinalCollimatorIORTDir;  
    delete beamLineDir; 
     
 
}




void Collimator80BeamLineMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{ 

 
    if( command == innerRadiusFinalCollimatorIORTCmd )
    { collimator80 -> SetInnerRadiusFinalCollimatorIORT
	(innerRadiusFinalCollimatorIORTCmd -> GetNewDoubleValue(newValue));}

    else if( command == OuterRadiusFinalCollimatorIORTCmd )
    { collimator80 -> SetOuterRadiusFinalCollimatorIORT
	(OuterRadiusFinalCollimatorIORTCmd -> GetNewDoubleValue(newValue));}

}

