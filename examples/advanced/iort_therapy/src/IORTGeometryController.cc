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

#include "IORTGeometryController.hh"
#include "IORTDetectorConstruction.hh"
#include "Collimator40BeamLine.hh"
#include "Collimator50BeamLine.hh"
#include "Collimator60BeamLine.hh"
#include "Collimator70BeamLine.hh"
#include "Collimator80BeamLine.hh"
#include "Collimator100BeamLine.hh"
#include "G4RunManager.hh"
#include "IORTGeometryMessenger.hh"

/////////////////////////////////////////////////////////////////////////////
IORTGeometryController::IORTGeometryController()
{}

/////////////////////////////////////////////////////////////////////////////
IORTGeometryController::~IORTGeometryController()
{}

/////////////////////////////////////////////////////////////////////////////
void IORTGeometryController::SetGeometry(G4String name)
{

    if (name == "coll100")
    {
	registerGeometry(new Collimator100BeamLine());
	G4cout <<"Collimator 100 geometry activated" << G4endl;
      
    }
    else if (name == "coll80")
    {
	registerGeometry(new Collimator80BeamLine());
	G4cout <<"Collimator 80 geometry activated" << G4endl;
      
    } 
     
    else if (name == "coll70")
    {
	registerGeometry(new Collimator70BeamLine());
	G4cout <<"Collimator 70 geometry activated" << G4endl;
      
    } 
    else if (name == "coll50")
    {
	registerGeometry(new Collimator50BeamLine());
	G4cout <<"Collimator 50 geometry activated" << G4endl;
      
    } 
    else if (name == "coll40")
    {
	registerGeometry(new Collimator40BeamLine());
	G4cout <<"Collimator 40 geometry activated" << G4endl;
      
    } 
    else if(name == "default") 
        
    {
	registerGeometry(new Collimator60BeamLine());
    } 
    else
    {
	G4cout <<"Unknown geometry: " << name << ". Geometry not changed." << G4endl;
    }
}
	
/////////////////////////////////////////////////////////////////////////////
void IORTGeometryController::registerGeometry(G4VUserDetectorConstruction *detector)
{
	G4RunManager *runManager = G4RunManager::GetRunManager();
	runManager->SetUserInitialization(detector);
	runManager->GeometryHasBeenModified();
}
/////////////////////////////////////////////////////////////////////////////

void IORTGeometryController::UpdateGeometry()
{
        G4RunManager::GetRunManager() -> GeometryHasBeenModified();

}




