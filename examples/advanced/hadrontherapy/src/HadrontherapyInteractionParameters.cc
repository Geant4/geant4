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
// $Id: HadrontherapyInteractionParameters.cc;
//

#include "HadrontherapyInteractionParameters.hh"
#include "HadrontherapyParameterMessenger.hh"
#include "HadrontherapyDetectorConstruction.hh"

#include "G4NistMaterialBuilder.hh"
#include "G4UnitsTable.hh"
#include "G4UImanager.hh"
#include "G4RunManager.hh"
#include "G4LossTableManager.hh"
#include "G4Material.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4NistManager.hh"
#include "G4Element.hh"

#include "G4StateManager.hh"

#include "globals.hh"
#include <fstream>
#include <iostream>
#include <sstream>
#include <math.h>
#include <unistd.h>

#include <vector>


HadrontherapyInteractionParameters::HadrontherapyInteractionParameters(): 
	data(std::cout.rdbuf()), emCal(new G4EmCalculator),
	pMessenger(new HadrontherapyParameterMessenger(this)) 
{
}

HadrontherapyInteractionParameters::~HadrontherapyInteractionParameters()
{
	delete emCal;
    delete pMessenger; 
}
bool HadrontherapyInteractionParameters::GetStoppingTable(G4String vararg)
{
	// load & check arguments
	if ( !ParseArg(vararg)) { return false; }

	// check data for writing and  kinEmax > kinEmin > 0 && npoints > 1, integer 
	// define two dimensional dinamic array for energy-stopping power
    
	std::vector<G4double> energy;
	std::vector<G4double> massDedx;
	G4double dedxtot ;


	// log scale (remember to add linear scale algo also)
    
	// compute data
	if (kinEmin != kinEmax){
        G4double logmin = std::log10(kinEmin);
	    G4double logmax = std::log10(kinEmax); 
	    G4double en;
		// uniform logarithm space
    for (G4double c = 0; c <npoints; c++ ){
		en = std::pow(10, logmin + ( c*(logmax-logmin)  / (npoints - 1)) );  
	    energy.push_back(en);
		dedxtot =  emCal -> ComputeTotalDEDX (en, particle, material);
		massDedx.push_back ( dedxtot / density );
	   }
	}
	
	//data <<  particle << " (into " << material << ")" << std::setw(12) << G4BestUnit(density,"Volumic Mass") << G4endl;
    data << std::left << std::setfill(' ');
    for (size_t i=0; i<energy.size(); i++){
		data << std::setw(12) << energy[i]/MeV << std::setw(12) << massDedx[i]/(MeV*cm2/g) << G4endl;
		//data << std::left << std::setw(10) << G4BestUnit(energy[i],"Energy") << 
		//std::right << std::setw(17) <<  G4BestUnit(massDedx[i],"Energy*Surface/Mass") << G4endl;
	}

	outfile.close();
	return true;

}


// search for user material choice inside G4NistManager database
G4Material* HadrontherapyInteractionParameters::GetNistMaterial(G4String material)
{
    Pmaterial = G4NistManager::Instance()->FindOrBuildMaterial(material);
    if (Pmaterial) {
		density = Pmaterial -> GetDensity();}
	return Pmaterial;
}


bool HadrontherapyInteractionParameters::ParseArg(const G4String& vararg)
{
  
  kinEmin = kinEmax = npoints = 0;
  particle = material = filename = "";
  // set internal variables
  std::istringstream strParam(vararg);
  // here check for number and parameter consistency!! 
  strParam >> std::skipws >> material >> kinEmin >> kinEmax >> npoints >> particle >> filename;
  // npoints must be an integer
  npoints = std::floor(npoints); 

		// check energy range 
	if (kinEmax == 0 && kinEmin == 0) {}// Set NIST points!
	else if (kinEmax == 0) kinEmax = kinEmin;
	else if (kinEmax < kinEmin) {
    G4cout << "WARNING: kinEmin must not exceed kinEmax!" << G4endl;
    //G4cout << "Usage: /parameter/command  material [kinetic Emin, kinetic Emax nPoints] [particle] [output filename]" << G4endl; 	
    G4cout << "Usage: /parameter/command  material kinetic Emin kinetic Emax nPoints [particle] [output filename]" << G4endl; 	
    return false;
	if (npoints < 2){
		    G4cout << "WARNING: you must request at least 2 points!" << G4endl;
		    return false;
	        }       
	}
	// check if material is here!
	
	if (!GetNistMaterial(material) ){
        G4cout << "WARNING: material \"" << material << "\" doesn't exist in NIST materials table [source/materials/src/G4NistMaterialBuilder.cc]" << G4endl; 
        G4cout << "Type /parameter/materials to see full materials list" << G4endl; 
		return false;}
    // check if particle is any of: alpha, proton, e-
    if (particle == "") particle = "proton";
		else if (particle != "proton" && particle != "alpha" && particle != "e-"){
		G4cout << "WARNING: Particle \"" << particle << "\" isn't supported" << G4endl;
		return false;
		}
    

    // start physic
	BeamOn();

    // Set output file
    // check data for writing and [ inf > kinEmax > kinEmin > 0 ] [ npoints > 1 integer ]
	// define two dimensional dinamic array for energy-stopping power
	if( filename != "" ) 
	   {
	      outfile.open(filename,std::ios_base::trunc); // overwrite existing file
	      data.rdbuf(outfile.rdbuf());
       }
	else data.rdbuf(std::cout.rdbuf());	// output is G4cout                

  G4cout << "User choice:\n";
  G4cout << "kinEmin= "<< G4BestUnit(kinEmin,"Energy") << " kinEmax= " << G4BestUnit(kinEmax,"Energy") << 
	         " npoints= "<< npoints << " particle=\"" << particle << "\" material=\"" << material <<
			 "\" filename=\"" << filename <<"\""<< G4endl;
  G4cout.precision(6);  
  return true;
}
// Initialize physic
//
bool HadrontherapyInteractionParameters::beamFlag(false);
void HadrontherapyInteractionParameters::BeamOn()
{
	// first check if RunManager is in G4State_Idle 
	//
	G4StateManager* mState = G4StateManager::GetStateManager();
    G4ApplicationState  aState = mState -> GetCurrentState(); 
	if ( aState <= G4State_Idle ){
    if (beamFlag == false){ 
//  G4cout << "Run State " << mState -> GetStateString( aState ) << G4endl; 
    G4RunManager::GetRunManager() -> BeamOn(0);
    //G4UImanager* UI = G4UImanager::GetUIpointer();  
    //UI -> ApplyCommand("/run/beamOn");
	beamFlag = true;
	}
	}

}
// print full list of Nist elements and materials
void HadrontherapyInteractionParameters::ListOfNistMaterials()
{
	
	//source/materials/src/G4NistMaterialBuilder.cc
	//~/cvs/geant4-09-02-ref-04/source/materials/src/G4NistElementBuilder.cc
    //./global/management/src/G4UnitsTable.cc	
	//
	G4NistMaterialBuilder* nistmat = new G4NistMaterialBuilder(new G4NistElementBuilder(0), 0);

	const std::vector<G4String>& vec =  nistmat -> GetMaterialNames(); 
    for (size_t i=0; i<vec.size(); i++){
		//G4cout << std::setw(24) << std::left << vec[i] ;
		G4cout << std::setw(24) << std::left << i+1 << vec[i] << G4endl;
	}
	G4cout << G4endl;
/*
    // full list (formula, density, etc)
    nistmat -> ListMaterials("simple");
	G4cout << "....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......\n";
    nistmat -> ListMaterials("compound");
	G4cout << "....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......\n";
    nistmat -> ListMaterials("hep");
	G4cout << "....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......\n";
*/
}


