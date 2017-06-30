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
// Hadrontherapy advanced example for Geant4
// See more at: https://twiki.cern.ch/twiki/bin/view/Geant4/AdvancedExamplesHadrontherapy

#include <fstream>
#include <iostream>
#include <sstream>
#include <cmath>
#include <vector>

#include "HadrontherapyInteractionParameters.hh"
#include "HadrontherapyParameterMessenger.hh"
#include "HadrontherapyDetectorConstruction.hh"

#include "globals.hh"
#include "G4SystemOfUnits.hh"
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

HadrontherapyInteractionParameters::HadrontherapyInteractionParameters(G4bool wantMessenger): 
    nistEle(new G4NistElementBuilder(0)),										  
    nistMat(new G4NistMaterialBuilder(nistEle, 0)),									  
    data(G4cout.rdbuf()), 
    pMessenger(0),
    beamFlag(false)

{
    if (wantMessenger) pMessenger = new HadrontherapyParameterMessenger(this); 
}

HadrontherapyInteractionParameters::~HadrontherapyInteractionParameters()
{
    if (pMessenger) delete pMessenger; 
    delete nistMat; 
    delete nistEle; 
}

G4double HadrontherapyInteractionParameters::GetStopping (G4double ene, 
	                                                  const G4ParticleDefinition* pDef, 
						          const G4Material* pMat,
							  G4double dens)
{
    if (dens) return ComputeTotalDEDX(ene, pDef, pMat)/dens;
    return ComputeTotalDEDX(ene, pDef, pMat);
}
bool HadrontherapyInteractionParameters::GetStoppingTable(const G4String& vararg)
{
	// Check arguments
	if ( !ParseArg(vararg)) return false;
	// Clear previous energy & mass sp vectors
        energy.clear(); 
        massDedx.clear();
	// log scale 
	if (kinEmin != kinEmax && npoints >1)
	{
	    G4double logmin = std::log10(kinEmin);
	    G4double logmax = std::log10(kinEmax); 
	    G4double en;
	    // uniform log space
            for (G4double c = 0.; c < npoints; c++)
	       {
		    en = std::pow(10., logmin + ( c*(logmax-logmin)  / (npoints - 1.)) );  
		    energy.push_back(en/MeV);
		    dedxtot =  ComputeTotalDEDX (en, particle, material);
	            massDedx.push_back ( (dedxtot / density)/(MeV*cm2/g) );
	       }
	}
	else // one point only
	{
	    energy.push_back(kinEmin/MeV);
	    dedxtot =  ComputeTotalDEDX (kinEmin, particle, material);
	    massDedx.push_back ( (dedxtot / density)/(MeV*cm2/g) );
	}

    G4cout.precision(6);  
    data <<  "MeV             " << "MeV*cm2/g      " << particle << " (into " << 
	    material << ", density = " << G4BestUnit(density,"Volumic Mass") << ")" << G4endl;
    data << G4endl;
    data << std::left << std::setfill(' ');
    for (size_t i=0; i<energy.size(); i++){
		data << std::setw(16) << energy[i] << massDedx[i] << G4endl;
	}
    outfile.close();
    // This will plot

// Info to user
    G4String ofName = (filename == "") ? "User terminal": filename;
    G4cout << "User choice:\n";
    G4cout << "Kinetic energy lower limit= "<< G4BestUnit(kinEmin,"Energy") << 
	      ", Kinetic energy upper limit= " << G4BestUnit(kinEmax,"Energy") << 
	         ", npoints= "<< npoints << ", particle= \"" << particle << 
		 "\", material= \"" << material << "\", filename= \""<< 
		 ofName << "\"" << G4endl;
    return true;
}
// Search for user material choice inside G4NistManager database
G4Material* HadrontherapyInteractionParameters::GetNistMaterial(G4String mat)
{
    Pmaterial = G4NistManager::Instance()->FindOrBuildMaterial(mat);
    if (Pmaterial) density = Pmaterial -> GetDensity(); 
    return Pmaterial;
}
// Parse arguments line
bool HadrontherapyInteractionParameters::ParseArg(const G4String& vararg)
{
  kinEmin = kinEmax = npoints = 0.;
  particle = material = filename = "";
  // set internal variables
  std::istringstream strParam(vararg);
  // TODO here check for number and parameters consistency 
  strParam >> std::skipws >> material >> kinEmin >> kinEmax >> npoints >> particle >> filename;
  // npoints must be an integer!
  npoints = std::floor(npoints); 

// Check that  kinEmax >= kinEmin > 0 &&  npoints >= 1 
// TODO NIST points and linear scale
   if (kinEmax == 0. && kinEmin > 0. ) kinEmax = kinEmin;
   if (kinEmax == 0. && kinEmin == 0. ) kinEmax = kinEmin = 1.*MeV;
   if (kinEmax < kinEmin) 
   {
     G4cout << "WARNING: kinEmin must not exceed kinEmax!" << G4endl;
     G4cout << "Usage: /parameter/command  material EkinMin EKinMax nPoints [particle] [output filename]" << G4endl; 	
     return false;
   }
  if (npoints < 1) npoints = 1;

    // check if element/material is into database
    if (!GetNistMaterial(material) )
	   {
	     G4cout << "WARNING: material \"" << material << "\" doesn't exist in NIST elements/materials"
	 	       " table [$G4INSTALL/source/materials/src/G4NistMaterialBuilder.cc]" << G4endl; 
	     G4cout << "Use command \"/parameter/nist\" to see full materials list" << G4endl; 
	     return false;
	   }
    // Check for particle
    if (particle == "") particle = "proton"; // default to "proton"
    else if ( !FindParticle(particle) )
	   {
		G4cout << "WARNING: Particle \"" << particle << "\" isn't supported." << G4endl;
		G4cout << "Try the command \"/particle/list\" to get full supported particles list." << G4endl;
		G4cout << "If you are interested in an ion that isn't in this list you must give it to the particle gun."
		          "\nTry the commands:\n/gun/particle ion"
			  "\n/gun/ion <atomic number> <mass number> <[charge]>" << G4endl << G4endl;
		return false;
	   }
    // start physics by forcing a G4RunManager::BeamOn(): 
    BeamOn();
    // Set output file
    if( filename != "" ) 
       {
          outfile.open(filename,std::ios_base::trunc); // overwrite existing file
          data.rdbuf(outfile.rdbuf());
       }
    else data.rdbuf(G4cout.rdbuf());	// output is G4cout!                
    return true;
}
// Force physics tables build
void HadrontherapyInteractionParameters::BeamOn()
{
    // first check if RunManager is above G4State_Idle 
    G4StateManager* mState = G4StateManager::GetStateManager();
    G4ApplicationState  aState = mState -> GetCurrentState(); 
    if ( aState <= G4State_Idle && beamFlag == false)
	 {
	    G4cout << "Issuing a G4RunManager::beamOn()... "; 
	    G4cout << "Current Run State is " << mState -> GetStateString( aState ) << G4endl; 
	    G4RunManager::GetRunManager() -> BeamOn(0);
	    beamFlag = true;
	 }

}
// print a list of Nist elements and materials
void HadrontherapyInteractionParameters::ListOfNistMaterials(const G4String& vararg)
{
/*
 $G4INSTALL/source/materials/src/G4NistElementBuilder.cc
 You can also construct a new material by the ConstructNewMaterial method:
 see $G4INSTALL/source/materials/src/G4NistMaterialBuilder.cc
*/
    // Get simplest full list
     if (vararg =="list")
	  {
	    const std::vector<G4String>& vec =  nistMat -> GetMaterialNames(); 
	    for (size_t i=0; i<vec.size(); i++)
	     {
	        G4cout << std::setw(12) << std::left << i+1 << vec[i] << G4endl;
	     }
	    G4cout << G4endl;
          }
    else if (vararg =="all" || vararg =="simple" || vararg =="compound" || vararg =="hep" )
	{
		nistMat -> ListMaterials(vararg);
	}
}


