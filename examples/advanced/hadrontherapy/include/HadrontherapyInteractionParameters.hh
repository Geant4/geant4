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
// $Id: HadrontherapyInteractionParameters.hh;
//



#ifndef HadrontherapyInteractionParameters_H
#define HadrontherapyInteractionParameters_H 1

#include "G4EmCalculator.hh"
#include "G4NistMaterialBuilder.hh"
#include "G4NistElementBuilder.hh"

class HadrontherapyDetectorConstruction;
class HadrontherapyParameterMessenger; 
class HadrontherapyInteractionParameters : public G4EmCalculator 
{
public:

    HadrontherapyInteractionParameters();
	~HadrontherapyInteractionParameters();

// Get data for Mass SP (MeV*cm2/g)   
// into G4NistMaterialBuilder class materials
// User must provide: material kinetic energy lower limit, kinetic energy upper limit, number of points to retrieve,
// [particle], [output filename].

	bool GetStoppingTable (const G4String& vararg);
    void ListOfNistMaterials (const G4String& vararg);
    void BeamOn();
    bool ParseArg (const G4String& vararg);	

private:
	G4Material* GetNistMaterial(G4String material);

	G4NistElementBuilder* nistEle;
	G4NistMaterialBuilder* nistMat;
	G4double kinEmin, kinEmax, npoints;
	G4String particle, material, filename; 
	std::ofstream outfile;
	std::ostream data;
    G4Material* Pmaterial;
	G4double density;
    G4EmCalculator* emCal;
    HadrontherapyParameterMessenger* pMessenger; 
	bool beamFlag;
};
#endif

