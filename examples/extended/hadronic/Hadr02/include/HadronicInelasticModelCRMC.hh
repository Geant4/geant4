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
/// \file hadronic/Hadr02/include/HadronicInelasticModelCRMC.hh
/// \brief Definition of the HadronicInelasticModelCRMC class
//
//
// ------------------------------------------------------------
//
//              CRMC interface to GEANT4
//              for more details on CRMC, see:
//              https://web.ikp.kit.edu/rulrich/crmc.html
//
//
// Author:      Andrii Tykhonov  (University of Geneva)
// Email:       andrii.tykhonov@cern.ch
// Created:     14.02.2018
//
//
// (
//   A few, trivial modifications made by A. Ribon in May 2021
//   in order to use it inside the Geant4 example Hadr02 :
//   -  Copied here, instead of using the original class
//      G4HadronicInelasticModelCRMC as it is distributed
//      with crmc-svn-geant4, to avoid some problems with CMake
//      which is unable to find the library libGeantCrmc.
//   -  Renamed the class as HadronicInelasticModelCRMC,
//      to follow the Geant4 convention that classes whose
//      names start with "G4" are only those distributed in
//      source/ .
//   -  Fixed a few compilation warnings .
//   -  Set the seed by hand, because the method
//        CLHEP::HepRandom::getTheSeed()
//      returns 0 which is not accepted.
// )
// ------------------------------------------------------------

#ifndef  HadronicInelasticModelCRMC_h
#define  HadronicInelasticModelCRMC_h

#include "G4HadronicInteraction.hh"
#include "G4HadFinalState.hh"
#include "G4SystemOfUnits.hh"

#include "CRMCinterface.h"
#include <string>


extern CRMCdata gCRMC_data;

class G4HadFinalState;
class G4ParticleTable;
class G4IonTable;


class  G4ParticleDefinition; //class G4DynamicParticle; 


class HadronicInelasticModelCRMC : public G4HadronicInteraction
{
public:
	//! model: 
	//!         0  :  EPOS LHC
	//!         1  :  EPOS 1.99
	//!         12 :  DPMJET3
	HadronicInelasticModelCRMC(int model, const G4String& modelName);
	~HadronicInelasticModelCRMC();

	G4HadFinalState * ApplyYourself (const G4HadProjectile &aTrack, G4Nucleus &targetNucleus);
	G4bool 	IsApplicable (const G4HadProjectile &, G4Nucleus &);

	void SetPrintDebug(bool printdebug) {fPrintDebug = printdebug;}
	G4ParticleDefinition* GetParticleDefinition(long particle_id,int& error_code); 
	void SplitMultiNeutrons(CRMCdata& CRMC_data);
	bool IsMultiNeutron(int Z, int A);

	virtual const std::pair<G4double, G4double> GetFatalEnergyCheckLevels() const {
		// possible energy non-coservations of up to 1 TeV are ignored
		return std::pair<G4double, G4double>( 10.0*perCent, 1000.0*GeV );
	}
private:
	CRMCinterface* fInterface;
	//CRMCdata fCRMCdata;
	int fTypeOutput;
	G4HadFinalState* finalState;
	G4ParticleTable* fParticleTable;
	G4IonTable*      fIonTable;

	//std::vector<G4ParticleDefinition*> fParticleDefinitions;
	//std::vector<G4DynamicParticle*>    fDynamicParticles;
	bool fPrintDebug;

	//
	std::string GetCrmcParamPath();


};

#endif
