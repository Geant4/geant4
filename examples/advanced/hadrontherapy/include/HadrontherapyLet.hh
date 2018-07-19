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

#ifndef HadrontherapyLet_h
#define HadrontherapyLet_h 1
#endif

#include "G4ParticleDefinition.hh"
#include "globals.hh"
#include <fstream>
#include <vector>
#include <string>

#include "g4csv.hh"
#include "HadrontherapyMatrix.hh"
struct ionLet 
{ 
    G4bool isPrimary;	    // True if particle is primary
    G4int PDGencoding;      // Particle data group id for the particle
    G4String fullName;      // AZ[excitation energy]: like He3[1277.4], He4[0.0], Li7[231.4], ...
    G4String name;          // simple name without excitation energy: He3, He4, Li7, ...
    G4int Z;                // atomic number
    G4int A;		    // mass number
    G4double *letDN , *letDD; // Track averaged LET and Dose averaged LET 
    //friend bool operator<(const ionLet& a, const ionLet& b) {return (a.Z == b.Z) ? b.A < a.A : b.Z < a.Z ;}
    G4bool operator<(const ionLet& a) const{return (this->Z == a.Z) ? this-> A < a.A : this->Z < a.Z ;}
    // For isotopes sort by the mass number, else sort by the atomic one.
};

class G4Material;
class HadrontherapyMatrix;
class HadrontherapyPrimaryGeneratorAction;
class HadrontherapyInteractionParameters;
class HadrontherapyDetectorConstruction;

class HadrontherapyLet
{
    private:
	HadrontherapyLet(HadrontherapyDetectorConstruction*); 

    public:
	~HadrontherapyLet();
	static HadrontherapyLet* GetInstance(HadrontherapyDetectorConstruction*);
	static HadrontherapyLet* GetInstance();
	static G4bool doCalculation;
	void Initialize();
	void Clear();

void Fill(G4int i, G4int j, G4int k, G4double DE, G4double DX);
void FillEnergySpectrum (G4int trackID,
		 		G4ParticleDefinition* particleDef,
				/*G4double kinEnergy,*/
				G4double DE,
				G4double DX,
				G4int i, G4int j, G4int k); 
	void LetOutput(); 
	void StoreLetAscii();


    private:
	static HadrontherapyLet *instance;
	HadrontherapyPrimaryGeneratorAction* pPGA;
	// Detector material
	G4Material* detectorMat;  
	G4double density;
	G4String filename;

	std::ofstream ofs;
        std::ofstream stopFile;
	HadrontherapyMatrix *matrix;
	G4int nVoxels, numberOfVoxelAlongX, numberOfVoxelAlongY, numberOfVoxelAlongZ ;
	G4double primaryEnergy, energyLimit, binWidth; 
	G4int nBins;
	G4double  nT, dT, nD, dD;
	G4double  nSecondaryT, nSecondaryD, dSecondaryT, dSecondaryD;
	G4double  nPrimaryT, nPrimaryD, dPrimaryT, dPrimaryD ;

	G4double *secondaryLetT,  *secondaryLetD, *totalLetT,  *DtotalLetD, *totalLetD;
	G4String nome_file;

	std::vector<ionLet> ionLetStore;
};
