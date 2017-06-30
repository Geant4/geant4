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

#ifndef HadrontherapyMatrix_H
#define HadrontherapyMatrix_H 1
#include <G4ParticleDefinition.hh>
#include "globals.hh"
#include <vector>
#include <fstream>
#include "g4csv.hh"


#ifndef HADRONTHERAPYANALYSISMANAGER_HH
#define HADRONTHERAPYANALYSISMANAGER_HH 1

class HadrontherapyAnalysisFileMessenger;

/**
 * A class for connecting the simulation to an analysis package.
 */
class HadrontherapyAnalysisManager
{
private:
    /**
     * Analysis manager is a singleton object (there is only one instance).
     * The pointer to this object is available through the use of the method GetInstance();
     *
     * @see GetInstance
     */
    HadrontherapyAnalysisManager();
    
    
    
public:
    ~HadrontherapyAnalysisManager();
    
    /**
     * Get the pointer to the analysis manager.
     */
    static HadrontherapyAnalysisManager* GetInstance();
    
    
    
    static HadrontherapyAnalysisManager* instance;
    HadrontherapyAnalysisFileMessenger* fMess;
    
};

#endif

// The information: energy deposit and position in the phantom
// is stored in a matrix

// type struct useful to store nucludes data


struct ion
{ 
  G4bool isPrimary;   // true if particle is primary
  G4int PDGencoding;  // Particle data group id for the particle
  //G4String extName; //  AZ[excitation energy]: like He3[1277.4], He4[0.0], Li7[231.4], ...
  G4String name;   	 // simple name without excitation energy: He3, He4, Li7, ...
  std::string::size_type len; 	 // name length
  G4int Z; 		 // atomic number
  G4int A; 		 // mass number
  G4double *dose; 	 // pointer to dose matrix
  unsigned int    *fluence;  // pointer to fluence matrix
  //friend bool operator<(const ion& a, const ion& b) {return (a.Z == b.Z) ? b.A < a.A : b.Z < a.Z ;}
  G4bool operator<(const ion& a) const{return (this->Z == a.Z) ? this-> A < a.A : this->Z < a.Z ;}
};

class HadrontherapyMatrix 
{
private:
  HadrontherapyMatrix(G4int numberOfVoxelAlongX, 
		      G4int numberOfVoxelAlongY, 
		      G4int numberOfVoxelAlongZ,
		      G4double massOfVoxel); //< this is supposed to be a singleton


public:

  ~HadrontherapyMatrix();
  // Get object instance only
  static HadrontherapyMatrix* GetInstance();
  // Make & Get instance
  static HadrontherapyMatrix* GetInstance(G4int nX, G4int nY, G4int nZ, G4double mass);

  static G4bool secondary;
  // Full list of generated nuclides
    
    
    
  void PrintNuclides(); 
  // Hit array marker (useful to avoid multiple counts of fluence)
  void ClearHitTrack();
  G4int* GetHitTrack(G4int i, G4int j, G4int k);

  // All the elements of the matrix are initialised to zero
  void Initialize(); 
  void Clear();
  // Fill DOSE/fluence matrix for particle: 
  // if fluence parameter is true then fluence at voxel (i, j, k) is increased 
  // else energyDeposit fill the dose matrix for voxel (i,j,k) 
  G4bool Fill(G4int, G4ParticleDefinition* particleDef, G4int i, G4int j, G4int k, G4double energyDeposit, G4bool fluence=false); 

  // Fill TOTAL DOSE matrix for primary particles only 
  void Fill(G4int i, G4int j, G4int k, G4double energyDeposit);
  // The matrix is filled with the energy deposit 
  // in the element corresponding to the voxel of the phantom where
  // the energy deposit was registered
  
  // Store the information of the matrix in a ntuple and in 
  // a 1D Histogram
  //void TotalEnergyDeposit();
   
  // Store single matrix data to filename 
  void StoreMatrix(G4String file, void* data,size_t psize);
  // Store all fluence data to filenames
  void StoreFluenceData();
  // Store all dose data to filenames
  void StoreDoseData();

  // Store all data (except the total dose) to ONE filename
  void StoreDoseFluenceAscii(G4String filename = "");
  

    
  inline G4int Index(G4int i, G4int j, G4int k) { return (i * numberOfVoxelAlongY + j) * numberOfVoxelAlongZ + k; } 
  // Get a unique index from  a three dimensional one 

  G4double * GetMatrix(){return matrix;}

  G4int GetNvoxel(){return numberOfVoxelAlongX*numberOfVoxelAlongY*numberOfVoxelAlongZ;}
  // Total number of voxels read only access  
  G4int GetNumberOfVoxelAlongX(){return numberOfVoxelAlongX;}
  G4int GetNumberOfVoxelAlongY(){return numberOfVoxelAlongY;}
  G4int GetNumberOfVoxelAlongZ(){return numberOfVoxelAlongZ;}
private:

  static HadrontherapyMatrix* instance;
  G4int numberOfVoxelAlongX;
  G4int numberOfVoxelAlongY;
  G4int numberOfVoxelAlongZ;
  G4double massOfVoxel;

  G4double* matrix;
  G4int* hitTrack;
  G4String stdFile, filename;
  std::ofstream ofs;

  // Dose&fluence data store 
  std::vector <ion> ionStore;
  // want secondary particles?
  G4double doseUnit;
};
#endif

