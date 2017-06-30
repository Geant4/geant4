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
#include <iomanip>

#include "HadrontherapyMatrix.hh"
#include "HadrontherapyPrimaryGeneratorAction.hh"
#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "G4RunManager.hh"
#include "G4ParticleGun.hh"
#include "HadrontherapySteppingAction.hh"


#include "HadrontherapyAnalysisFileMessenger.hh"
#include "G4SystemOfUnits.hh"
#include <time.h>

HadrontherapyAnalysisManager* HadrontherapyAnalysisManager::instance = 0;
HadrontherapyAnalysisManager::HadrontherapyAnalysisManager()


{
    fMess = new HadrontherapyAnalysisFileMessenger(this);
}
HadrontherapyAnalysisManager::~HadrontherapyAnalysisManager()
{
    delete fMess;
}

HadrontherapyAnalysisManager* HadrontherapyAnalysisManager::GetInstance(){
    
    if (instance == 0) instance = new HadrontherapyAnalysisManager;
    return instance;
}



HadrontherapyMatrix* HadrontherapyMatrix::instance = NULL;
G4bool HadrontherapyMatrix::secondary = false;


// Only return a pointer to matrix
HadrontherapyMatrix* HadrontherapyMatrix::GetInstance()
{
    return instance;
}
// This STATIC method delete (!) the old matrix and rewrite a new object returning a pointer to it
// TODO A check on the parameters is required!
HadrontherapyMatrix* HadrontherapyMatrix::GetInstance(G4int voxelX, G4int voxelY, G4int voxelZ, G4double mass)
{
    if (instance) delete instance;
    instance = new HadrontherapyMatrix(voxelX, voxelY, voxelZ, mass);
    instance -> Initialize();
    return instance;
}


HadrontherapyMatrix::HadrontherapyMatrix(G4int voxelX, G4int voxelY, G4int voxelZ, G4double mass):
stdFile("Dose.out"),
doseUnit(gray)

{
    // Number of the voxels of the phantom
    // For Y = Z = 1 the phantom is divided in slices (and not in voxels)
    // orthogonal to the beam axis
    numberOfVoxelAlongX = voxelX;
    numberOfVoxelAlongY = voxelY;
    numberOfVoxelAlongZ = voxelZ;
    massOfVoxel = mass;
    
    
    // Create the dose matrix
    matrix = new G4double[numberOfVoxelAlongX*numberOfVoxelAlongY*numberOfVoxelAlongZ];
    if (matrix)
    {
        G4cout << "HadrontherapyMatrix: Memory space to store physical dose into " <<
        numberOfVoxelAlongX*numberOfVoxelAlongY*numberOfVoxelAlongZ <<
        " voxels has been allocated " << G4endl;
    }
    
    
    else G4Exception("HadrontherapyMatrix::HadrontherapyMatrix()", "Hadrontherapy0005", FatalException, "Can't allocate memory to store physical dose!");
    
    
    // Hit voxel (TrackID) marker
    // This array mark the status of voxel, if a hit occur, with the trackID of the particle
    // Must be initialized
    
    hitTrack = new G4int[numberOfVoxelAlongX*numberOfVoxelAlongY*numberOfVoxelAlongZ];
    ClearHitTrack();
}

/////////////////////////////////////////////////////////////////////////////
HadrontherapyMatrix::~HadrontherapyMatrix()
{
    delete[] matrix;
    delete[] hitTrack;
    Clear();
}

/////////////////////////////////////////////////////////////////////////////
void HadrontherapyMatrix::Clear()
{
    for (size_t i=0; i<ionStore.size(); i++)
    {
        delete[] ionStore[i].dose;
        delete[] ionStore[i].fluence;
    }
    ionStore.clear();
}

/////////////////////////////////////////////////////////////////////////////
// Initialise the elements of the matrix to zero

void HadrontherapyMatrix::Initialize()
{
    // Clear ions store
    Clear();
    // Clear dose
    for(int i=0;i<numberOfVoxelAlongX*numberOfVoxelAlongY*numberOfVoxelAlongZ;i++)
    {
        matrix[i] = 0;
    }
}

/////////////////////////////////////////////////////////////////////////////
// Print generated nuclides list


void HadrontherapyMatrix::PrintNuclides()
{
    for (size_t i=0; i<ionStore.size(); i++)
    {
        G4cout << ionStore[i].name << G4endl;
    }
}

/////////////////////////////////////////////////////////////////////////////
// Clear Hit voxel (TrackID) markers

void HadrontherapyMatrix::ClearHitTrack()
{
    for(G4int i=0; i<numberOfVoxelAlongX*numberOfVoxelAlongY*numberOfVoxelAlongZ; i++) hitTrack[i] = 0;
}


// Return Hit status
G4int* HadrontherapyMatrix::GetHitTrack(G4int i, G4int j, G4int k)
{
    return &(hitTrack[Index(i,j,k)]);
}



/////////////////////////////////////////////////////////////////////////////
// Dose methods...
// Fill DOSE/fluence matrix for secondary particles:
// If fluence parameter is true (default value is FALSE) then fluence at voxel (i, j, k) is increased.
// The energyDeposit parameter fill the dose matrix for voxel (i,j,k)
/////////////////////////////////////////////////////////////////////////////

G4bool HadrontherapyMatrix::Fill(G4int trackID,
                                 G4ParticleDefinition* particleDef,
                                 G4int i, G4int j, G4int k,
                                 G4double energyDeposit,
                                 G4bool fluence)
{
    
    if ( (energyDeposit <=0. && !fluence) || !secondary) return false;
    
    // Get Particle Data Group particle ID
    G4int PDGencoding = particleDef -> GetPDGEncoding();
    PDGencoding -= PDGencoding%10;
    
    // Search for already allocated data...
    for (size_t l=0; l < ionStore.size(); l++)
    {
        if (ionStore[l].PDGencoding == PDGencoding )
        {   // Is it a primary or a secondary particle?
            
            if ( (trackID ==1 && ionStore[l].isPrimary) || (trackID !=1 && !ionStore[l].isPrimary))
            {
                if (energyDeposit > 0.)
                    
                    ionStore[l].dose[Index(i, j, k)] += energyDeposit;
                
                // Fill a matrix per each ion with the fluence
                
                if (fluence) ionStore[l].fluence[Index(i, j, k)]++;
                return true;
            }
        }
    }
    
    G4int Z = particleDef-> GetAtomicNumber();
    G4int A = particleDef-> GetAtomicMass();
    G4String fullName = particleDef -> GetParticleName();
    G4String name = fullName.substr (0, fullName.find("[") ); // cut excitation energy
    
    // Let's put a new particle in our store...
   
    
    ion newIon =
    {
        (trackID == 1) ? true:false,
        PDGencoding,
        name,
        name.length(),
        Z,
        A,
        new G4double[numberOfVoxelAlongX * numberOfVoxelAlongY * numberOfVoxelAlongZ],
        new unsigned int[numberOfVoxelAlongX * numberOfVoxelAlongY * numberOfVoxelAlongZ]
    };
    
    
    // Initialize data
    if (newIon.dose && newIon.fluence)
    {
        for(G4int q=0; q<numberOfVoxelAlongX*numberOfVoxelAlongY*numberOfVoxelAlongZ; q++)
        {
            newIon.dose[q] = 0.;
            newIon.fluence[q] = 0;
        }
        
        if (energyDeposit > 0.) newIon.dose[Index(i, j, k)] += energyDeposit;
        if (fluence) newIon.fluence[Index(i, j, k)]++;
        
        ionStore.push_back(newIon);
        return true;
    }
    
    else // XXX Out of memory! XXX
    
    {
        return false;
    }
    
}

/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
// Methods to store data to filenames...
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
//
// General method to store matrix data to filename


void HadrontherapyMatrix::StoreMatrix(G4String file, void* data, size_t psize)
{
    if (data)
    {
        ofs.open(file, std::ios::out);
        if (ofs.is_open())
        {
            for(G4int i = 0; i < numberOfVoxelAlongX; i++)
                for(G4int j = 0; j < numberOfVoxelAlongY; j++)
                    for(G4int k = 0; k < numberOfVoxelAlongZ; k++)
                    {
                        G4int n = Index(i, j, k);
                        
                        if (psize == sizeof(unsigned int))
                        {
                            unsigned int* pdata = (unsigned int*)data;
                            
                            if (pdata[n])
                                
                                ofs << i << '\t' << j << '\t' << k << '\t' << pdata[n] << G4endl;
                            
                        }
                       
                        else if (psize == sizeof(G4double))
                        
                        {
                            G4double* pdata = (G4double*)data;
                            if (pdata[n]) ofs << i << '\t' << j << '\t' << k << '\t' << pdata[n] << G4endl;
                        }
                    }
            
            ofs.close();
        }
    }
}

// Store fluence per single ion in multiple files
void HadrontherapyMatrix::StoreFluenceData()
{
    for (size_t i=0; i < ionStore.size(); i++){
        StoreMatrix(ionStore[i].name + "_Fluence.out", ionStore[i].fluence, sizeof(unsigned int));
    }
}
// Store dose per single ion in multiple files
void HadrontherapyMatrix::StoreDoseData()
{
    
    for (size_t i=0; i < ionStore.size(); i++){
        StoreMatrix(ionStore[i].name + "_Dose.out", ionStore[i].dose, sizeof(G4double));
    }
}


/////////////////////////////////////////////////////////////////////////
// Store dose for all ions into a single file and into ntuples.
// Please note that this function is called via messenger commands
// defined in the HadrontherapyAnalysisFileMessenger.cc class file


void HadrontherapyMatrix::StoreDoseFluenceAscii(G4String file)
{
    
    
#define width 15L
    filename = (file=="") ? stdFile:file;
    
    // Sort like periodic table
    
    std::sort(ionStore.begin(), ionStore.end());
    G4cout << "Dose is being written to " << filename << G4endl;
    ofs.open(filename, std::ios::out);
    
    
    if (ofs.is_open())
    {
        // Write the voxels index and the list of particles/ions
        
        ofs << std::setprecision(6) << std::left <<
        "i\tj\tk\t";
        
        
        // Total dose
        ofs << std::setw(width) << "Dose(Gy)";
        
        
        if (secondary)
        {
            for (size_t l=0; l < ionStore.size(); l++)
            {
                G4String a = (ionStore[l].isPrimary) ? "_1":""; // is it a primary?
                ofs << std::setw(width) << ionStore[l].name + a <<
                std::setw(width) << ionStore[l].name  + a;
                
                
            }
            ofs << G4endl;
        }
        
  
        G4AnalysisManager* Analysis = G4AnalysisManager::Instance();
        
        Analysis ->SetVerboseLevel(1);
        Analysis ->SetFirstHistoId(1);
        Analysis ->SetFirstNtupleId(1);
        Analysis ->OpenFile("Dose");
        
        
        Analysis ->CreateNtuple("coordinate", "dose");
        
        Analysis ->CreateNtupleIColumn("i");//1
        Analysis ->CreateNtupleIColumn("j");//2
        Analysis ->CreateNtupleIColumn("k");//3
        Analysis ->CreateNtupleDColumn("totaldose");//6
        Analysis ->CreateNtupleIColumn("A");//4
        Analysis ->CreateNtupleIColumn("Z");//5
        Analysis ->CreateNtupleDColumn("Iondose");//6
        Analysis ->CreateNtupleDColumn("fluence");//7
        Analysis ->FinishNtuple();
        
        
        
        // Write data
        for(G4int i = 0; i < numberOfVoxelAlongX; i++)
            for(G4int j = 0; j < numberOfVoxelAlongY; j++)
                for(G4int k = 0; k < numberOfVoxelAlongZ; k++)
                {
                    G4int n = Index(i, j, k);
                    // Write only not identically null data lines
                    
                    
                    Analysis->FillNtupleIColumn(1,0, i);
                    Analysis->FillNtupleIColumn(1,1, j);
                    Analysis->FillNtupleIColumn(1,2, k);
                    if (matrix[n])
                    {
                        ofs << G4endl;
                        ofs << i << '\t' << j << '\t' << k << '\t';
                        // Total dose
                        ofs << std::setw(width) << (matrix[n]/massOfVoxel)/doseUnit;
                       
                        
                        Analysis->FillNtupleDColumn(1,3, (matrix[n]/massOfVoxel)/doseUnit);
                        if (secondary)
                        {
                            for (size_t l=0; l < ionStore.size(); l++)
                            {
                                // Fill ASCII file rows
                                ofs << std::setw(width) << ionStore[l].dose[n]/massOfVoxel/doseUnit <<
                                std::setw(width) << ionStore[l].fluence[n];
                                
                                
                                Analysis->FillNtupleIColumn(1,4, ionStore[l].A);
                                Analysis->FillNtupleIColumn(1,5, ionStore[l].Z);
                                
                                Analysis->FillNtupleDColumn(1,6, ionStore[l].dose[n]/massOfVoxel/doseUnit);
                                Analysis->FillNtupleDColumn(1,7, ionStore[l].fluence[n]);
                                Analysis->AddNtupleRow(1);
                                
                                
                                
                            }
                        }
                    }
                }
        ofs.close();
        
        Analysis->Write();
        Analysis->CloseFile();
    }
    
    

}
void HadrontherapyMatrix::Fill(G4int i, G4int j, G4int k, 
                               G4double energyDeposit)
{
    if (matrix)
        matrix[Index(i,j,k)] += energyDeposit;
    
    // Store the energy deposit in the matrix element corresponding 
    // to the phantom voxel  
}




