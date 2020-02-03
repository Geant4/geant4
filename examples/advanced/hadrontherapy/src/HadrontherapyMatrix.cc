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
#include "HadrontherapyAnalysisManager.hh"
#include "HadrontherapyPrimaryGeneratorAction.hh"
#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "G4RunManager.hh"
#include "G4ParticleGun.hh"
<<<<<<< HEAD

// Units definition: CLHEP/Units/SystemOfUnits.h
//
=======
#include "HadrontherapySteppingAction.hh"
#include "HadrontherapyAnalysisFileMessenger.hh"
#include "G4SystemOfUnits.hh"
#include <time.h>

HadrontherapyAnalysis* HadrontherapyAnalysis::instance = 0;
/////////////////////////////////////////////////////////////////////////////

HadrontherapyAnalysis::HadrontherapyAnalysis()
{
    fMess = new HadrontherapyAnalysisFileMessenger(this);
}

/////////////////////////////////////////////////////////////////////////////
HadrontherapyAnalysis::~HadrontherapyAnalysis()
{
    delete fMess;
}

/////////////////////////////////////////////////////////////////////////////
HadrontherapyAnalysis* HadrontherapyAnalysis::GetInstance(){
    
  if (instance == 0) instance = new HadrontherapyAnalysis;
    return instance;
}

>>>>>>> 5baee230e93612916bcea11ebf822756cfa7282c
HadrontherapyMatrix* HadrontherapyMatrix::instance = NULL;
G4bool HadrontherapyMatrix::secondary = false;

// Only return a pointer to matrix
HadrontherapyMatrix* HadrontherapyMatrix::GetInstance() 
{
	return instance;
}
<<<<<<< HEAD
	// This STATIC method delete (!) the old matrix and rewrite a new object returning a pointer to it
	// TODO A check on the parameters is required!
HadrontherapyMatrix* HadrontherapyMatrix::GetInstance(G4int voxelX, G4int voxelY, G4int voxelZ, G4double mass)  
=======

/////////////////////////////////////////////////////////////////////////////
// This STATIC method delete (!) the old matrix and rewrite a new object returning a pointer to it
// TODO A check on the parameters is required!
HadrontherapyMatrix* HadrontherapyMatrix::GetInstance(G4int voxelX, G4int voxelY, G4int voxelZ, G4double mass)
>>>>>>> 5baee230e93612916bcea11ebf822756cfa7282c
{
	if (instance) delete instance;
	instance = new HadrontherapyMatrix(voxelX, voxelY, voxelZ, mass);
	instance -> Initialize(); 
	return instance;
}
<<<<<<< HEAD
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
=======

/////////////////////////////////////////////////////////////////////////////
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
>>>>>>> 5baee230e93612916bcea11ebf822756cfa7282c
}

/////////////////////////////////////////////////////////////////////////////
HadrontherapyMatrix::~HadrontherapyMatrix()
{
    delete[] matrix;
    delete[] hitTrack;
    // free fluences/dose data memory
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
<<<<<<< HEAD
	/////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////
	// Print generated nuclides list
=======

/////////////////////////////////////////////////////////////////////////////
// Print generated nuclides list

/////////////////////////////////////////////////////////////////////////////
>>>>>>> 5baee230e93612916bcea11ebf822756cfa7282c
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
				if (energyDeposit > 0.) ionStore[l].dose[Index(i, j, k)] += energyDeposit;
				
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
		
		// TODO Put some verbosity check
		/*
		 G4cout << "Memory space to store the DOSE/FLUENCE into " <<  
		 numberOfVoxelAlongX*numberOfVoxelAlongY*numberOfVoxelAlongZ << 
		 " voxels has been allocated for the nuclide " << newIon.name << 
		 " (Z = " << Z << ", A = " << A << ")" << G4endl ;
		 */
		return true;
    }
    else // XXX Out of memory! XXX
    {
		return false;
    }
	
}
<<<<<<< HEAD
	
	/////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////
	// Methods to store data to filenames...
	////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////
	//
	// General method to store matrix data to filename
=======

/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
// Methods to store data to filenames...
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
//
// General method to store matrix data to filename
>>>>>>> 5baee230e93612916bcea11ebf822756cfa7282c
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
							// Check for data type: u_int, G4double, XXX 
						if (psize == sizeof(unsigned int))
						{
							unsigned int* pdata = (unsigned int*)data;
							if (pdata[n]) ofs << i << '\t' << j << '\t' <<
								k << '\t' << pdata[n] << G4endl;
						}
						else if (psize == sizeof(G4double))
						{
							G4double* pdata = (G4double*)data;
							if (pdata[n]) ofs << i << '\t' << j << '\t' <<
								k << '\t' << pdata[n] << G4endl;
						}
					}
			ofs.close();
		}
    }
}

<<<<<<< HEAD
	// Store fluence per single ion in multiple files
=======
/////////////////////////////////////////////////////////////////////////////
// Store fluence per single ion in multiple files
>>>>>>> 5baee230e93612916bcea11ebf822756cfa7282c
void HadrontherapyMatrix::StoreFluenceData()
{
    for (size_t i=0; i < ionStore.size(); i++){
		StoreMatrix(ionStore[i].name + "_Fluence.out", ionStore[i].fluence, sizeof(unsigned int));
    }
}
<<<<<<< HEAD
	// Store dose per single ion in multiple files
=======

/////////////////////////////////////////////////////////////////////////////
// Store dose per single ion in multiple files
>>>>>>> 5baee230e93612916bcea11ebf822756cfa7282c
void HadrontherapyMatrix::StoreDoseData()
{
	
    for (size_t i=0; i < ionStore.size(); i++){
		StoreMatrix(ionStore[i].name + "_Dose.out", ionStore[i].dose, sizeof(G4double));
    }
}
<<<<<<< HEAD
	/////////////////////////////////////////////////////////////////////////
	// Store dose for all ions into a single file and into ntuples.
	// Please note that this function is called via messenger commands
	// defined in the HadrontherapyAnalysisFileMessenger.cc class file
=======


////////////////////////////////////////////////////////////////////////
// Store dose into a single file
// or in histograms. Please note that this function is called via
// messenger commands
// defined in the HadrontherapyAnalysisFileMessenger.cc class file
>>>>>>> 5baee230e93612916bcea11ebf822756cfa7282c
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
<<<<<<< HEAD
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

	    /*
	     * PDGencondig
	     */
	    /*
	    ofs << std::setprecision(6) << std::left <<
	    "0\t0\t0\t"; 

	    // Total dose 
	    ofs << std::setw(width) << '0';
	    for (size_t l=0; l < ionStore.size(); l++)
	    {
	    ofs << std::setw(width) << ionStore[l].PDGencoding  <<
	    std::setw(width) << ionStore[l].PDGencoding;
	    }
	    ofs << G4endl;
	    */
	}
	// Write data
	for(G4int i = 0; i < numberOfVoxelAlongX; i++) 
	    for(G4int j = 0; j < numberOfVoxelAlongY; j++) 
		for(G4int k = 0; k < numberOfVoxelAlongZ; k++) 
		{
		    G4int n = Index(i, j, k);
		    // Write only not identically null data lines
		    if (matrix[n])
		    {
			ofs << G4endl;
			ofs << i << '\t' << j << '\t' << k << '\t';
			// Total dose 
			ofs << std::setw(width) << (matrix[n]/massOfVoxel)/doseUnit; 
			if (secondary)
			{
			    for (size_t l=0; l < ionStore.size(); l++)
			    {
				// Fill ASCII file rows
				ofs << std::setw(width) << ionStore[l].dose[n]/massOfVoxel/doseUnit <<
				    std::setw(width) << ionStore[l].fluence[n]; 
			    }
			}
		    }
		}
	ofs.close();
=======
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
        
        // Write data
        for(G4int i = 0; i < numberOfVoxelAlongX; i++)
            for(G4int j = 0; j < numberOfVoxelAlongY; j++)
                for(G4int k = 0; k < numberOfVoxelAlongZ; k++)
                {
                    G4int n = Index(i, j, k);
                    
                    if (matrix[n])
                    {
                        ofs << G4endl;
                        ofs << i << '\t' << j << '\t' << k << '\t';
                        
                        // Total dose
                        ofs << std::setw(width) << (matrix[n]/massOfVoxel)/doseUnit;
                       
                        
                        if (secondary)
                        {
                            for (size_t l=0; l < ionStore.size(); l++)
                            {
                                // Fill ASCII file rows
                                ofs << std::setw(width) << ionStore[l].dose[n]/massOfVoxel/doseUnit <<
                                std::setw(width) << ionStore[l].fluence[n];
                                
                            }
                        }
                    }
                }
        ofs.close();
>>>>>>> 5baee230e93612916bcea11ebf822756cfa7282c
    }
}
/////////////////////////////////////////////////////////////////////////////

#ifdef G4ANALYSIS_USE_ROOT
void HadrontherapyMatrix::StoreDoseFluenceRoot()
{
    HadrontherapyAnalysisManager* analysis = HadrontherapyAnalysisManager::GetInstance();
    if (analysis -> IsTheTFile())
    {
	for(G4int i = 0; i < numberOfVoxelAlongX; i++) 
	    for(G4int j = 0; j < numberOfVoxelAlongY; j++) 
		for(G4int k = 0; k < numberOfVoxelAlongZ; k++) 
		{
		    G4int n = Index(i, j, k);
		    for (size_t l=0; l < ionStore.size(); l++)

		    {
			// Do the same work for .root file: fill dose/fluence ntuple  
			analysis -> FillVoxelFragmentTuple( i, j, k, 
				ionStore[l].A, 
				ionStore[l].Z, 
				ionStore[l].dose[n]/massOfVoxel/doseUnit, 
				ionStore[l].fluence[n] );


		    }
		}
    }
}
<<<<<<< HEAD
#endif

void HadrontherapyMatrix::Fill(G4int i, G4int j, G4int k, 
			       G4double energyDeposit)
=======
//////////////////////////////////////////////////////////////////////////////
void HadrontherapyMatrix::Fill(G4int i, G4int j, G4int k,
                               G4double energyDeposit)
>>>>>>> 5baee230e93612916bcea11ebf822756cfa7282c
{
    if (matrix)
		matrix[Index(i,j,k)] += energyDeposit;
	
    // Store the energy deposit in the matrix element corresponding 
    // to the phantom voxel  
}
void HadrontherapyMatrix::TotalEnergyDeposit()
{
    // Convert energy deposited to dose.
    // Store the information of the matrix in a ntuple and in 
    // a 1D Histogram
#ifdef G4ANALYSIS_USE_ROOT
      HadrontherapyAnalysisManager* analysis = HadrontherapyAnalysisManager::GetInstance();
#endif

    if (matrix)
    {  
	for(G4int i = 0; i < numberOfVoxelAlongX; i++) 
	    for(G4int j = 0; j < numberOfVoxelAlongY; j++) 
		for(G4int k = 0; k < numberOfVoxelAlongZ; k++)
		{
#ifdef G4ANALYSIS_USE_ROOT
		    G4int n = Index(i,j,k);
		    if (analysis -> IsTheTFile() )
		    {
			analysis -> FillEnergyDeposit(i, j, k, matrix[n]/massOfVoxel/doseUnit);
			analysis -> BraggPeak(i, matrix[n]/massOfVoxel/doseUnit);
		    }
#endif
		}
    }
}

