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
// $Id: HadrontherapyLet.cc,v 1.0, Dec 2009

#include "HadrontherapyLet.hh"
#include "HadrontherapyInteractionParameters.hh"
#include "HadrontherapyPrimaryGeneratorAction.hh"
#include "HadrontherapyMatrix.hh"
#include "HadrontherapyDetectorConstruction.hh"
#include "G4RunManager.hh"

#include <cmath>

HadrontherapyLet* HadrontherapyLet::instance = NULL;

HadrontherapyLet* HadrontherapyLet::GetInstance(HadrontherapyDetectorConstruction *pDet)  
{
    if (instance) delete instance;
    instance = new HadrontherapyLet(pDet);
    return instance;
}

HadrontherapyLet* HadrontherapyLet::GetInstance()  
{
    return instance;
}

HadrontherapyLet::HadrontherapyLet(HadrontherapyDetectorConstruction* pDet)
   :pParam(0)
{
    pParam = new HadrontherapyInteractionParameters(false); // no messenger
    //  letMessenger = new HadrontherapyLetMessenger(this);
    matrix = HadrontherapyMatrix::GetInstance();
    if (!matrix) G4Exception("HadrontherapyMatrix not found. Firstly instance it.");

    nVoxels = matrix -> GetNvoxel();
    numberOfVoxelAlongX = matrix -> GetNumberOfVoxelAlongX();
    numberOfVoxelAlongY = matrix -> GetNumberOfVoxelAlongY();
    numberOfVoxelAlongZ = matrix -> GetNumberOfVoxelAlongZ();

    G4RunManager *runManager = G4RunManager::GetRunManager();
    pPGA = (HadrontherapyPrimaryGeneratorAction*)runManager -> GetUserPrimaryGeneratorAction();
    // Pointer to the detector material
    detectorMat = pDet -> GetDetectorLogicalVolume() -> GetMaterial();
    density = detectorMat -> GetDensity();

}
HadrontherapyLet::~HadrontherapyLet()
{
    delete pParam;
    Clear();
}

// Fill energy spectrum for every voxel (local energy spectrum)
void HadrontherapyLet::Initialize()
{
    primaryEnergy = pPGA -> GetmeanKineticEnergy();
    energyLimit =   trunc((primaryEnergy /10.)+1.)*10.;// round toward zero
    binWidth =      0.25*MeV;
    nBins = 	    (G4int)ceil(energyLimit/binWidth); //round up toward nearest integer
    // Clear data, if any 
    Clear();
}
void HadrontherapyLet::Clear()
{
    for (size_t i=0; i < ionLetStore.size(); i++)
    {
	for(G4int v=0; v < nVoxels; v++) 
	    if(ionLetStore[i].spectrum[v]) delete[] ionLetStore[i].spectrum[v]; 
	delete[] ionLetStore[i].spectrum;
    }
    ionLetStore.clear();

}
void  HadrontherapyLet::FillEnergySpectrum(G4int trackID,
					   G4ParticleDefinition* particleDef, 
					   G4double kinEnergy, 
					   G4int i, G4int j, G4int k) 
{
    // First step energy
    if (kinEnergy<=0) return;
    G4int Z = particleDef -> GetAtomicNumber();
    G4int A = particleDef -> GetAtomicMass();

    G4String fullName = particleDef -> GetParticleName();
    G4String name = fullName.substr (0, fullName.find("[") ); // cut excitation energy  
    // Is it a primary particle?
    //if (trackID == 1) name +="_1"; 

    G4int voxel = matrix -> Index(i,j,k);
    G4int enBin = lround(trunc(kinEnergy/binWidth)); 
    // bins are [n.binWidth, n.binWidth + binWidth), n natural
    // for example for 0.25 binWidth we have [0, 0.25) [0.25, 0.5) [0.5, 0.75) [0.75, 1) [1, 1.25) ...
    
    // Search for already allocated data...
    size_t l;
    for (l=0; l < ionLetStore.size(); l++) 
    {
	if (ionLetStore[l].name == name) 
	    if ( trackID ==1 && ionLetStore[l].isPrimary || trackID !=1 && !ionLetStore[l].isPrimary)
		break;
    }
    //for (vector<ionLet>::const_iterator iter = ionLetStore.begin(); iter!=ionLetStore.end(); ++iter)
    //if ((*iter).name == name) break; 
    
    if (l == ionLetStore.size()) // Just another type of ion/particle for our store...
    {
    ionLet ion =
	{
	    (trackID == 1) ? true:false, // is it the primary particle? 
	    fullName,
	    name,
	    Z,
	    A,
	    new G4double[nBins],   // Stopping Powers table
	    new G4int*  [nVoxels], // Array of pointers to Energy Spectrum
	    new G4double[nVoxels], // Let_T
	    new G4double[nVoxels]  // Let_D
	};

	// Get stopping powers table (keV/um)
	G4int bin = 0;
	for(G4double E = binWidth/2; E < energyLimit ; E += binWidth )
	{
	    ion.stop[bin++] = pParam -> GetStopping(E, particleDef, detectorMat)*(keV/um); // Total linear stopping power (keV/um)
	    //G4cout << E/MeV << '\t' << ion.stop[bin-1] << '\t' << bin-1 << '\t' << nBins << '\n';
	}
	// Clear array of pointer to histograms
	for(G4int v=0; v < nVoxels; v++) ion.spectrum[v] = NULL;
	// Initialize let
	for(G4int v=0; v < nVoxels; v++) ion.letT[v] = ion.letD[v] = 0.;
	ionLetStore.push_back(ion);
	//G4cout << "Allocated LET data for " << ion.name << G4endl;

    }

    if (!ionLetStore[l].spectrum[voxel]) 
    {
	ionLetStore[l].spectrum[voxel] = new G4int[nBins];// allocate new histogram for every hit voxel!
	for(G4int bin=0; bin < nBins; bin++) ionLetStore[l].spectrum[voxel][bin] = 0; // clear it
    }

    ionLetStore[l].spectrum[voxel][enBin]++; // fill spectrum
}

// LET calculation
// Must be issued at endOfRunAction!
void HadrontherapyLet::LetOutput()
{
    std::sort(ionLetStore.begin(), ionLetStore.end());
    for (size_t l=0; l < ionLetStore.size(); l++)
    {
	for(G4int v=0; v < nVoxels; v++)  
	{
	    // Numeric calculation to get LET Track & LET Dose (keV/um)   
	    // Med. Phys. 30(5), May 2003. Equations (13) and (14)
	    if ( ionLetStore[l].spectrum[v] )
	    {
		nT = dT = nD = dD = 0.;
		for(G4int bin=0; bin < nBins; bin++) // histogram bin index
		{
		    // numerator and denominator for Let_Track(nT,dT), Let_Dose(nD,dD)
		    nT += ionLetStore[l].spectrum[v][bin]*ionLetStore[l].stop[bin];
		    dT += ionLetStore[l].spectrum[v][bin];

		    nD += ionLetStore[l].spectrum[v][bin]*(ionLetStore[l].stop[bin]*ionLetStore[l].stop[bin]);
		    dD += ionLetStore[l].spectrum[v][bin]*ionLetStore[l].stop[bin];

		} 
		//G4cout << "LetT " << nT/dT << " LetD " << nD/dD << '\n';  
		ionLetStore[l].letT[v] = nT/dT; //(keV/um)
		ionLetStore[l].letD[v] = nD/dD; //(keV/um)
	    }
	}
    }
}
void HadrontherapyLet::StoreData(G4String filename)
{
#define width 15L
    if(ionLetStore.size())
    {
	if (ofs)
	{
	    ofs.open(filename, std::ios::out);
	    // Write the voxels index and the list of particles/ions 
	    ofs << std::setprecision(6) << std::left <<
		"i\tj\tk\t"; 
	    for (size_t l=0; l < ionLetStore.size(); l++)
	    {
		G4String a = (ionLetStore[l].isPrimary) ? "_1":"";
		ofs << std::setw(width) << ionLetStore[l].name + "_lT" + a <<
		       std::setw(width) << ionLetStore[l].name + "_lD" + a;
	    }
	    ofs << G4endl;
	    ofs << std::setfill('_');
	    for (size_t l=0; l < 2*ionLetStore.size(); l++)
	    {
		ofs << std::setw(width) <<  "_";
	    }

	    ofs << std::setfill(' ');
	    // Write data
	    for(G4int i = 0; i < numberOfVoxelAlongX; i++) 
		for(G4int j = 0; j < numberOfVoxelAlongY; j++) 
		    for(G4int k = 0; k < numberOfVoxelAlongZ; k++) 
		    {
			G4int v = matrix -> Index(i, j, k);
			for (size_t l=0; l < ionLetStore.size(); l++)
			{
			    // Write only not identically null data lines
			    if(ionLetStore[l].letT[v] || ionLetStore[l].letD[v])
			    {
				ofs << G4endl;
				ofs << i << '\t' << j << '\t' << k << '\t';
				for (size_t l=0; l < ionLetStore.size(); l++)
				{
				    ofs << std::setw(width) << ionLetStore[l].letT[v] <<
					std::setw(width) << ionLetStore[l].letD[v]; 
				}
				break;
			    }
			}
		    }
	    ofs.close();
	}
    }
}


