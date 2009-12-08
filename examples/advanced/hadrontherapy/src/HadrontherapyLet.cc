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
// $Id: HadrontherapyLet.cc,v 1.0, May 2007;
// 

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
    pParam = new HadrontherapyInteractionParameters(false);
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
	    if(ionLetStore[i].pSpectrum[v]) delete[] ionLetStore[i].pSpectrum[v]; 
	delete[] ionLetStore[i].pSpectrum;
    }
    ionLetStore.clear();

}
void  HadrontherapyLet::FillEnergySpectrum(G4ParticleDefinition* particleDef, 
					   G4double kinEnergy, 
					   G4int i, G4int j, G4int k) 
{
    // First step energy
    if (kinEnergy<=0) return;
    //G4int Z = particleDef-> GetAtomicNumber();
    //G4int A = particleDef-> GetAtomicMass();
    G4String fullName = particleDef -> GetParticleName();
    G4String name = fullName.substr (0, fullName.find("[") ); // cut excitation energy  
    G4int voxel = matrix -> Index(i,j,k);
    G4int enBin = lround(trunc(kinEnergy/binWidth)); 
    // bins are [n.binWidth, n.binWidth + binWidth), n natural
    // for example for 0.25 binWidth we have [0, 0.25) [0.25, 0.5) [0.5, 0.75) [0.75, 1) [1, 1.25) ...
    //if (i > 219)
    //G4cout << "Slab= " << i << ", E= " << kinEnergy << ", enBin= " << enBin << G4endl;
    // Search for already allocated data...
    for (size_t l=0; l < ionLetStore.size(); l++)
    {
	// Search for particle/ion name...
	if (ionLetStore[l].name == name) 
	{

	    if (!ionLetStore[l].pSpectrum[voxel]) 
	    {
		ionLetStore[l].pSpectrum[voxel] = new G4int[nBins];// allocate new histogram for every hit voxel!
		for(G4int bin=0; bin < nBins; bin++) ionLetStore[l].pSpectrum[voxel][bin] = 0; // clear it
	    }
	    ionLetStore[l].pSpectrum[voxel][enBin]++; // fill spectrum
	    return ; // exit
	}
    }
    // Just another type of ion/particle for our store...
    ionLet ion =
    {
	fullName,
	name,
	new G4double[nBins],   // Stopping Powers table
	new G4int*  [nVoxels], // Energy Spectrum
	new G4double[nVoxels], // Let_T
	new G4double[nVoxels]  // Let_D
    };

    // Stopping powers table
    G4int bin = 0;
    for(G4double E = binWidth/2; E < energyLimit ; E += binWidth )
    {
	ion.stop[bin++] = pParam -> GetStopping(E, particleDef, detectorMat)*(keV/um);
	G4cout << E/MeV << '\t' << ion.stop[bin-1] << '\t' << bin-1 << '\t' << nBins << '\n';
    }
    // Clear array of pointer to histograms
    for(G4int v=0; v < nVoxels; v++) ion.pSpectrum[v] = 0;
    ion.pSpectrum[voxel] = new G4int[nBins];// allocate new histogram
    for(int bin=0; bin<nBins; bin++) ion.pSpectrum[voxel][bin] = 0; // clear it
    
    ion.pSpectrum[voxel][enBin]++; // fill spectrum
    // Initialize let
    for(G4int v=0; v < nVoxels; v++) 
    {
	ion.letT[v] = 0.; 
	ion.letD[v] = 0.;
    }
    G4cout << "Allocated LET data for " << ion.fullName << '\n';
    ionLetStore.push_back(ion);
}

// Issued at endOfRunAction!
// LET calculation
void HadrontherapyLet::LetOutput()
{
    for (size_t l=0; l < ionLetStore.size(); l++)
    {
	for(G4int v=0; v < nVoxels; v++) // indice sui voxels in cui calcolare il LET
	{
	    // Take LET for voxel v  
	    // Med. Phys. 30(5), May 2003. Equations (13) and (14)
	    // only if a Spectrum exists
	    if ( ionLetStore[l].pSpectrum[v] )
	    {
		nT = dT = nD = dD = 0.;
		for(G4int bin=0; bin < nBins; bin++) // m histogram bin index
		{
		    // numeratore e denominatore per Let_Track(nT,dT) e Let_Dose(nD,dD)
		    nT += ionLetStore[l].pSpectrum[v][bin]*ionLetStore[l].stop[bin];
		    dT += ionLetStore[l].pSpectrum[v][bin];

		    nD += ionLetStore[l].pSpectrum[v][bin]*(ionLetStore[l].stop[bin]*ionLetStore[l].stop[bin]);
		    dD += ionLetStore[l].pSpectrum[v][bin]*ionLetStore[l].stop[bin];

		} 
		//G4cout << "LetT " << nT/dT << " LetD " << nD/dD << '\n';  
		ionLetStore[l].letT[v] = (nT/dT);
		ionLetStore[l].letD[v] = (nD/dD);
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
		ofs << std::setw(width) << ionLetStore[l].name + "_LetT" <<
		    std::setw(width) << ionLetStore[l].name + "_LetD";
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


