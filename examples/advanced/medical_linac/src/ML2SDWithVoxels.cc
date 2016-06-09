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
// The code was written by :
//	^Claudio Andenna claudio.andenna@iss.infn.it, claudio.andenna@ispesl.it
//      *Barbara Caccia barbara.caccia@iss.it
//      with the support of Pablo Cirrone (LNS, INFN Catania Italy)
//
// ^ISPESL and INFN Roma, gruppo collegato Sanità, Italy
// *Istituto Superiore di Sanità and INFN Roma, gruppo collegato Sanità, Italy
//  Viale Regina Elena 299, 00161 Roma (Italy)
//  tel (39) 06 49902246
//  fax (39) 06 49387075
//
// more information:
// http://g4advancedexamples.lngs.infn.it/Examples/medical-linac
//
//*******************************************************//


#include "ML2SDWithVoxels.hh"
#include "ML2ExpVoxels.hh"


// *************************************************************************************

CML2SDWithVoxels::CML2SDWithVoxels(G4String name, G4int saving_in_ROG_Voxels_every_events, G4int seed, G4String ROGOutFile, G4bool bSaveROG, G4ThreeVector halfSize, G4int NumberOfVoxelsAlongX, G4int NumberOfVoxelsAlongY, G4int NumberOfVoxelsAlongZ)
: G4VSensitiveDetector(name),voxels(0)
{
	this->saving_in_ROG_Voxels_every_events=saving_in_ROG_Voxels_every_events;
	this->bSaveROG=bSaveROG;
	this->bActive=true;
	this->nParticle=0;
	this->nParticleValatile=0;
	this->nTotalEvents=0;
	this->density=1.;
	this->voxelVolume=0.;
	this->voxelMass=0.;
	G4String seedName;
	char a[10];
	sprintf(a,"%d", seed);
	seedName=(G4String)a;


	this->fullOutFileData=ROGOutFile+"_"+seedName+".txt";
#ifdef ML2FILEOUT
	char *MyDirOut=new char[1000];
	MyDirOut=getenv("ML2FILEOUT");
	G4String myDirOut=(G4String) MyDirOut;
	this->fullOutFileData=myDirOut+"/"+ROGOutFile+"_"+seedName+".txt";
#endif

	if(this->bSaveROG)
	{
		this->halfSize=halfSize;
		this->NumberOfVoxelsAlongX=NumberOfVoxelsAlongX;
		this->NumberOfVoxelsAlongY=NumberOfVoxelsAlongY;
		this->NumberOfVoxelsAlongZ=NumberOfVoxelsAlongZ;
		this->halfXVoxelDimensionX=this->halfSize.getX()/this->NumberOfVoxelsAlongX;
		this->halfXVoxelDimensionY=this->halfSize.getY()/this->NumberOfVoxelsAlongY;
		this->halfXVoxelDimensionZ=this->halfSize.getZ()/this->NumberOfVoxelsAlongZ;

		this->voxelVolume=this->halfXVoxelDimensionX*this->halfXVoxelDimensionY*this->halfXVoxelDimensionZ;

		this->voxels=new Svoxel**[this->NumberOfVoxelsAlongX];
		for (int ix=0; ix< this->NumberOfVoxelsAlongX; ix++)
		{
			this->voxels[ix]=new Svoxel*[this->NumberOfVoxelsAlongY];
			for (int iy=0; iy< this->NumberOfVoxelsAlongY; iy++)
			{
				this->voxels[ix][iy]=new Svoxel[this->NumberOfVoxelsAlongZ];
				for (int iz=0; iz< this->NumberOfVoxelsAlongZ; iz++)
				{
					this->voxels[ix][iy][iz].volumeName="noData";
					this->voxels[ix][iy][iz].depEnergy=0.;
					this->voxels[ix][iy][iz].depEnergy2=0.;
					this->voxels[ix][iy][iz].depEnergyNorm=0.;
					this->voxels[ix][iy][iz].depEnergyNormError=0.;
					this->voxels[ix][iy][iz].expDose=0.;
					this->voxels[ix][iy][iz].halfSize.set(this->halfXVoxelDimensionX, this->halfXVoxelDimensionY, this->halfXVoxelDimensionZ);
					this->voxels[ix][iy][iz].pos.set(2.*(ix)*this->halfXVoxelDimensionX  -this->halfSize.getX()+this->halfXVoxelDimensionX, 
			2.*(iy)*this->halfXVoxelDimensionY  -this->halfSize.getY()+this->halfXVoxelDimensionY, 
			2.*(iz)*this->halfXVoxelDimensionZ  -this->halfSize.getZ()+this->halfXVoxelDimensionZ);
					this->voxels[ix][iy][iz].nEvents=0;
				}
			}
		}
	}
}

CML2SDWithVoxels::~CML2SDWithVoxels()
{
	if(this->bSaveROG)
	{
		delete [] this->voxels;
		delete this->voxels;
	}
}
G4bool CML2SDWithVoxels::ProcessHits(G4Step *aStep, G4TouchableHistory *ROHist)
{
	
	if (this->bActive)
	{
		G4double energyDep = aStep->GetTotalEnergyDeposit();
		if (this->bSaveROG && energyDep>0.) 
		{
			G4int ix, iy, iz;

			ix=ROHist->GetReplicaNumber(2);
			iy=ROHist->GetReplicaNumber(0);
			iz=ROHist->GetReplicaNumber(1);

			this->density=aStep->GetPreStepPoint()->GetPhysicalVolume()->GetLogicalVolume()->GetMaterial()->GetDensity();
			this->voxelMass=this->voxelVolume*this->density;
			energyDep/=this->voxelMass;
			this->voxels[ix][iy][iz].volumeName=aStep->GetPreStepPoint()->GetPhysicalVolume()->GetName();
			this->voxels[ix][iy][iz].depEnergy+=energyDep;
			this->voxels[ix][iy][iz].depEnergy2+=energyDep*energyDep;
			this->voxels[ix][iy][iz].nEvents++;
			this->nTotalEvents++;
			if (this->nTotalEvents%this->saving_in_ROG_Voxels_every_events==0 && this->nTotalEvents>0)
			{
				this->saveDataInVoxels();
			}
		}
	}
	return true;
}

void CML2SDWithVoxels::saveHeaderDataInVoxels()
{
	std::ofstream out;
	out.open(this->fullOutFileData, std::ios::out);
	out << "Sensitive Detector-Voxels. Total number of events: "<<this->nTotalEvents<<G4endl;
	out << "Phys Volume x [mm], y [mm], z [mm], ix, iy, iz, Dose [Gy], Dose2 [Gy^2], nEvents" << G4endl;
	out.close();
}
void CML2SDWithVoxels::saveDataInVoxels()
{
	this->saveHeaderDataInVoxels();
	std::ofstream out;
	out.open(this->fullOutFileData, std::ios::app);
	for (int ix=0; ix< this->NumberOfVoxelsAlongX; ix++)
	{
		for (int iy=0; iy< this->NumberOfVoxelsAlongY; iy++)
		{
			for (int iz=0; iz< this->NumberOfVoxelsAlongZ; iz++)
			{
				if (this->voxels[ix][iy][iz].nEvents>0)
				{
					out << this->voxels[ix][iy][iz].volumeName << '\t';
					out << this->voxels[ix][iy][iz].pos.getX()/mm << '\t';
					out << this->voxels[ix][iy][iz].pos.getY()/mm << '\t';
					out << this->voxels[ix][iy][iz].pos.getZ()/mm << '\t';
					out << ix << '\t';
					out << iy << '\t';
					out << iz << '\t';
					out << this->voxels[ix][iy][iz].depEnergy/(joule/kg) << '\t';
					out << this->voxels[ix][iy][iz].depEnergy2/((joule/kg)*(joule/kg)) << '\t';
					out << this->voxels[ix][iy][iz].nEvents << G4endl;
				}
			}
		}
	}
	out.close();
}
