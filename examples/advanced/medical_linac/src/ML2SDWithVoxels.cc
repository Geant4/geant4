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
//	^Claudio Andenna  claudio.andenna@ispesl.it, claudio.andenna@iss.infn.it
//      *Barbara Caccia barbara.caccia@iss.it
//      with the support of Pablo Cirrone (LNS, INFN Catania Italy)
//	with the contribute of Alessandro Occhigrossi*
//
// ^INAIL DIPIA - ex ISPESL and INFN Roma, gruppo collegato Sanità, Italy
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

CML2SDWithVoxels::CML2SDWithVoxels(G4String name, G4int saving_in_ROG_Voxels_every_events, G4int seed, G4String ROGOutFile, G4bool bSaveROG, G4ThreeVector centre, G4ThreeVector halfSize, G4int NumberOfVoxelsAlongX, G4int NumberOfVoxelsAlongY, G4int NumberOfVoxelsAlongZ)
: G4VSensitiveDetector(name),voxelsSum(0), voxelsSingle(0)
{
	this->saving_in_ROG_Voxels_every_events=saving_in_ROG_Voxels_every_events;
	this->bSaveROG=bSaveROG;
	this->bActive=true;
	this->nParticle=0;
	this->nParticleValatile=0;
	this->nTotalEvents=0;
	this->nSingleTotalEvents=0;
	this->density=1.;
	this->voxelVolume=0.;
	this->voxelMass=0.;
	this->nRecycling=1;

	G4String seedName;
	char a[10];
	sprintf(a,"%d", seed);
	seedName=(G4String)a;

	this->fullOutFileData=ROGOutFile+"_"+seedName+".txt";
	this->fullOutFileDataSingle="";
	if(this->bSaveROG)
	{
		this->centre=centre;
		this->halfSize=halfSize;
		this->NumberOfVoxelsAlongX=NumberOfVoxelsAlongX;
		this->NumberOfVoxelsAlongY=NumberOfVoxelsAlongY;
		this->NumberOfVoxelsAlongZ=NumberOfVoxelsAlongZ;
		this->halfXVoxelDimensionX=this->halfSize.getX()/this->NumberOfVoxelsAlongX;
		this->halfXVoxelDimensionY=this->halfSize.getY()/this->NumberOfVoxelsAlongY;
		this->halfXVoxelDimensionZ=this->halfSize.getZ()/this->NumberOfVoxelsAlongZ;

		this->voxelVolume=this->halfXVoxelDimensionX*this->halfXVoxelDimensionY*this->halfXVoxelDimensionZ*8.;

// voxels to store and save the sum of the geometry configurations
		this->voxelsSum=new Svoxel**[this->NumberOfVoxelsAlongX];
		for (int ix=0; ix< this->NumberOfVoxelsAlongX; ix++)
		{
			this->voxelsSum[ix]=new Svoxel*[this->NumberOfVoxelsAlongY];
			for (int iy=0; iy< this->NumberOfVoxelsAlongY; iy++)
			{
				this->voxelsSum[ix][iy]=new Svoxel[this->NumberOfVoxelsAlongZ];
				for (int iz=0; iz< this->NumberOfVoxelsAlongZ; iz++)
				{
					this->voxelsSum[ix][iy][iz].volumeId=-1;
					this->voxelsSum[ix][iy][iz].depEnergy=0.;
					this->voxelsSum[ix][iy][iz].depEnergy2=0.;
					this->voxelsSum[ix][iy][iz].depEnergyNorm=0.;
					this->voxelsSum[ix][iy][iz].depEnergyNormError=0.;
					this->voxelsSum[ix][iy][iz].expDose=0.;
					this->voxelsSum[ix][iy][iz].halfSize.set(this->halfXVoxelDimensionX, this->halfXVoxelDimensionY, this->halfXVoxelDimensionZ);
					this->voxelsSum[ix][iy][iz].pos.set(2.*(ix)*this->halfXVoxelDimensionX  -this->halfSize.getX()+this->halfXVoxelDimensionX + this->centre.getX(), 
			2.*(iy)*this->halfXVoxelDimensionY  -this->halfSize.getY()+this->halfXVoxelDimensionY + this->centre.getY(), 
			2.*(iz)*this->halfXVoxelDimensionZ  -this->halfSize.getZ()+this->halfXVoxelDimensionZ + this->centre.getZ());
					this->voxelsSum[ix][iy][iz].nEvents=0;
				}
			}
		}


// voxels to store and save the single geometry configuration
		this->voxelsSingle=new Svoxel**[this->NumberOfVoxelsAlongX];
		for (int ix=0; ix< this->NumberOfVoxelsAlongX; ix++)
		{
			this->voxelsSingle[ix]=new Svoxel*[this->NumberOfVoxelsAlongY];
			for (int iy=0; iy< this->NumberOfVoxelsAlongY; iy++)
			{
				this->voxelsSingle[ix][iy]=new Svoxel[this->NumberOfVoxelsAlongZ];
				for (int iz=0; iz< this->NumberOfVoxelsAlongZ; iz++)
				{
					this->voxelsSingle[ix][iy][iz].volumeId=-1;
					this->voxelsSingle[ix][iy][iz].depEnergy=0.;
					this->voxelsSingle[ix][iy][iz].depEnergy2=0.;
					this->voxelsSingle[ix][iy][iz].depEnergyNorm=0.;
					this->voxelsSingle[ix][iy][iz].depEnergyNormError=0.;
					this->voxelsSingle[ix][iy][iz].expDose=0.;
					this->voxelsSingle[ix][iy][iz].halfSize.set(this->halfXVoxelDimensionX, this->halfXVoxelDimensionY, this->halfXVoxelDimensionZ);
					this->voxelsSingle[ix][iy][iz].pos.set(2.*(ix)*this->halfXVoxelDimensionX  -this->halfSize.getX()+this->halfXVoxelDimensionX + this->centre.getX(), 
			2.*(iy)*this->halfXVoxelDimensionY  -this->halfSize.getY()+this->halfXVoxelDimensionY + this->centre.getY(), 
			2.*(iz)*this->halfXVoxelDimensionZ  -this->halfSize.getZ()+this->halfXVoxelDimensionZ + this->centre.getZ());
					this->voxelsSingle[ix][iy][iz].nEvents=0;
				}
			}
		}
	}
}

CML2SDWithVoxels::~CML2SDWithVoxels()
{
	if(this->bSaveROG)
	{
		delete [] this->voxelsSum;
		delete this->voxelsSum;

		delete [] this->voxelsSingle;
		delete this->voxelsSingle;
	}
}
void CML2SDWithVoxels::resetVoxelsSingle()
{
	for (int ix=0; ix< this->NumberOfVoxelsAlongX; ix++)
	{
		for (int iy=0; iy< this->NumberOfVoxelsAlongY; iy++)
		{
			for (int iz=0; iz< this->NumberOfVoxelsAlongZ; iz++)
			{
				this->voxelsSingle[ix][iy][iz].volumeId=-1;
				this->voxelsSingle[ix][iy][iz].depEnergy=0.;
				this->voxelsSingle[ix][iy][iz].depEnergy2=0.;
				this->voxelsSingle[ix][iy][iz].depEnergyNorm=0.;
				this->voxelsSingle[ix][iy][iz].depEnergyNormError=0.;
				this->voxelsSingle[ix][iy][iz].expDose=0.;
				this->voxelsSingle[ix][iy][iz].halfSize.set(this->halfXVoxelDimensionX, this->halfXVoxelDimensionY, this->halfXVoxelDimensionZ);
				this->voxelsSingle[ix][iy][iz].pos.set(2.*(ix)*this->halfXVoxelDimensionX  -this->halfSize.getX()+this->halfXVoxelDimensionX + this->centre.getX(), 
		2.*(iy)*this->halfXVoxelDimensionY  -this->halfSize.getY()+this->halfXVoxelDimensionY + this->centre.getY(), 
		2.*(iz)*this->halfXVoxelDimensionZ  -this->halfSize.getZ()+this->halfXVoxelDimensionZ + this->centre.getZ());
				this->voxelsSingle[ix][iy][iz].nEvents=0;
			}
		}
	}
	this->nSingleTotalEvents=0;
}
G4bool CML2SDWithVoxels::ProcessHits(G4Step *aStep, G4TouchableHistory *ROHist)
{
	
	if (this->bActive)
	{
		G4double energyDep = aStep->GetTotalEnergyDeposit();
		if (this->bSaveROG && energyDep>0.) 
		{
			G4int ix, iy, iz;
			G4String volumeName;

			ix=ROHist->GetReplicaNumber(2);
			iy=ROHist->GetReplicaNumber(0);
			iz=ROHist->GetReplicaNumber(1);

			this->density=aStep->GetPreStepPoint()->GetPhysicalVolume()->GetLogicalVolume()->GetMaterial()->GetDensity();

			this->voxelMass=this->voxelVolume*this->density;
			energyDep/=this->voxelMass*this->nRecycling;
			this->voxelsSum[ix][iy][iz].volumeId=this->getIdFromVolumeName(aStep->GetPreStepPoint()->GetPhysicalVolume()->GetLogicalVolume()->GetMaterial()->GetName());
			this->voxelsSum[ix][iy][iz].depEnergy+=energyDep;
			this->voxelsSum[ix][iy][iz].depEnergy2+=energyDep*energyDep;
			this->voxelsSum[ix][iy][iz].nEvents++;

			this->voxelsSingle[ix][iy][iz].volumeId=this->getIdFromVolumeName(aStep->GetPreStepPoint()->GetPhysicalVolume()->GetLogicalVolume()->GetMaterial()->GetName());
			this->voxelsSingle[ix][iy][iz].depEnergy+=energyDep;
			this->voxelsSingle[ix][iy][iz].depEnergy2+=energyDep*energyDep;
			this->voxelsSingle[ix][iy][iz].nEvents++;
			this->nTotalEvents++;
			this->nSingleTotalEvents++;
			if (this->nTotalEvents%this->saving_in_ROG_Voxels_every_events==0 && this->nTotalEvents>0)
			{
				this->save();
			}
		}
	}
	return true;
}
G4int CML2SDWithVoxels::getIdFromVolumeName(G4String name)
{
	for (int i=0; i<(int)this->volumeNameIdLink.size(); i++)
	{
		if (this->volumeNameIdLink[i].volumeName==name)
		{
			return this->volumeNameIdLink[i].volumeId;
			break;
		}
	}
	return -1;
}
void CML2SDWithVoxels::save()
{
std::cout<< "n. of events collected in the whole ROG phantom for all geometries: "<< this->nTotalEvents<< G4endl;
std::cout<< "n. of events collected in the whole ROG phantom for the current geometry: "<< this->nSingleTotalEvents<< G4endl;
	if (this->nTotalEvents>0)
	{this->saveData(this->fullOutFileData, this->voxelsSum);}
	if (this->nSingleTotalEvents>0)
	{this->saveData(this->fullOutFileDataSingle, this->voxelsSingle);}
}
void CML2SDWithVoxels::saveData(G4String Filename, Svoxel ***voxels)
{
	std::ofstream out;
	out.open(Filename, std::ios::out);
	out << "Sensitive Detector-Voxels. Total number of events, [mm]->centreX centreY centreZ HalfSizeX HalfSizeY HalfSizeZ minX maxX, minY maxY, minZ maxZ, Dx, Dy, Dz, nX, nY, nZ: \n";
	out <<this->nTotalEvents<<'\t';
	out <<this->centre.getX()/mm << '\t' << this->centre.getY()/mm<< '\t'<< this->centre.getZ()/mm<<'\t';
	out <<this->halfSize.getX()/mm << '\t' << this->halfSize.getY()/mm <<'\t'<< this->halfSize.getZ()/mm<<'\t';
	out <<(this->centre.getX()-this->halfSize.getX())/mm<<'\t'<<(this->centre.getX()+this->halfSize.getX())/mm<<'\t';
	out <<(this->centre.getY()-this->halfSize.getY())/mm<<'\t'<<(this->centre.getY()+this->halfSize.getY())/mm<<'\t';
	out <<(this->centre.getZ()-this->halfSize.getZ())/mm<<'\t'<<(this->centre.getZ()+this->halfSize.getZ())/mm<<'\t';
	out <<this->halfXVoxelDimensionX/mm<<'\t'<<this->halfXVoxelDimensionY/mm<<'\t'<<this->halfXVoxelDimensionZ/mm<<'\t';
	out <<this->NumberOfVoxelsAlongX <<'\t'<<this->NumberOfVoxelsAlongY <<'\t'<<this->NumberOfVoxelsAlongZ <<'\n';
	out << "Number of physical volumes: "<< this->volumeNameIdLink.size() << '\n';
	for (int i=0; i<(int)this->volumeNameIdLink.size(); i++)
	{
		out << this->volumeNameIdLink[i].volumeName  <<'\t'<< this->volumeNameIdLink[i].volumeId << G4endl;
	}


	out << "Phys Volume x [mm], y [mm], z [mm], ix, iy, iz, Dose [Gy], Dose2 [Gy^2], nEvents" << G4endl;

	for (int ix=0; ix< this->NumberOfVoxelsAlongX; ix++)
	{
		for (int iy=0; iy< this->NumberOfVoxelsAlongY; iy++)
		{
			for (int iz=0; iz< this->NumberOfVoxelsAlongZ; iz++)
			{
				if (voxels[ix][iy][iz].nEvents>0)
				{
					out << voxels[ix][iy][iz].volumeId << '\t';
					out << voxels[ix][iy][iz].pos.getX()/mm << '\t';
					out << voxels[ix][iy][iz].pos.getY()/mm << '\t';
					out << voxels[ix][iy][iz].pos.getZ()/mm << '\t';
					out << ix << '\t';
					out << iy << '\t';
					out << iz << '\t';
					out << voxels[ix][iy][iz].depEnergy/(joule/kg) << '\t';
					out << voxels[ix][iy][iz].depEnergy2/((joule/kg)*(joule/kg)) << '\t';
					out << voxels[ix][iy][iz].nEvents << G4endl;
				}
			}
		}
	}
	out.close();

}
void CML2SDWithVoxels::setFullOutFileDataSingle(G4String val)
{
	unsigned int ind = this->fullOutFileData.find(".txt");
	G4String onlyName=this->fullOutFileData.substr( 0, ind);
	if (val=="")
	{
		static unsigned int indGeom=0;
		char cT[5];
		sprintf(cT,"%d",indGeom);
		this->fullOutFileDataSingle=onlyName+"Single_"+G4String(cT)+".txt";
		indGeom++;
	}
	else
	{
		this->fullOutFileDataSingle=onlyName+val+".txt";
	}
}

