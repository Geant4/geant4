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
#include "G4SystemOfUnits.hh"

// *************************************************************************************

CML2SDWithVoxels::CML2SDWithVoxels(G4String name, G4int voxSave,
                                   G4int seed, G4String ROGOutFile,
                                   G4bool bROG, G4ThreeVector ctr, G4ThreeVector hSiz,
                                   G4int NumVX, G4int NumVY, G4int NumVZ)
: G4VSensitiveDetector(name),voxelsSum(0), voxelsSingle(0)
{
	saving_in_ROG_Voxels_every_events=voxSave;
	bSaveROG=bROG;
	bActive=true;
	nParticle=0;
	nParticleValatile=0;
	nTotalEvents=0;
	nSingleTotalEvents=0;
	density=1.;
	voxelVolume=0.;
	voxelMass=0.;
	nRecycling=1;

	G4String seedName;
	char a[10];
	sprintf(a,"%d", seed);
	seedName=(G4String)a;

	fullOutFileData=ROGOutFile+"_"+seedName+".txt";
	fullOutFileDataSingle="";
	if(bSaveROG)
	{
		centre=ctr;
		halfSize=hSiz;
		NumberOfVoxelsAlongX=NumVX;
		NumberOfVoxelsAlongY=NumVY;
		NumberOfVoxelsAlongZ=NumVZ;
		halfXVoxelDimensionX=halfSize.getX()/NumberOfVoxelsAlongX;
		halfXVoxelDimensionY=halfSize.getY()/NumberOfVoxelsAlongY;
		halfXVoxelDimensionZ=halfSize.getZ()/NumberOfVoxelsAlongZ;

		voxelVolume=halfXVoxelDimensionX*halfXVoxelDimensionY*halfXVoxelDimensionZ*8.;

// voxels to store and save the sum of the geometry configurations
		voxelsSum=new Svoxel**[NumberOfVoxelsAlongX];
		for (int ix=0; ix< NumberOfVoxelsAlongX; ix++)
		{
			voxelsSum[ix]=new Svoxel*[NumberOfVoxelsAlongY];
			for (int iy=0; iy< NumberOfVoxelsAlongY; iy++)
			{
				voxelsSum[ix][iy]=new Svoxel[NumberOfVoxelsAlongZ];
				for (int iz=0; iz< NumberOfVoxelsAlongZ; iz++)
				{
					voxelsSum[ix][iy][iz].volumeId=-1;
					voxelsSum[ix][iy][iz].depEnergy=0.;
					voxelsSum[ix][iy][iz].depEnergy2=0.;
					voxelsSum[ix][iy][iz].depEnergyNorm=0.;
					voxelsSum[ix][iy][iz].depEnergyNormError=0.;
					voxelsSum[ix][iy][iz].expDose=0.;
					voxelsSum[ix][iy][iz].halfSize.set(halfXVoxelDimensionX, halfXVoxelDimensionY, halfXVoxelDimensionZ);
					voxelsSum[ix][iy][iz].pos.set(2.*(ix)*halfXVoxelDimensionX  -halfSize.getX()+halfXVoxelDimensionX + centre.getX(), 
			2.*(iy)*halfXVoxelDimensionY  -halfSize.getY()+halfXVoxelDimensionY + centre.getY(), 
			2.*(iz)*halfXVoxelDimensionZ  -halfSize.getZ()+halfXVoxelDimensionZ + centre.getZ());
					voxelsSum[ix][iy][iz].nEvents=0;
				}
			}
		}


// voxels to store and save the single geometry configuration
		voxelsSingle=new Svoxel**[NumberOfVoxelsAlongX];
		for (int ix=0; ix< NumberOfVoxelsAlongX; ix++)
		{
			voxelsSingle[ix]=new Svoxel*[NumberOfVoxelsAlongY];
			for (int iy=0; iy< NumberOfVoxelsAlongY; iy++)
			{
				voxelsSingle[ix][iy]=new Svoxel[NumberOfVoxelsAlongZ];
				for (int iz=0; iz< NumberOfVoxelsAlongZ; iz++)
				{
					voxelsSingle[ix][iy][iz].volumeId=-1;
					voxelsSingle[ix][iy][iz].depEnergy=0.;
					voxelsSingle[ix][iy][iz].depEnergy2=0.;
					voxelsSingle[ix][iy][iz].depEnergyNorm=0.;
					voxelsSingle[ix][iy][iz].depEnergyNormError=0.;
					voxelsSingle[ix][iy][iz].expDose=0.;
					voxelsSingle[ix][iy][iz].halfSize.set(halfXVoxelDimensionX, halfXVoxelDimensionY, halfXVoxelDimensionZ);
					voxelsSingle[ix][iy][iz].pos.set(2.*(ix)*halfXVoxelDimensionX  -halfSize.getX()+halfXVoxelDimensionX + centre.getX(), 
			2.*(iy)*halfXVoxelDimensionY  -halfSize.getY()+halfXVoxelDimensionY + centre.getY(), 
			2.*(iz)*halfXVoxelDimensionZ  -halfSize.getZ()+halfXVoxelDimensionZ + centre.getZ());
					voxelsSingle[ix][iy][iz].nEvents=0;
				}
			}
		}
	}
}

CML2SDWithVoxels::~CML2SDWithVoxels()
{
	if(bSaveROG)
	{
		delete [] voxelsSum;
		delete [] voxelsSingle;
	}
}
void CML2SDWithVoxels::resetVoxelsSingle()
{
	for (int ix=0; ix< NumberOfVoxelsAlongX; ix++)
	{
		for (int iy=0; iy< NumberOfVoxelsAlongY; iy++)
		{
			for (int iz=0; iz< NumberOfVoxelsAlongZ; iz++)
			{
				voxelsSingle[ix][iy][iz].volumeId=-1;
				voxelsSingle[ix][iy][iz].depEnergy=0.;
				voxelsSingle[ix][iy][iz].depEnergy2=0.;
				voxelsSingle[ix][iy][iz].depEnergyNorm=0.;
				voxelsSingle[ix][iy][iz].depEnergyNormError=0.;
				voxelsSingle[ix][iy][iz].expDose=0.;
				voxelsSingle[ix][iy][iz].halfSize.set(halfXVoxelDimensionX, halfXVoxelDimensionY, halfXVoxelDimensionZ);
				voxelsSingle[ix][iy][iz].pos.set(2.*(ix)*halfXVoxelDimensionX  -halfSize.getX()+halfXVoxelDimensionX + centre.getX(), 
		2.*(iy)*halfXVoxelDimensionY  -halfSize.getY()+halfXVoxelDimensionY + centre.getY(), 
		2.*(iz)*halfXVoxelDimensionZ  -halfSize.getZ()+halfXVoxelDimensionZ + centre.getZ());
				voxelsSingle[ix][iy][iz].nEvents=0;
			}
		}
	}
	nSingleTotalEvents=0;
}
G4bool CML2SDWithVoxels::ProcessHits(G4Step *aStep, G4TouchableHistory *ROHist)
{
	
	if (bActive)
	{
		G4double energyDep = aStep->GetTotalEnergyDeposit();
		if (bSaveROG && energyDep>0.) 
		{
			G4int ix, iy, iz;
			G4String volumeName;

			ix=ROHist->GetReplicaNumber(2);
			iy=ROHist->GetReplicaNumber(0);
			iz=ROHist->GetReplicaNumber(1);

			density=aStep->GetPreStepPoint()->GetPhysicalVolume()->GetLogicalVolume()->GetMaterial()->GetDensity();

			voxelMass=voxelVolume*density;
			energyDep/=voxelMass*nRecycling;
			voxelsSum[ix][iy][iz].volumeId=getIdFromVolumeName(aStep->GetPreStepPoint()->GetPhysicalVolume()->GetLogicalVolume()->GetMaterial()->GetName());
			voxelsSum[ix][iy][iz].depEnergy+=energyDep;
			voxelsSum[ix][iy][iz].depEnergy2+=energyDep*energyDep;
			voxelsSum[ix][iy][iz].nEvents++;

			voxelsSingle[ix][iy][iz].volumeId=getIdFromVolumeName(aStep->GetPreStepPoint()->GetPhysicalVolume()->GetLogicalVolume()->GetMaterial()->GetName());
			voxelsSingle[ix][iy][iz].depEnergy+=energyDep;
			voxelsSingle[ix][iy][iz].depEnergy2+=energyDep*energyDep;
			voxelsSingle[ix][iy][iz].nEvents++;
			nTotalEvents++;
			nSingleTotalEvents++;
			if (nTotalEvents%saving_in_ROG_Voxels_every_events==0 && nTotalEvents>0)
			{
				save();
			}
		}
	}
	return true;
}
G4int CML2SDWithVoxels::getIdFromVolumeName(G4String name)
{
	for (int i=0; i<(int)volumeNameIdLink.size(); i++)
	{
		if (volumeNameIdLink[i].volumeName==name)
		{
			return volumeNameIdLink[i].volumeId;
			break;
		}
	}
	return -1;
}
void CML2SDWithVoxels::save()
{
std::cout<< "n. of events collected in the whole ROG phantom for all geometries: "<< nTotalEvents<< G4endl;
std::cout<< "n. of events collected in the whole ROG phantom for the current geometry: "<< nSingleTotalEvents<< G4endl;
	if (nTotalEvents>0)
	{saveData(fullOutFileData, voxelsSum);}
	if (nSingleTotalEvents>0)
	{saveData(fullOutFileDataSingle, voxelsSingle);}
}
void CML2SDWithVoxels::saveData(G4String Filename, Svoxel ***voxels)
{
	std::ofstream out;
	out.open(Filename, std::ios::out);
	out << "Sensitive Detector-Voxels. Total number of events, [mm]->centreX centreY centreZ HalfSizeX HalfSizeY HalfSizeZ minX maxX, minY maxY, minZ maxZ, Dx, Dy, Dz, nX, nY, nZ: \n";
	out <<nTotalEvents<<'\t';
	out <<centre.getX()/mm << '\t' << centre.getY()/mm<< '\t'<< centre.getZ()/mm<<'\t';
	out <<halfSize.getX()/mm << '\t' << halfSize.getY()/mm <<'\t'<< halfSize.getZ()/mm<<'\t';
	out <<(centre.getX()-halfSize.getX())/mm<<'\t'<<(centre.getX()+halfSize.getX())/mm<<'\t';
	out <<(centre.getY()-halfSize.getY())/mm<<'\t'<<(centre.getY()+halfSize.getY())/mm<<'\t';
	out <<(centre.getZ()-halfSize.getZ())/mm<<'\t'<<(centre.getZ()+halfSize.getZ())/mm<<'\t';
	out <<halfXVoxelDimensionX/mm<<'\t'<<halfXVoxelDimensionY/mm<<'\t'<<halfXVoxelDimensionZ/mm<<'\t';
	out <<NumberOfVoxelsAlongX <<'\t'<<NumberOfVoxelsAlongY <<'\t'<<NumberOfVoxelsAlongZ <<'\n';
	out << "Number of physical volumes: "<< volumeNameIdLink.size() << '\n';
	for (int i=0; i<(int)volumeNameIdLink.size(); i++)
	{
		out << volumeNameIdLink[i].volumeName  <<'\t'<< volumeNameIdLink[i].volumeId << G4endl;
	}


	out << "Phys Volume x [mm], y [mm], z [mm], ix, iy, iz, Dose [Gy], Dose2 [Gy^2], nEvents" << G4endl;

	for (int ix=0; ix< NumberOfVoxelsAlongX; ix++)
	{
		for (int iy=0; iy< NumberOfVoxelsAlongY; iy++)
		{
			for (int iz=0; iz< NumberOfVoxelsAlongZ; iz++)
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
	unsigned int ind = fullOutFileData.find(".txt");
	G4String onlyName=fullOutFileData.substr( 0, ind);
	if (val=="")
	{
		static unsigned int indGeom=0;
		char cT[5];
		sprintf(cT,"%d",indGeom);
		fullOutFileDataSingle=onlyName+"Single_"+G4String(cT)+".txt";
		indGeom++;
	}
	else
	{
		fullOutFileDataSingle=onlyName+val+".txt";
	}
}

