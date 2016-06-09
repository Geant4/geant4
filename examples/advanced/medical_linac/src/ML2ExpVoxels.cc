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


#include "ML2ExpVoxels.hh"

#include <fstream>

CML2ExpVoxels::CML2ExpVoxels(G4bool bHasExperimentalData, G4int saving_in_Selected_Voxels_every_events, G4int seed, G4String FileExperimentalData):startCurve(0), stopCurve(0),chi2Factor(0)
{
	char a[10];
	sprintf(a,"%d", seed);
	this->seedName=(G4String)a;
	this->saving_in_Selected_Voxels_every_events=saving_in_Selected_Voxels_every_events;

	this->fullFileOut="EXP_seed_"+this->seedName+".m";
#ifdef ML2FILEOUT
// full path of the output data
	char *MyDir=new char[1000];
	MyDir=getenv("ML2FILEOUT");
	G4String myDir=(G4String) MyDir;
	this->fullFileOut=myDir+"/EXP_seed_"+this->seedName+".m";
#endif

	this->fullFileIn=FileExperimentalData;
#ifdef ML2FILEIN
// full path of the experimental data
	MyDir=getenv("ML2FILEIN");
	myDir=(G4String) MyDir;
	this->fullFileIn=myDir+"/"+FileExperimentalData;
#endif

	this->nParticle=this->nTotalEvents=0;

// define the extremes of global-volume containing all experimental voxels
	G4double extr=100000000000.;
	this->minZone.set(extr, extr, extr);
	this->maxZone.set(-extr, -extr, -extr);
	this->bHasExperimentalData=bHasExperimentalData;
}

CML2ExpVoxels::~CML2ExpVoxels(void)
{
	delete [] this->startCurve;
	delete [] this->stopCurve;
	delete [] this->chi2Factor;
	delete this->startCurve;
	delete this->stopCurve;
	delete this->chi2Factor;
}
G4bool CML2ExpVoxels::loadData(void)
{
	this->bHasExperimentalData=true;
	std::ifstream in;

	Svoxel voxel;
	G4ThreeVector pos, halfSize;
	G4double expDose;

	in.open(this->fullFileIn, std::ios::in);
	if (in !=0)
	{
		G4String appo;
		char a[1000];
		in.getline(a,1000,'\n'); this->headerText1=(G4String)a;
		in.getline(a,1000,'\n');
		in >> this->nCurves;
		this->startCurve=new G4int[this->nCurves];
		this->stopCurve=new G4int[this->nCurves]; 
		this->chi2Factor=new G4double[this->nCurves];
		for (int i=0; i< this->nCurves; i++)
		{
			this->chi2Factor[i]=0.;
			in >> this->startCurve[i]; 
			in >> this->stopCurve[i]; 
			in >> this->chi2Factor[i];
		}
		in.getline(a,1000,'\n');
		in.getline(a,1000,'\n'); this->headerText2=(G4String)a;

		while (!in.eof())
		{
			in >> pos;
			in >> halfSize; 
			if (this->bHasExperimentalData)
			{
				in >> expDose;
				voxel.expDose=expDose/100.*(joule/kg); // input data in cGy
			}
			else
			{
				voxel.expDose=0.;
			}
			voxel.pos=pos;
			voxel.halfSize=halfSize;
			voxel.depEnergy=0.;
			voxel.depEnergy2=0.;
			voxel.nEvents=0;
			voxel.depEnergyNorm=0.;
			voxel.depEnergyNormError=0.;
			this->voxels.push_back(voxel);

// calculate the actual extremes of the global-volume containing all the experimental data
			if (this->minZone.getX()>pos.getX()-halfSize.getX())
			{this->minZone.setX(pos.getX()-halfSize.getX());}
			if (this->maxZone.getX()<pos.getX()+halfSize.getX())
			{this->maxZone.setX(pos.getX()+halfSize.getX());}

			if (this->minZone.getY()>pos.getY()-halfSize.getY())
			{this->minZone.setY(pos.getY()-halfSize.getY());}
			if (this->maxZone.getY()<pos.getY()+halfSize.getY())
			{this->maxZone.setY(pos.getY()+halfSize.getY());}

			if (this->minZone.getZ()>pos.getZ()-halfSize.getZ())
			{this->minZone.setZ(pos.getZ()-halfSize.getZ());}
			if (this->maxZone.getZ()<pos.getZ()+halfSize.getZ())
			{this->maxZone.setZ(pos.getZ()+halfSize.getZ());}

		}
	}
	else
	{
		std::cout << "ERROR I can't find the experimental data file" << G4endl;
		return false;
	}
	in.close();
	return true;
}
void CML2ExpVoxels::add(G4ThreeVector pos, G4double depEnergy, G4double density)
{
	G4ThreeVector minPos, maxPos;
	G4bool newEvent=false;
	G4double voxelMass, dose;

// check if the event is inside the global-volume 
	if (this->minZone.getX()<= pos.getX() && pos.getX()<this->maxZone.getX() && 
		this->minZone.getY()<= pos.getY() && pos.getY()<this->maxZone.getY() && 
		this->minZone.getZ()<= pos.getZ() && pos.getZ()<this->maxZone.getZ())
	{
// look for the voxel containing the event
		for (int i=0; i<(int)this->voxels.size(); i++)
		{
			minPos=this->voxels[i].pos-this->voxels[i].halfSize;
			maxPos=this->voxels[i].pos+this->voxels[i].halfSize;
			if (minPos.getX()<= pos.getX() && pos.getX()<maxPos.getX() && 
				minPos.getY()<= pos.getY() && pos.getY()<maxPos.getY() && 
				minPos.getZ()<= pos.getZ() && pos.getZ()<maxPos.getZ())
			{
				G4double a1, a2, a3;
				a1=this->voxels[i].halfSize.getX();
				a2=this->voxels[i].halfSize.getY();
				a3=this->voxels[i].halfSize.getZ();

				voxelMass=density*this->voxels[i].halfSize.getX()*this->voxels[i].halfSize.getY()*this->voxels[i].halfSize.getZ();
// calculate the dose 
				dose=depEnergy/voxelMass;
				this->voxels[i].nEvents++;
				this->voxels[i].depEnergy+=dose;
				this->voxels[i].depEnergy2+=dose*dose;
				newEvent=true;
			}
		}
		if (newEvent)
		{
// save data
			this->nTotalEvents++;
			if (this->nTotalEvents%this->saving_in_Selected_Voxels_every_events==0 && this->nTotalEvents>0)
			{
				this->saveResults(this->fullFileOut, this->voxels);
			}
		}
	}
}

G4int CML2ExpVoxels::getMinNumberOfEvents()
{
	int n=10000000;
	for (int i=this->startCurve[0];i<this->stopCurve[this->nCurves];i++)
	{ 
		if (n>voxels[i].nEvents){n = voxels[i].nEvents;}
	}
	return n;
}
void CML2ExpVoxels::saveHeader(G4String fullOutFileName)
{
	std::ofstream out;
	out.open(fullOutFileName, std::ios::out);
	out <<"% "<< this->headerText1 << G4endl;
	out <<"n"<< this->seedName<<"="<< this->nCurves<<";" << G4endl;
	out <<"fh"<< this->seedName<<"=["<< G4endl;
	for (int i=0; i< this->nCurves; i++)
	{
		out << this->startCurve[i] << '\t';
		out << this->stopCurve[i] << '\t';
		out << this->chi2Factor[i]<< G4endl;
	}
	out << "];"<<G4endl;
	out <<"% x [mm], y [mm], z [mm], expDose [Gy], Calculated dose [Gy], Calculated dose2 [Gy^2], nEvents, normDose [Gy], normDoseError [Gy]"<< G4endl;
	out.close();	
}

void CML2ExpVoxels::saveResults(G4String fullOutFileName, std::vector <Svoxel> voxels)
{
	this->calculateNormalizedEd(voxels);
	this->saveHeader(fullOutFileName);
	std::ofstream out;
	out.open(fullOutFileName, std::ios::app);
	out <<"d"<< this->seedName<<"=["<< G4endl;
	for (int i=0; i<(int)voxels.size(); i++)
	{
		out <<voxels[i].pos.getX()/mm<<'\t'<<voxels[i].pos.getY()/mm<<'\t'<<voxels[i].pos.getZ()/mm<<'\t';
		out <<voxels[i].expDose/(joule/kg)<<'\t'<<voxels[i].depEnergy/(joule/kg)<<'\t'<<voxels[i].depEnergy2/((joule/kg)*(joule/kg))<<'\t'<<voxels[i].nEvents<< '\t';
		out <<voxels[i].depEnergyNorm/(joule/kg)<<'\t'<<voxels[i].depEnergyNormError/(joule/kg)<<G4endl;
	}
	out << "];"<<G4endl;
	out.close();
}
void CML2ExpVoxels::calculateNormalizedEd(std::vector <Svoxel> &voxels)
{
	int i,j;
	G4double cs, cc;
	int n;
	G4double d2, dd;
	G4double v, appoo;
	for (j=0;j<this->nCurves;j++)
	{
		cs=cc=0.;
		for (i=this->startCurve[j];i<this->stopCurve[j];i++)
		{
			cs+=voxels[i].depEnergy*voxels[i].expDose;
			cc+=voxels[i].depEnergy*voxels[i].depEnergy;
		}
		if (cc>0.)
		{
			this->chi2Factor[j]=cs/cc; 
		}
		for (i=this->startCurve[j];i<this->stopCurve[j];i++)
		{
			dd=voxels[i].depEnergy*voxels[i].depEnergy;
			d2=voxels[i].depEnergy2;
			n=voxels[i].nEvents;
			voxels[i].depEnergyNorm=chi2Factor[j]*voxels[i].depEnergy;
			appoo=chi2Factor[j];
			appoo=voxels[i].depEnergy;
			appoo=voxels[i].depEnergyNorm;
			v=n*d2-dd;
			if (v<0.){v=0;}
			if (n>1){voxels[i].depEnergyNormError=this->chi2Factor[j]*std::sqrt(v/(n-1));}
			if (n==1){voxels[i].depEnergyNormError=voxels[i].depEnergyNorm;}
		}
	}
}


