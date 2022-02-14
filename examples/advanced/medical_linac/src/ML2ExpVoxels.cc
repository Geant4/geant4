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

#include <fstream>

#include "ML2ExpVoxels.hh"
#include "G4SystemOfUnits.hh"

CML2ExpVoxels::CML2ExpVoxels(G4bool bData, G4int saveEvents, G4int seed,
          G4String FileExperimentalData, G4String FileExperimentalDataOut):startCurve(0),
          stopCurve(0),chi2Factor(0)
{
	char a[10];
	sprintf(a,"%d", seed);
	seedName = (G4String)a;
	saving_in_Selected_Voxels_every_events = saveEvents;
	nRecycling = 1;

	
	fullFileOut = FileExperimentalDataOut+seedName+".m";
	fullFileIn = FileExperimentalData;
	nParticle = nTotalEvents = 0;

// define the extremes of global-volume containing all experimental voxels
	G4double extr = 100000000000.;
	minZone.set(extr, extr, extr);
	maxZone.set(-extr, -extr, -extr);
	bHasExperimentalData = bData;
}

CML2ExpVoxels::~CML2ExpVoxels(void)
{
	delete [] startCurve;
	delete [] stopCurve;
	delete [] chi2Factor;
	delete [] nVoxelsgeometry;
}
G4bool CML2ExpVoxels::loadData(void)
{
	bHasExperimentalData = true;
	std::ifstream in;

	Svoxel voxel;
	voxel.volumeId = 0;
	G4ThreeVector pos, halfSize;
	G4double expDose;

	in.open(fullFileIn, std::ios::in);
	if (in)
	{
		G4String appo;
		char a[1000];
		in.getline(a,1000,'\n');
		headerText1 = (G4String)a;
		in.getline(a,1000,'\n');
		in >> nCurves;
		startCurve = new G4int[nCurves];
		stopCurve = new G4int[nCurves];
		chi2Factor = new G4double[nCurves];
		for (int i = 0; i < nCurves; i++)
		{
			chi2Factor[i] = 0.;
			in >> startCurve[i]; 
			in >> stopCurve[i]; 
			in >> chi2Factor[i];
		}
		in.getline(a,1000,'\n');
		in.getline(a,1000,'\n');
		headerText2 = (G4String)a;
		std::string line;

		while ( !in.eof() )
		{
			in >> pos;
			in >> halfSize;

			if (bHasExperimentalData)
			{
				in >> expDose;
				voxel.expDose = expDose/100.*(joule/kg); // input data in cGy
			}
			else
			{
				voxel.expDose = 0.;
			}
			voxel.pos=pos;
			voxel.halfSize = halfSize;
			voxel.depEnergy = 0.;
			voxel.depEnergy2 = 0.;
			voxel.nEvents = 0;
			voxel.depEnergyNorm = 0.;
			voxel.depEnergyNormError = 0.;
			vec_voxels.push_back(voxel);

			// calculate the actual extremes of the global-volume containing all the experimental data
			if ( minZone.getX()>pos.getX()-halfSize.getX() )
			{ minZone.setX(pos.getX()-halfSize.getX()); }
			if ( maxZone.getX()<pos.getX()+halfSize.getX() )
			{ maxZone.setX(pos.getX()+halfSize.getX()); }

			if ( minZone.getY()>pos.getY()-halfSize.getY() )
			{ minZone.setY(pos.getY()-halfSize.getY()); }
			if ( maxZone.getY()<pos.getY()+halfSize.getY() )
			{ maxZone.setY(pos.getY()+halfSize.getY()); }

			if ( minZone.getZ()>pos.getZ()-halfSize.getZ() )
			{ minZone.setZ(pos.getZ()-halfSize.getZ()); }
			if ( maxZone.getZ()<pos.getZ()+halfSize.getZ() )
			{ maxZone.setZ(pos.getZ()+halfSize.getZ()); }
		}
	}
	else
	{
		G4cout << "ERROR I can't find the experimental data file" << G4endl;
		return false;
	}
	in.close();

        nVoxelsgeometry = new G4int[(G4int) vec_voxels.size()];
        resetNEventsInVoxels();

	return true;
}
void CML2ExpVoxels::resetNEventsInVoxels()
{
    for (int i=0; i<(int) vec_voxels.size(); i++ )
    {nVoxelsgeometry[i] = 0;}
}

void CML2ExpVoxels::add(const G4Step* aStep)
{
	G4ThreeVector pos;
	G4double depEnergy, density, voxelVolume;
	
	pos = aStep->GetPreStepPoint()->GetPosition();
	depEnergy = aStep->GetTotalEnergyDeposit();
	density = aStep->GetPreStepPoint()->GetPhysicalVolume()->GetLogicalVolume()->GetMaterial()->GetDensity();

	G4ThreeVector minPos, maxPos;
	G4bool newEvent=false;
	G4double voxelMass, dose;

	// check if the event is inside the global-volume
	if (minZone.getX() <= pos.getX() && pos.getX() < maxZone.getX() &&
		minZone.getY() <= pos.getY() && pos.getY() < maxZone.getY() &&
		minZone.getZ() <= pos.getZ() && pos.getZ() < maxZone.getZ())
	{
		// look for the voxel containing the event
		for (int i = 0; i < (int)vec_voxels.size(); i++)
		{
			minPos = vec_voxels[i].pos-vec_voxels[i].halfSize;
			maxPos = vec_voxels[i].pos+vec_voxels[i].halfSize;
			if ( minPos.getX() <= pos.getX() && pos.getX() < maxPos.getX() &&
				 minPos.getY() <= pos.getY() && pos.getY() < maxPos.getY() &&
				 minPos.getZ() <= pos.getZ() && pos.getZ() < maxPos.getZ() )
			{
				voxelVolume = vec_voxels[i].halfSize.getX()*vec_voxels[i].halfSize.getY()*vec_voxels[i].halfSize.getZ()*8.;
				voxelMass = density*voxelVolume;
				// calculate the dose
				dose=depEnergy/(voxelMass*nRecycling);
				vec_voxels[i].nEvents++;
				nVoxelsgeometry[i]++;
				vec_voxels[i].depEnergy += dose;
				vec_voxels[i].depEnergy2 += dose*dose;
				newEvent = true;

				Sparticle *particle = new Sparticle;
				particle -> dir = aStep -> GetPreStepPoint() -> GetMomentumDirection();
				particle -> pos = aStep -> GetPreStepPoint() -> GetPosition();
				particle -> kinEnergy = dose; // I use the same kinEnergy name to store the dose
				particle -> nPrimaryPart = -1;
				particle -> partPDGE = aStep -> GetTrack() -> GetDefinition() -> GetPDGEncoding();
				particle -> primaryParticlePDGE = -1;
				particle -> volumeId = i; // voxel index where the dose is accumulating
				particle -> volumeName = "-1";
			}
		}
		if ( newEvent )
		{
			// save data
			nTotalEvents++;
			if ( nTotalEvents%saving_in_Selected_Voxels_every_events == 0 && nTotalEvents > 0 )
			{
				saveResults();
			}
		}
	}
}

G4int CML2ExpVoxels::getMinNumberOfEvents()
{
	int n = vec_voxels[0].nEvents;
	for (int i = 0; i < (int)vec_voxels.size(); i++)
	{ 
		if ( n > vec_voxels[i].nEvents )
		{ n = vec_voxels[i].nEvents; }
	}
	return n;
}
G4int CML2ExpVoxels::getMaxNumberOfEvents()
{
	int n = nVoxelsgeometry[0];
	for ( int i = 0; i < (int)vec_voxels.size(); i++)
	{ 
		if ( n < nVoxelsgeometry[i] )
		{ n = nVoxelsgeometry[i]; }
	}
	return n;
}
void CML2ExpVoxels::saveHeader()
{
	std::ofstream out;
	out.open(fullFileOut, std::ios::out);
	out << "% " << headerText1 << G4endl;
	out << "n"  << seedName	<< "=" << nCurves << ";" << G4endl;
	out << "fh" << seedName	<< "=[" << G4endl;
	for (int i = 0; i< nCurves; i++)
	{
		out << startCurve[i] << '\t';
		out << stopCurve[i] << '\t';
		out << chi2Factor[i] << G4endl;
	}
	out << "];" << G4endl;
	out << "% x [mm], y [mm], z [mm], Dx [mm], Dy [mm], Dz [mm], expDose [Gy], Calculated dose [Gy], Calculated dose2 [Gy^2], nEvents, normDose [Gy], normDoseError [Gy]";
	out << G4endl;
	out.close();	
}

void CML2ExpVoxels::saveResults()
{
	if (nTotalEvents > 0)
	{
		calculateNormalizedEd(vec_voxels);
		saveHeader();
		std::ofstream out;
		out.open(fullFileOut, std::ios::app);
		out << "d" << seedName << "=[" << G4endl;
		for (int i=0; i<(int)vec_voxels.size(); i++)
		{
			out << vec_voxels[i].pos.getX()/mm << '\t' << vec_voxels[i].pos.getY()/mm << '\t' << vec_voxels[i].pos.getZ()/mm << '\t';
			out << vec_voxels[i].halfSize.getX()/mm << '\t' << vec_voxels[i].halfSize.getY()/mm << '\t' << vec_voxels[i].halfSize.getZ()/mm << '\t';
			out << vec_voxels[i].expDose/(joule/kg) << '\t' << vec_voxels[i].depEnergy/(joule/kg) << '\t' << vec_voxels[i].depEnergy2/((joule/kg)*(joule/kg)) << '\t' << vec_voxels[i].nEvents << '\t';
			out << vec_voxels[i].depEnergyNorm/(joule/kg) << '\t' << vec_voxels[i].depEnergyNormError/(joule/kg);
			out << G4endl;
		}
		out << "];" << G4endl;
		out.close();
	}
}

void CML2ExpVoxels::calculateNormalizedEd(std::vector <Svoxel> &vox)
{
	int i,j;
	G4double cs, cc;
	int n;
	G4double d2, dd;
        G4double v;
	for (j = 0; j < nCurves; j++)
	{
		cs = cc = 0.;
		for (i = startCurve[j]-1;i<stopCurve[j];i++)
		{
			cs += vox[i].depEnergy*vox[i].expDose;
			cc += vox[i].depEnergy*vox[i].depEnergy;
		}
		if (cc>0.)
		{
			chi2Factor[j] = cs/cc;
		}
		for (i = startCurve[j]-1; i < stopCurve[j]; i++)
		{
			dd = vox[i].depEnergy*vox[i].depEnergy;
			d2 = vox[i].depEnergy2;
			n = vox[i].nEvents;
			vox[i].depEnergyNorm = chi2Factor[j]*vox[i].depEnergy;
			v = n*d2-dd;
			if (v < 0.) { v=0; }
			if (n > 1) 	{ vox[i].depEnergyNormError = chi2Factor[j]*std::sqrt(v/(n-1)); }
			if (n == 1) { vox[i].depEnergyNormError = vox[i].depEnergyNorm; }
		}
	}
}
