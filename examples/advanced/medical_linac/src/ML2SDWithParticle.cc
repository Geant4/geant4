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

#include "ML2SDWithParticle.hh"
#include "ML2ExpVoxels.hh"
#include "ML2AcceleratorConstruction.hh"
#include "G4SystemOfUnits.hh"

CML2SDWithParticle::CML2SDWithParticle()
: G4VSensitiveDetector("killer_plane"),particles(0)
{
	idType=idSD_KillerPlane;
	bStopAtPhaseSpace=true;
	nTotalParticles=0;
	bActive=true;
}

CML2SDWithParticle::CML2SDWithParticle(G4int id, G4int maxPartFile,
                                       G4int seed, G4int nMaxPart, G4String name,
                                       G4String PhaseSpaceOutFile, G4bool bSave, G4bool bStop,
                                       SPrimaryParticle *pData, G4double  ZPos)
: G4VSensitiveDetector(name),particles(0)
{
	accTargetZPosition=ZPos;
	max_N_particles_in_PhSp_File=maxPartFile;
	nMaxParticlesInRamPhaseSpace = nMaxPart;
	idType=id;
	primaryParticleData=pData;
	bActive=true;
	nParticle=0;
	nTotalParticles=0;
	bSavePhaseSpace=bSave;
	bStopAtPhaseSpace=bStop;
	if (bSavePhaseSpace)
	{particles=new Sparticle[nMaxPart];}
	G4String seedName;
	char a[10];
	sprintf(a,"%d", seed);
	seedName=(G4String)a;

	fullOutFileData=PhaseSpaceOutFile+"_"+seedName+".txt";
	fullOutFileData=PhaseSpaceOutFile+"_"+seedName+".txt";
}
CML2SDWithParticle::~CML2SDWithParticle()
{
        if (particles!=0)
	{
		delete [] particles;
	}
}
void CML2SDWithParticle::saveHeaderParticles()
{
	std::ofstream out;
	out.open(fullOutFileData, std::ios::out);
	out << "Sensitive Detector-Particles"<<G4endl;
	out << "n Total Events,\t x [mm],\t y [mm],\t z [mm],\t dirX,\t dirY,\t dirZ,\t KinEnergy [MeV],\t part Type,\t primary part type,\t nPrimaryPart" << G4endl;
	out.close();
}
void CML2SDWithParticle::saveDataParticles(G4int nPart)
{
	std::ofstream out;
	out.open(fullOutFileData, std::ios::app);
	static G4int nTotParticles=0;
	for (int i=0; i< nPart; i++)
	{
		out << nTotParticles++ << '\t';
//		out << particles[i].volumeName << '\t';
		out << particles[i].pos.getX()/mm  << '\t';
		out << particles[i].pos.getY()/mm  << '\t';
		out << (accTargetZPosition + particles[i].pos.getZ())/mm  << '\t'; // it translates the current z value in global coordinates to the accelerator local coordinates (only z)
		out << particles[i].dir.getX() << '\t';
		out << particles[i].dir.getY() << '\t';
		out << particles[i].dir.getZ() << '\t';
		out << particles[i].kinEnergy/MeV << '\t';
		out << particles[i].partPDGE<< '\t';
		out << particles[i].primaryParticlePDGE<< '\t';
		out << particles[i].nPrimaryPart<< G4endl;
	}
	out.close();
}

G4bool CML2SDWithParticle::ProcessHits(G4Step *aStep, G4TouchableHistory *)
{
	if (bActive &&  (CML2AcceleratorConstruction::GetInstance()->getPhysicalVolume()->GetRotation()->isIdentity()))
	{
		G4double energyKin= aStep->GetTrack()->GetKineticEnergy();
		static bool bFirstTime=true;
		if (idType==idSD_KillerPlane)
		{
			nTotalParticles++;
			aStep->GetTrack()->SetTrackStatus(fStopAndKill);
		}
		else
		{
			if (energyKin>0.)
			{
				particles[nParticle].volumeName="";
				particles[nParticle].pos=aStep->GetPreStepPoint()->GetPosition();
				particles[nParticle].dir=aStep->GetPreStepPoint()->GetMomentumDirection();
				particles[nParticle].kinEnergy=energyKin;
				particles[nParticle].partPDGE=aStep->GetTrack()->GetDefinition()->GetPDGEncoding();
				particles[nParticle].primaryParticlePDGE=primaryParticleData->partPDGE;
				particles[nParticle].nPrimaryPart=primaryParticleData->nPrimaryParticle;
				nParticle++;
				nTotalParticles++;
				if (nTotalParticles==max_N_particles_in_PhSp_File)
				{
					if (bFirstTime)
					{
						bFirstTime=false;
						saveHeaderParticles();
					}
					saveDataParticles(nParticle);
					nParticle=0;
					bActive =false;// to stop the phase space creation
				}
				if (nParticle==nMaxParticlesInRamPhaseSpace)
				{
					if (bFirstTime)
					{
						bFirstTime=false;
						saveHeaderParticles();
					}
					saveDataParticles(nParticle);
					nParticle=0;
				}

				Sparticle *particle=new Sparticle;
				particle->dir=aStep->GetPreStepPoint()->GetMomentumDirection();
				particle->pos=aStep->GetPreStepPoint()->GetPosition();
				particle->kinEnergy=energyKin; 
				particle->nPrimaryPart=nTotalParticles; // to pass the id of this phase space particle 
				particle->partPDGE=aStep->GetTrack()->GetDefinition()->GetPDGEncoding();
				particle->primaryParticlePDGE=primaryParticleData->partPDGE;
				particle->volumeId=-1;
				particle->volumeName="-1";
			}
			if (bStopAtPhaseSpace)
			{aStep->GetTrack()->SetTrackStatus(fStopAndKill);}
		}
	}
	else
	{
		if (bStopAtPhaseSpace)
		{aStep->GetTrack()->SetTrackStatus(fStopAndKill);}
	}

	return true;
}

void CML2SDWithParticle::save()
{
	if ((bActive) && (nParticle>0) && (CML2AcceleratorConstruction::GetInstance()->getPhysicalVolume()->GetRotation()->isIdentity()))
	{saveDataParticles(nParticle);nParticle=0;}
}
