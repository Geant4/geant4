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

CML2SDWithParticle::CML2SDWithParticle()
: G4VSensitiveDetector("killer_plane"),particles(0)
{
	this->idType=idSD_KillerPlane;
	this->bStopAtPhaseSpace=true;
	this->nTotalParticles=0;
	this->bActive=true;
}
CML2SDWithParticle::CML2SDWithParticle(G4int idType, G4int max_N_particles_in_PhSp_File, G4int seed, G4int nMaxParticlesInRamPhaseSpace, G4String name, G4String PhaseSpaceOutFile, G4bool bSavePhaseSpace, G4bool bStopAtPhaseSpace, SPrimaryParticle *primaryParticleData, G4double  accTargetZPosition)
: G4VSensitiveDetector(name),particles(0)
{
	this->accTargetZPosition=accTargetZPosition;
	this->max_N_particles_in_PhSp_File=max_N_particles_in_PhSp_File;
	this->nMaxParticlesInRamPhaseSpace = nMaxParticlesInRamPhaseSpace;
	this->idType=idType;
	this->primaryParticleData=primaryParticleData;
	this->bActive=true;
	this->nParticle=0;
	this->nTotalParticles=0;
	this->bSavePhaseSpace=bSavePhaseSpace;
	this->bStopAtPhaseSpace=bStopAtPhaseSpace;
	if (this->bSavePhaseSpace)
	{this->particles=new Sparticle[nMaxParticlesInRamPhaseSpace];}
	G4String seedName;
	char a[10];
	sprintf(a,"%d", seed);
	seedName=(G4String)a;

	this->fullOutFileData=PhaseSpaceOutFile+"_"+seedName+".txt";
	this->fullOutFileData=PhaseSpaceOutFile+"_"+seedName+".txt";
}
CML2SDWithParticle::~CML2SDWithParticle()
{
	if (this->bSavePhaseSpace)
	{
		delete [] this->particles;
		delete this->particles;
	}
}
void CML2SDWithParticle::saveHeaderParticles()
{
	std::ofstream out;
	out.open(this->fullOutFileData, std::ios::out);
	out << "Sensitive Detector-Particles"<<G4endl;
	out << "n Total Events,\t x [mm],\t y [mm],\t z [mm],\t dirX,\t dirY,\t dirZ,\t KinEnergy [MeV],\t part Type,\t primary part type,\t nPrimaryPart" << G4endl;
	out.close();
}
void CML2SDWithParticle::saveDataParticles(G4int nParticle)
{
	std::ofstream out;
	out.open(this->fullOutFileData, std::ios::app);
	static G4int nTotalParticles=0;
	for (int i=0; i< nParticle; i++)
	{
		out << nTotalParticles++ << '\t';
//		out << this->particles[i].volumeName << '\t';
		out << this->particles[i].pos.getX()/mm  << '\t';
		out << this->particles[i].pos.getY()/mm  << '\t';
		out << (this->accTargetZPosition + this->particles[i].pos.getZ())/mm  << '\t'; // it translates the current z value in global coordinates to the accelerator local coordinates (only z)
		out << this->particles[i].dir.getX() << '\t';
		out << this->particles[i].dir.getY() << '\t';
		out << this->particles[i].dir.getZ() << '\t';
		out << this->particles[i].kinEnergy/MeV << '\t';
		out << this->particles[i].partPDGE<< '\t';
		out << this->particles[i].primaryParticlePDGE<< '\t';
		out << this->particles[i].nPrimaryPart<< G4endl;
	}
	out.close();
}

G4bool CML2SDWithParticle::ProcessHits(G4Step *aStep, G4TouchableHistory *)
{
	if (this->bActive &&  (CML2AcceleratorConstruction::GetInstance()->getPhysicalVolume()->GetRotation()->isIdentity()))
	{
		G4double energyKin= aStep->GetTrack()->GetKineticEnergy();
		static bool bFirstTime=true;
		if (this->idType==idSD_KillerPlane)
		{
			this->nTotalParticles++;
			aStep->GetTrack()->SetTrackStatus(fStopAndKill);
		}
		else
		{
			if (energyKin>0.)
			{
				this->particles[this->nParticle].volumeName="";
				this->particles[this->nParticle].pos=aStep->GetPreStepPoint()->GetPosition();
				this->particles[this->nParticle].dir=aStep->GetPreStepPoint()->GetMomentumDirection();
				this->particles[this->nParticle].kinEnergy=energyKin;
				this->particles[this->nParticle].partPDGE=aStep->GetTrack()->GetDefinition()->GetPDGEncoding();
				this->particles[this->nParticle].primaryParticlePDGE=this->primaryParticleData->partPDGE;
				this->particles[this->nParticle].nPrimaryPart=this->primaryParticleData->nPrimaryParticle;
				this->nParticle++;
				this->nTotalParticles++;
				if (this->nTotalParticles==this->max_N_particles_in_PhSp_File)
				{
					if (bFirstTime)
					{
						bFirstTime=false;
						this->saveHeaderParticles();
					}
					this->saveDataParticles(this->nParticle);
					this->nParticle=0;
					this->bActive =false;// to stop the phase space creation
				}
				if (this->nParticle==this->nMaxParticlesInRamPhaseSpace)
				{
					if (bFirstTime)
					{
						bFirstTime=false;
						this->saveHeaderParticles();
					}
					this->saveDataParticles(this->nParticle);
					this->nParticle=0;
				}

				Sparticle *particle=new Sparticle;
				particle->dir=aStep->GetPreStepPoint()->GetMomentumDirection();
				particle->pos=aStep->GetPreStepPoint()->GetPosition();
				particle->kinEnergy=energyKin; 
				particle->nPrimaryPart=this->nTotalParticles; // to pass the id of this phase space particle 
				particle->partPDGE=aStep->GetTrack()->GetDefinition()->GetPDGEncoding();
				particle->primaryParticlePDGE=this->primaryParticleData->partPDGE;
				particle->volumeId=-1;
				particle->volumeName="-1";
			}
			if (this->bStopAtPhaseSpace)
			{aStep->GetTrack()->SetTrackStatus(fStopAndKill);}
		}
	}
	else
	{
		if (this->bStopAtPhaseSpace)
		{aStep->GetTrack()->SetTrackStatus(fStopAndKill);}
	}

	return true;
}

void CML2SDWithParticle::save()
{
	if ((this->bActive) && (this->nParticle>0) && (CML2AcceleratorConstruction::GetInstance()->getPhysicalVolume()->GetRotation()->isIdentity()))
	{this->saveDataParticles(this->nParticle);this->nParticle=0;}
}


