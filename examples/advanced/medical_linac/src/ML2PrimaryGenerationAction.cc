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

#include "ML2PrimaryGenerationAction.hh"

using namespace CLHEP;

CML2PrimaryGenerationAction::CML2PrimaryGenerationAction(void)
:particleGun(0),gamma(0),electron(0),positron(0),primaryParticleData(0),particles(0)
{
}
CML2PrimaryGenerationAction* CML2PrimaryGenerationAction::instance = 0;

CML2PrimaryGenerationAction* CML2PrimaryGenerationAction::GetInstance(void)
{
  if (instance == 0)
    {
      instance = new CML2PrimaryGenerationAction();
    }
  return instance;
}
void CML2PrimaryGenerationAction::inizialize(SPrimaryParticle *primaryParticleData)
{
	this->rm=new G4RotationMatrix();
	this->PrimaryGenerationActionMessenger=new CML2PrimaryGenerationActionMessenger(this);
	this->particle=new Sparticle;
	this->nParticle=this->nPhSpParticles=this->nRandomParticles=0;

	G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();

	this->gamma=particleTable->FindParticle("gamma");
	this->electron=particleTable->FindParticle("e-");
	this->positron=particleTable->FindParticle("e+");
	this->particleGun=new G4ParticleGun();

	this->primaryParticleData=primaryParticleData;
	this->primaryParticleData->nPrimaryParticle=0;
	this->primaryParticleData->partPDGE=0;
}
void CML2PrimaryGenerationAction::design(G4double accTargetZPosition)
{
	this->accTargetZPosition=accTargetZPosition;
	switch (this->idParticleSource)
	{
	case id_randomTarget:
		this->setGunRandom();
		break;
	case id_phaseSpace:
		this->setGunCalculatedPhaseSpace();
		break;
	}
}
void CML2PrimaryGenerationAction::setGunRandom()
{
	this->particleGun->SetParticleDefinition(this->electron);
	this->particleGun->SetNumberOfParticles(1);
	this->idCurrentParticleSource=this->idParticleSource;
}
void CML2PrimaryGenerationAction::setGunCalculatedPhaseSpace()
{
	this->particles=new Sparticle[this->nMaxParticlesInRamPhaseSpace];
	this->particleGun->SetNumberOfParticles(1);
	this->idCurrentParticleSource=this->idParticleSource;
}

CML2PrimaryGenerationAction::~CML2PrimaryGenerationAction(void)
{
	delete this->particleGun;
	delete [] this->particles;
	delete this->particles;
}
void CML2PrimaryGenerationAction::GeneratePrimaries(G4Event *anEvent)
{
	static int currentRecycle=this->nRecycling;
	static G4ThreeVector pos0, dir0;
	if (currentRecycle==this->nRecycling)
	{
		currentRecycle=0;
		switch (this->idCurrentParticleSource)
		{
		case id_randomTarget:
				this->GenerateFromRandom();
			break;
		case id_phaseSpace:
				this->GenerateFromCalculatedPhaseSpace();
			break;
		}
		pos0=this->pos;
		dir0=this->dir;
	}
	currentRecycle++;
	this->pos=pos0;
	this->dir=dir0;
	this->applySourceRotation(); // to follow the accelerator rotation

	this->primaryParticleData->partPDGE=this->particleGun->GetParticleDefinition()->GetPDGEncoding();
	this->primaryParticleData->nPrimaryParticle++;

	this->particleGun->SetParticleEnergy(this->ek*MeV);
	this->particleGun->SetParticlePosition(this->pos*mm);
	this->particleGun->SetParticleMomentumDirection((G4ParticleMomentum)this->dir);
	this->particleGun->GeneratePrimaryVertex(anEvent);
}
void CML2PrimaryGenerationAction::GenerateFromRandom()
{
	this->sinTheta=RandGauss::shoot(0., 0.003);
	this->cosTheta=std::sqrt(1 - this->sinTheta*this->sinTheta);
	this->phi=twopi*G4UniformRand();
	this->dir.set(this->sinTheta*std::cos(this->phi), this->sinTheta*std::sin(this->phi), this->cosTheta);

	this->ro=G4UniformRand()*this->GunRadious; 
	this->alfa=G4UniformRand()*twopi;

	this->pos.setX(this->ro*std::sin(this->alfa));
	this->pos.setY(this->ro*std::cos(this->alfa));
	this->pos.setZ(-(this->accTargetZPosition +5.)*mm); // the primary electrons are generated 5 mm before the target
	this->ek=RandGauss::shoot(this->GunMeanEnegy, this->GunStdEnegy);
	this->nRandomParticles++;
}
void CML2PrimaryGenerationAction::GenerateFromCalculatedPhaseSpace()
{
	static bool bFirstTime=true;
	if (bFirstTime)
	{bFirstTime=false;this->fillParticlesContainer();}
	if (this->nParticle==this->nMaxParticlesInRamPhaseSpace) // once all the particles stored in RAM hae been processed a new set is loaded
	{
		this->fillParticlesContainer();
		this->nParticle=0;
	}

	this->pos=this->particles[this->nParticle].pos;
	this->dir=this->particles[this->nParticle].dir;
	this->ek=this->particles[this->nParticle].kinEnergy;
	switch (this->particles[this->nParticle].partPDGE)
	{
	case -11:
		this->particleGun->SetParticleDefinition(this->positron);
		break;
	case 11:
		this->particleGun->SetParticleDefinition(this->electron);
		break;
	case 22:
		this->particleGun->SetParticleDefinition(this->gamma);
		break;
	}
	this->nPhSpParticles++;
	this->nParticle++;
}
void CML2PrimaryGenerationAction::applySourceRotation()
{
	this->pos=*this->rm*this->pos;
	this->dir=*this->rm*this->dir;
}
void CML2PrimaryGenerationAction::fillParticlesContainer()
{
	static int currentFilePosition=0;
	static int currentFileSize=0;
	int startDataFilePosition;
	std::ifstream in;
	in.open(this->calculatedPhaseSpaceFileIN, std::ios::in);
	if (in==0)
	{
		std::cout <<"ERROR phase space file: "<< this->calculatedPhaseSpaceFileIN << " NOT found. Run abort "<< G4endl;
		G4RunManager::GetRunManager()->AbortRun(true);
	}

	static bool bFirstTime=true;
	if (bFirstTime)
	{
		in.seekg(-1,std::ios::end);
		currentFileSize=in.tellg();
		in.seekg(0,std::ios::beg);
		bFirstTime=false;
	}

	char a[1000];
	in.getline(a,1000,'\n');
	in.getline(a,1000,'\n');
	startDataFilePosition=in.tellg();
	if (currentFilePosition>0)
	{in.seekg(currentFilePosition, std::ios::beg);}
	int i;
	G4double x,y,z;
	G4int d;
	G4String s;

	static bool checkFileRewind=false;
	static bool bRewindTheFile=false;
	static int nPhSpFileRewind=0;

	for (i=0;i<this->nMaxParticlesInRamPhaseSpace;i++)
	{
		if (bRewindTheFile) // to read the phase space file again to fill the container
		{
			in.close();
		 	in.open(this->calculatedPhaseSpaceFileIN, std::ios::in);
			in.seekg(startDataFilePosition, std::ios::beg);
			checkFileRewind=true;
			bRewindTheFile=false;
			std::cout<<"\n################\nI have reached the end of the phase space file "<<++nPhSpFileRewind <<" times, I rewind the file\n" << G4endl;
			std::cout <<"loaded " <<i <<"/"<< this->nMaxParticlesInRamPhaseSpace<<" particles" << G4endl;
		}
		in >> d; 
		in >> x; in >>y; in >> z; 
/*			std::cout <<"x:" <<x << G4endl;
			std::cout <<"y:" <<y << G4endl;
			std::cout <<"z:" <<z << G4endl;*/
		this->particles[i].pos.set(x,y,z-this->accTargetZPosition);
		in >> x; in >>y; in >> z; 
		this->particles[i].dir.set(x,y,z);
		in >> x; 
		this->particles[i].kinEnergy=x;
		in >> d; 
		this->particles[i].partPDGE=d;
		in >> d; in >> d; 
		if (in.eof())	{bRewindTheFile=true;}
		if (checkFileRewind)	{checkFileRewind=false;}
	}
	std::cout <<"loaded " <<i <<"/"<< this->nMaxParticlesInRamPhaseSpace<<" particles" << G4endl;
	currentFilePosition=in.tellg(); // to remind the actual position in the phase space file
	if (currentFilePosition>=currentFileSize) // to read the phase space file again
	{currentFilePosition=startDataFilePosition;} 
	in.close();
}

