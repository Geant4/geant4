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


#include "ML2PrimaryGenerationAction.hh"
#include "ML2PrimaryGenerationActionMessenger.hh"

using namespace CLHEP;

CML2PrimaryGenerationAction::CML2PrimaryGenerationAction(SPrimaryParticle *primaryParticleData)
:particleGun(0),gamma(0),electron(0),positron(0),primaryParticleData(0),particles(0),firstFileParticle(0),lastLoadedParticle(0)
{
	this->PrimaryGenerationActionMessenger=new CML2PrimaryGenerationActionMessenger(this);
	this->particle=new Sparticle;
	this->firstFileParticle=new Sparticle;
	this->lastLoadedParticle=new Sparticle;
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
void CML2PrimaryGenerationAction::design()
{
	this->calculatedPhaseSpaceFileIN=this->calculatedPhaseSpaceFileIN;
#ifdef ML2FILEIN
	char *MyDirIn=new char[1000];
	MyDirIn=getenv("ML2FILEIN");
	G4String myDirIn=(G4String) MyDirIn;
	this->calculatedPhaseSpaceFileIN=myDirIn+"/"+this->calculatedPhaseSpaceFileIN;
#endif

// std::cout <<"this->calculatedPhaseSpaceFileIN " << this->calculatedPhaseSpaceFileIN<< G4endl;
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
	this->particleGun->SetNumberOfParticles(this->nIdenticalParticles);
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
	switch (this->idCurrentParticleSource)
	{
	case id_randomTarget:
			this->GenerateFromRandom();
		break;
	case id_phaseSpace:
			this->GenerateFromCalculatedPhaseSpace();
		break;
	}


// std::cout <<"this->pos " <<this->pos << '\t' <<"this->dir " <<this->dir << '\t'<<"this->ek " << this->ek<< G4endl;


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
	this->pos.setZ(-1305.*mm);

	this->ek=RandGauss::shoot(this->GunMeanEnegy, this->GunStdEnegy);
	this->nRandomParticles++;
}
void CML2PrimaryGenerationAction::GenerateFromCalculatedPhaseSpace()
{
	static bool bFirstTime=true;
	if (bFirstTime)
	{bFirstTime=false;this->fillParticlesContainer();}
	if (nParticle==this->nMaxParticlesInRamPhaseSpace)
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

void CML2PrimaryGenerationAction::fillParticlesContainer()
{
	static int currentFilePosition=0;
	static int currentFileSize=0;
	int startDataFilePosition;
	std::ifstream in;
	in.open(this->calculatedPhaseSpaceFileIN, std::ios::in);
	static bool bFirstTime=true;
	if (bFirstTime)
	{
		in.seekg(-1,std::ios::end);
		currentFileSize=in.tellg();
		in.seekg(0,std::ios::beg);
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
	for (i=0;i<this->nMaxParticlesInRamPhaseSpace;i++)
	{
		in >> d; 
		in >> x; in >>y; in >> z; 
		this->particles[i].pos.set(x,y,z);
		in >> x; in >>y; in >> z; 
		this->particles[i].dir.set(x,y,z);
		in >> x; 
		this->particles[i].kinEnergy=x;
		in >> d; 
		this->particles[i].partPDGE=d;
		in >> d; in >> d; 
		if (in.tellg() >= currentFileSize-2)
		{
			in.seekg(startDataFilePosition, std::ios::beg);
			checkFileRewind=true;
			in >> d; 
			in >> x; in >>y; in >> z; 
			this->particles[i].pos.set(x,y,z);
			in >> x; in >>y; in >> z; 
			this->particles[i].dir.set(x,y,z);
			in >> x; 
			this->particles[i].kinEnergy=x;
			in >> d; 
			this->particles[i].partPDGE=d;
			in >> d; in >> d; 
		}
		if (checkFileRewind)
		{
			checkFileRewind=false;
		}
		if (i==0 && bFirstTime)
		{
			this->firstFileParticle->pos=this->particles[i].pos;
			this->firstFileParticle->dir=this->particles[i].dir;
			this->firstFileParticle->kinEnergy=this->particles[i].kinEnergy;
			this->firstFileParticle->nPrimaryPart=this->particles[i].nPrimaryPart;
			this->firstFileParticle->partPDGE=this->particles[i].partPDGE;
			this->firstFileParticle->primaryParticlePDGE=this->particles[i].primaryParticlePDGE;
		}
	}
	currentFilePosition=in.tellg();
	if (currentFilePosition>=currentFileSize)
	{currentFilePosition=startDataFilePosition;}
	in.close();
	if (bFirstTime)
	{
		bFirstTime=false;
	}
}
bool CML2PrimaryGenerationAction::itIsTheSameParticle(Sparticle *p1, Sparticle *p2)
{
	if (p1->nPrimaryPart==p2->nPrimaryPart &&
		p1->kinEnergy==p2->kinEnergy &&
		p1->dir==p2->dir &&
		p1->pos==p2->pos &&
		p1->partPDGE==p2->partPDGE &&
		p1->primaryParticlePDGE==p2->primaryParticlePDGE)
	{return true;}
	return false;
}

