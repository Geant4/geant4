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
void CML2PrimaryGenerationAction::inizialize(SPrimaryParticle *pData)
{
	rm=new G4RotationMatrix();
	PrimaryGenerationActionMessenger=new CML2PrimaryGenerationActionMessenger(this);
	particle=new Sparticle;
	nParticle=nPhSpParticles=nRandomParticles=0;

	G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();

	gamma=particleTable->FindParticle("gamma");
	electron=particleTable->FindParticle("e-");
	positron=particleTable->FindParticle("e+");
	particleGun=new G4ParticleGun();

	primaryParticleData=pData;
	primaryParticleData->nPrimaryParticle=0;
	primaryParticleData->partPDGE=0;
}

void CML2PrimaryGenerationAction::design(G4double aTZ)
{
	accTargetZPosition=aTZ;
	switch (idParticleSource)
	{
	case id_randomTarget:
		setGunRandom();
		break;
	case id_phaseSpace:
		setGunCalculatedPhaseSpace();
		break;
	}
}

void CML2PrimaryGenerationAction::setGunRandom()
{
	particleGun->SetParticleDefinition(electron);
	particleGun->SetNumberOfParticles(1);
	idCurrentParticleSource=idParticleSource;
}

void CML2PrimaryGenerationAction::setGunCalculatedPhaseSpace()
{
	particles=new Sparticle[nMaxParticlesInRamPhaseSpace];
	particleGun->SetNumberOfParticles(1);
	idCurrentParticleSource=idParticleSource;
}

CML2PrimaryGenerationAction::~CML2PrimaryGenerationAction(void)
{
	delete particleGun;
	delete [] particles;
	delete particles;
}
void CML2PrimaryGenerationAction::GeneratePrimaries(G4Event *anEvent)
{
	static int currentRecycle=nRecycling;
	static G4ThreeVector pos0, dir0;
	if (currentRecycle==nRecycling)
	{
		currentRecycle=0;
		switch (idCurrentParticleSource)
		{
		case id_randomTarget:
				GenerateFromRandom();
			break;
		case id_phaseSpace:
				GenerateFromCalculatedPhaseSpace();
			break;
		}
		pos0=pos;
		dir0=dir;
	}
	currentRecycle++;
	pos=pos0;
	dir=dir0;
	applySourceRotation(); // to follow the accelerator rotation

	primaryParticleData->partPDGE=particleGun->GetParticleDefinition()->GetPDGEncoding();
	primaryParticleData->nPrimaryParticle++;

	particleGun->SetParticleEnergy(ek*MeV);
	particleGun->SetParticlePosition(pos*mm);
	particleGun->SetParticleMomentumDirection((G4ParticleMomentum)dir);
	particleGun->GeneratePrimaryVertex(anEvent);
}
void CML2PrimaryGenerationAction::GenerateFromRandom()
{
	sinTheta=RandGauss::shoot(0., 0.003);
	cosTheta=std::sqrt(1 - sinTheta*sinTheta);
	phi=twopi*G4UniformRand();
	dir.set(sinTheta*std::cos(phi), sinTheta*std::sin(phi), cosTheta);

	ro=G4UniformRand()*GunRadious; 
	alfa=G4UniformRand()*twopi;

	pos.setX(ro*std::sin(alfa));
	pos.setY(ro*std::cos(alfa));
	pos.setZ(-(accTargetZPosition +5.)*mm); // the primary electrons are generated 5 mm before the target
	ek=RandGauss::shoot(GunMeanEnegy, GunStdEnegy);
	nRandomParticles++;
}
void CML2PrimaryGenerationAction::GenerateFromCalculatedPhaseSpace()
{
	static bool bFirstTime=true;
	if (bFirstTime)
	{bFirstTime=false;fillParticlesContainer();}
	if (nParticle==nMaxParticlesInRamPhaseSpace) // once all the particles stored in RAM hae been processed a new set is loaded
	{
		fillParticlesContainer();
		nParticle=0;
	}

	pos=particles[nParticle].pos;
	dir=particles[nParticle].dir;
	ek=particles[nParticle].kinEnergy;
	switch (particles[nParticle].partPDGE)
	{
	case -11:
		particleGun->SetParticleDefinition(positron);
		break;
	case 11:
		particleGun->SetParticleDefinition(electron);
		break;
	case 22:
		particleGun->SetParticleDefinition(gamma);
		break;
	}
	nPhSpParticles++;
	nParticle++;
}
void CML2PrimaryGenerationAction::applySourceRotation()
{
	pos=*rm*pos;
	dir=*rm*dir;
}
void CML2PrimaryGenerationAction::fillParticlesContainer()
{
	static int currentFilePosition=0;
	static int currentFileSize=0;
	int startDataFilePosition;
	std::ifstream in;
	in.open(calculatedPhaseSpaceFileIN, std::ios::in);
	if (in)
	{
		std::cout <<"ERROR phase space file: "<< calculatedPhaseSpaceFileIN << " NOT found. Run abort "<< G4endl;
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

	static bool checkFileRewind=false;
	static bool bRewindTheFile=false;
	static int nPhSpFileRewind=0;

	for (i=0;i<nMaxParticlesInRamPhaseSpace;i++)
	{
		if (bRewindTheFile) // to read the phase space file again to fill the container
		{
			in.close();
		 	in.open(calculatedPhaseSpaceFileIN, std::ios::in);
			in.seekg(startDataFilePosition, std::ios::beg);
			checkFileRewind=true;
			bRewindTheFile=false;
			std::cout<<"\n################\nI have reached the end of the phase space file "<<++nPhSpFileRewind <<" times, I rewind the file\n" << G4endl;
			std::cout <<"loaded " <<i <<"/"<< nMaxParticlesInRamPhaseSpace<<" particles" << G4endl;
		}
		in >> d; 
		in >> x; in >>y; in >> z; 
/*			std::cout <<"x:" <<x << G4endl;
			std::cout <<"y:" <<y << G4endl;
			std::cout <<"z:" <<z << G4endl;*/
		particles[i].pos.set(x,y,z-accTargetZPosition);
		in >> x; in >>y; in >> z; 
		particles[i].dir.set(x,y,z);
		in >> x; 
		particles[i].kinEnergy=x;
		in >> d; 
		particles[i].partPDGE=d;
		in >> d; in >> d; 
		if (in.eof())	{bRewindTheFile=true;}
		if (checkFileRewind)	{checkFileRewind=false;}
	}
	std::cout <<"loaded " <<i <<"/"<< nMaxParticlesInRamPhaseSpace<<" particles" << G4endl;
	currentFilePosition=in.tellg(); // to remind the actual position in the phase space file
	if (currentFilePosition>=currentFileSize) // to read the phase space file again
	{currentFilePosition=startDataFilePosition;} 
	in.close();
}

