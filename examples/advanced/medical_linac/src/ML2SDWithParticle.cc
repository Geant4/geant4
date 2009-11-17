#include "ML2SDWithParticle.h"
#include "ML2ExpVoxels.h"

CML2SDWithParticle::CML2SDWithParticle()
: G4VSensitiveDetector("killer_plane"),particles(0)
{
	this->idType=idSD_KillerPlane;
	this->bStopAtPhaseSpace=true;
	this->nTotalParticles=0;
	this->bActive=true;
}
CML2SDWithParticle::CML2SDWithParticle(G4int idType, G4int max_N_particles_in_PhSp_File, G4int seed, G4int nMaxParticlesInRamPhaseSpace, G4String name, G4String PhaseSpaceOutFile, G4bool bSavePhaseSpace, G4bool bStopAtPhaseSpace, SPrimaryParticle *primaryParticleData)
: G4VSensitiveDetector(name),particles(0)
{
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
	this->bContinueRun=true;
	G4String seedName;
	char a[10];
	sprintf(a,"%d", seed);
	seedName=(G4String)a;

	char *MyDirOut=new char[1000];
	MyDirOut=getenv("G4MYFILEOUT");
	G4String myDirOut=(G4String) MyDirOut;
	this->fullOutFileData=myDirOut+"/"+PhaseSpaceOutFile+"_"+seedName+".txt";
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
		out << this->particles[i].volumeName << '\t';
		out << this->particles[i].pos.getX()/mm  << '\t';
		out << this->particles[i].pos.getY()/mm  << '\t';
		out << this->particles[i].pos.getZ()/mm  << '\t';
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

G4bool CML2SDWithParticle::ProcessHits(G4Step *aStep, G4TouchableHistory *ROHist)
{
	if (this->bActive)
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
				this->primaryParticleData->nParticlesInPhSp++;
				if (this->nTotalParticles==this->max_N_particles_in_PhSp_File)
				{
					this->bContinueRun=false;
					this->bSavePhaseSpace=false;
					if (bFirstTime)
					{
						bFirstTime=false;
						this->saveHeaderParticles();
					}
					this->saveDataParticles(this->nParticle);
					this->nParticle=0;
				}
				if (this->nParticle==this->nMaxParticlesInRamPhaseSpace)
				{
					if (this->bSavePhaseSpace)
					{
						if (bFirstTime)
						{
							bFirstTime=false;
							this->saveHeaderParticles();
						}
						this->saveDataParticles(this->nParticle);
					}
					if (this->nTotalParticles>=max_N_particles_in_PhSp_File)
					{
						this->bSavePhaseSpace=false;
					}
					this->nParticle=0;
				}
			}
		if (this->bStopAtPhaseSpace)
		{aStep->GetTrack()->SetTrackStatus(fStopAndKill);}
		}
	}
	return true;
}

