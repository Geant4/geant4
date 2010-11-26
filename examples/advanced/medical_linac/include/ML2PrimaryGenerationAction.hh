//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
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

#ifndef CML2PrimaryGenerationActionH
#define CML2PrimaryGenerationActionH


#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleGun.hh"
#include "G4Event.hh"
#include "G4Timer.hh"
#include "Randomize.hh" 
#include "G4RunManager.hh"

#include "G4ParticleDefinition.hh"
#include "ML2SinputData.hh"
#include "ML2SDWithParticle.hh"
#include "ML2SDWithVoxels.hh"
#include "ML2PrimaryGenerationActionMessenger.hh"

class G4ParticleGun;
class G4ParticleDefinition;
class CML2PrimaryGenerationActionMessenger;

class CML2PrimaryGenerationAction : public G4VUserPrimaryGeneratorAction
{
public:
	CML2PrimaryGenerationAction(void);
	static CML2PrimaryGenerationAction* GetInstance(void);

	~CML2PrimaryGenerationAction(void);
	void design(G4double accTargetZPosition);
	void GeneratePrimaries(G4Event *anEvent);
	void inizialize(SPrimaryParticle *primaryParticleData);
	inline void setNRecycling(G4int val){this->nRecycling=val;};
	inline void setNLoopsPhSpParticles(G4int val){this->nLoopsPhSpParticles=val;};
	inline void setNMaxParticlesInRamPhaseSpace(G4int val){this->nMaxParticlesInRamPhaseSpace=val;std::cout<<"Current nMaxParticlesInRamPhaseSpace: " << this->nMaxParticlesInRamPhaseSpace<< G4endl;};

	inline void setGunMeanEnergy(G4double val){this->GunMeanEnegy=val;};
	inline void setGunStdEnergy(G4double val){this->GunStdEnegy=val;};
	inline void setGunRadious(G4double val){this->GunRadious=val;};
	inline void setCalculatedPhaseSpaceFileIN(G4String val){this->calculatedPhaseSpaceFileIN=val;};
	inline void setSourceTypeName(G4String val)
	{
		this->sourceTypeName=val;
		if (this->sourceTypeName=="randomTarget")
		{
			this->idParticleSource=id_randomTarget;
		}
		else if (this->sourceTypeName=="phaseSpace")
		{
			this->idParticleSource=id_phaseSpace;
		}
	};
	inline void setRotation(G4RotationMatrix *val){this->rm=val;};
	inline G4double getNrecycling(){return this->nRecycling;};
	inline G4int getSourceTypeName(){return this->idParticleSource;};

private:
	void setGunRandom();
	void setGunCalculatedPhaseSpace();
	void GenerateFromRandom();
	void GenerateFromCalculatedPhaseSpace();
	void fillParticlesContainer();
	void applySourceRotation();

	static CML2PrimaryGenerationAction * instance;

	G4int nBeam, nRecycling, nLoopsPhSpParticles, idGunType, nMaxParticlesInRamPhaseSpace, idParticleSource;
	G4double GunMeanEnegy, GunStdEnegy, GunRadious;
	G4String calculatedPhaseSpaceFileIN;

	CML2PrimaryGenerationActionMessenger *PrimaryGenerationActionMessenger;
	G4double accTargetZPosition;

	G4ThreeVector dir, pos;
	G4double ek;
	G4RotationMatrix *rm;

	G4Timer myTime;
	G4double sinTheta, cosTheta, phi;
	G4double ro, alfa;
	G4ParticleGun *particleGun;
	G4ParticleDefinition *gamma;
	G4ParticleDefinition *electron;
	G4ParticleDefinition *positron;
	G4int nEventsPerRun;
	SPrimaryParticle *primaryParticleData;
	Sparticle *particles, *particle;
	int nParticle, nPhSpParticles, nRandomParticles, idCurrentParticleSource;
	G4String sourceTypeName;
};

#endif
