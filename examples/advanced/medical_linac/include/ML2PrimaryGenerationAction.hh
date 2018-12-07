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
	inline void setNRecycling(G4int val){nRecycling=val;}
	inline void setNLoopsPhSpParticles(G4int val){nLoopsPhSpParticles=val;}
	inline void setNMaxParticlesInRamPhaseSpace(G4int val){nMaxParticlesInRamPhaseSpace=val;G4cout<<"Current nMaxParticlesInRamPhaseSpace: " << nMaxParticlesInRamPhaseSpace<< G4endl;}

	inline void setGunMeanEnergy(G4double val){GunMeanEnergy=val;}
	inline void setGunStdEnergy(G4double val){GunStdEnergy=val;}
	inline void setGunRadius(G4double val){GunRadius=val;}
	inline void setCalculatedPhaseSpaceFileIN(G4String val){calculatedPhaseSpaceFileIN=val;}
	inline void setSourceTypeName(G4String val)
	{
		sourceTypeName=val;
		if (sourceTypeName=="randomTarget")
		{
			idParticleSource=id_randomTarget;
		}
		else if (sourceTypeName=="phaseSpace")
		{
			idParticleSource=id_phaseSpace;
		}
	}
	inline void setRotation(G4RotationMatrix *val){rm=val;};
	inline G4int getNrecycling(){return nRecycling;};
	inline G4int getSourceTypeName(){return idParticleSource;};

private:
	void setGunRandom();
	void setGunCalculatedPhaseSpace();
	void GenerateFromRandom();
	void GenerateFromCalculatedPhaseSpace();
	void fillParticlesContainer();
	void applySourceRotation();

	static CML2PrimaryGenerationAction * instance;

	G4int nRecycling, nLoopsPhSpParticles, nMaxParticlesInRamPhaseSpace, idParticleSource;
	G4double GunMeanEnergy, GunStdEnergy, GunRadius;
	G4String calculatedPhaseSpaceFileIN;

	CML2PrimaryGenerationActionMessenger *PrimaryGenerationActionMessenger;
	G4double accTargetZPosition;

	G4ThreeVector dir, pos;
	G4double ek;
	G4RotationMatrix *rm;

	G4Timer myTime;
	G4double sinTheta, cosTheta, phi;
	G4double rho, alpha;
	G4ParticleGun *particleGun;
	G4ParticleDefinition *gamma;
	G4ParticleDefinition *electron;
	G4ParticleDefinition *positron;
	SPrimaryParticle *primaryParticleData;
	Sparticle *particles, *particle;
	int nParticle, nPhSpParticles, nRandomParticles, idCurrentParticleSource;
	G4String sourceTypeName;
};

#endif
