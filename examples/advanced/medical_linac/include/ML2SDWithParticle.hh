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


#ifndef CML2SensitiveDetectorParticleH
#define CML2SensitiveDetectorParticleH

#include "G4VSensitiveDetector.hh"
#include "G4TouchableHistory.hh"
#include "G4VTouchable.hh"
#include "ML2SinputData.hh"


class CML2ReadOutGeometryVoxels;

class CML2SDWithParticle : public G4VSensitiveDetector
{
public:
	CML2SDWithParticle();
	CML2SDWithParticle(G4int idType, G4int max_N_particles_in_PhSp_File, G4int seed, G4int nMaxParticlesInRamPhaseSpace, G4String name, G4String PhaseSpaceOutFile, G4bool bSavePhaseSpace, G4bool bStopAtVolatilePhaseSpace, SPrimaryParticle *primaryParticleData, G4double  accTargetZPosition);
	~CML2SDWithParticle(void);
	G4bool ProcessHits(G4Step *aStep, G4TouchableHistory *ROHist);
	G4int getTotalNumberOfParticles(){return nTotalParticles;}
	inline CML2SDWithParticle* getCML2SensitiveDetectorParticle(){return this;}
	inline Sparticle getParticle(int i){return particles[i];}
	inline void setActive(G4bool act){bActive=act;}
	void save();
private:
	void saveDataParticles(G4int nParticle);
	void saveHeaderParticles();

	G4ThreeVector halfSize;
	G4ThreeVector pos;
	SPrimaryParticle *primaryParticleData;
	Sparticle *particles;
	G4String fullOutFileData;
	G4int nTotalParticles, nParticle;
	G4int idType, nMaxParticlesInRamPhaseSpace, max_N_particles_in_PhSp_File;
	G4bool bStopAtPhaseSpace, bSavePhaseSpace, bActive; 

	G4double accTargetZPosition;
};

#endif
