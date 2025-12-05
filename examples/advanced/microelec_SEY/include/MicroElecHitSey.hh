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
/// \file MicroElecHitSey.cc
/// \brief Implementation of the MicroElecHitSey class
//
//-

#ifndef __MICROELECHITSEY_H__
#define __MICROELECHITSEY_H__ 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"

class MicroElecHitSey : public G4VHit
{
	
	public:

		MicroElecHitSey();
		~MicroElecHitSey();

		MicroElecHitSey(const MicroElecHitSey&);

		const MicroElecHitSey& operator=(const MicroElecHitSey&);
		int operator==(const MicroElecHitSey&) const;

		inline void* operator new(size_t);
		inline void  operator delete(void*);

		void Print();

	public:

		void SetNbPrim			(G4int nbprim)		{ NbPrim=nbprim;				};
		void SetNbSec			(G4int nbsec)		{ NbSec = nbsec;				};
		void SetNbSup50			(G4int nbsup50)		{ NbSup50 = nbsup50;			};
		void SetPDGEncoding		(G4int pdgencoding)	{ PDGEncoding = pdgencoding;	};
		void SetParentID		(G4int parent)      { ParentID = parent;			};
		void SetTrackID			(G4int track)       { trackID = track;				};
		void SetParticleType	(G4String nameT)	{ ParticleType = nameT;			};
		void SetParticleName	(G4String name)		{ ParticleName = name;			};
		void SetVolumeName		(G4String nameV)	{ VolumeName = nameV;			};
		void SetZ   			(G4int z)        	{ Z = z; 						};
		void SetA   			(G4int a)        	{ A = a; 						};
		void SetVertexKineticEnergy	(G4double VQ)	{ VertexKineticEnergy = VQ; 	};
		void SetVertexPos(G4ThreeVector xyz)		{ VertexPos = xyz; 				};
		void SetVertexMomentum(G4ThreeVector xyz)	{ VertexMomentum = xyz; 		};
		void SetPreStepKineticEnergy (G4double Q)	{ PreStepKineticEnergy = Q; 	};
		void SetPostStepKineticEnergy (G4double Q)	{ PostStepKineticEnergy = Q; 	};
		void SetEdep      		(G4double de)       { Edep = de; 					};
		void SetNi_Edep      	(G4double NIde)     { Ni_Edep = NIde; 				};
		void SetStepLength		(G4double length)   { StepLength = length;			};
		void SetPrePos       		(G4ThreeVector xyz) { PrePos = xyz; 			};
		void SetPostPos       		(G4ThreeVector xyz) { PostPos = xyz; 			};
		void SetPreStepMomentum		(G4ThreeVector xyz) { PreStepMomentum = xyz; 	};
		void SetPostStepMomentum	(G4ThreeVector xyz) { PostStepMomentum = xyz;	};
			
		
		G4int			GetNbPrim()				{ return NbPrim;				};
		G4int			GetNbSec()				{ return NbSec;					};
		G4int			GetNbSup50()			{ return NbSup50;				};
		G4int			GetPDGEncoding ()		{ return PDGEncoding; 			};
		G4int			GetParentID ()			{ return ParentID; 				};
		G4int			GetTrackID ()			{ return trackID; 				};
		G4String		GetParticleType ()		{ return ParticleType;			};
		G4String		GetParticleName ()		{ return ParticleName;			};
		G4String		GetVolumeName ()		{ return VolumeName;			};
		G4int			GetZ ()					{ return Z;						};
		G4int			GetA ()					{ return A; 					};
		G4double		GetVertexKineticEnergy (){ return VertexKineticEnergy;	};
		G4ThreeVector GetVertexPos()			{ return VertexPos; 			};
		G4ThreeVector GetVertexMomentum()		{ return VertexMomentum ; 		};
		G4double		GetPreStepKineticEnergy(){ return PreStepKineticEnergy; };
		G4double		GetPostStepKineticEnergy(){ return PostStepKineticEnergy;};
		G4double		GetEdep ()				{ return Edep; 					};
		G4double		GetNi_Edep ()			{ return Ni_Edep; 				};
		G4ThreeVector	GetPrePos ()  			{ return PrePos; 				};
		G4ThreeVector	GetPostPos ()  			{ return PostPos; 				};
		G4ThreeVector	GetPreStepMomentum ()  	{ return PreStepMomentum; 		};
		G4ThreeVector	GetPostStepMomentum ()  { return PostStepMomentum; 		};
		G4double		GetStepLength ()		{ return StepLength; 			};
			
	//-------------------
	//New due to the use of touchables to identify the different sensitive detectors
		void SetVolumeCopyNumber(int value){VolumeCopyNumber = value;};
		G4int GetVolumeCopyNumber(){return VolumeCopyNumber;};
			
	private:

		G4int			NbPrim;
		G4int			NbSec;
		G4int			NbSup50;
		G4int			PDGEncoding;
		G4int			ParentID;
		G4int			trackID;
		G4String		ParticleType;
		G4String		ParticleName;
		G4String		VolumeName;
		G4int			Z;
		G4int			A;
		G4double		VertexKineticEnergy;
		G4double		PreStepKineticEnergy;
		G4double		PostStepKineticEnergy;
		G4double		Edep;
		G4double		Ni_Edep;
		G4double		StepLength;
		G4ThreeVector	PrePos;
		G4ThreeVector	PostPos;
		G4ThreeVector	VertexPos;
		G4ThreeVector	PreStepMomentum;
		G4ThreeVector	PostStepMomentum;
		G4ThreeVector	VertexMomentum;

	//-------------------
	//New due to the use of touchables to identify the different sensitive detectors
		G4int VolumeCopyNumber;
    
};

typedef G4THitsCollection<MicroElecHitSey> MicroElecHitSeyCollection;
extern G4ThreadLocal G4Allocator<MicroElecHitSey>* MicroElecHitSeyAllocator;
inline void* MicroElecHitSey::operator new(size_t){
    void *aHit;
	if (!MicroElecHitSeyAllocator) { MicroElecHitSeyAllocator = new G4Allocator<MicroElecHitSey>; }
	aHit = (void *) MicroElecHitSeyAllocator->MallocSingle();
    return aHit;
}
inline void MicroElecHitSey::operator delete(void *aHit){MicroElecHitSeyAllocator->FreeSingle((MicroElecHitSey*) aHit);}
#endif