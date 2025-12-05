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
// INCL++ intra-nuclear cascade model
// Alain Boudard, CEA-Saclay, France
// Joseph Cugnon, University of Liege, Belgium
// Jean-Christophe David, CEA-Saclay, France
// Pekka Kaitaniemi, CEA-Saclay, France, and Helsinki Institute of Physics, Finland
// Sylvie Leray, CEA-Saclay, France
// Davide Mancusi, CEA-Saclay, France
//
#define INCLXX_IN_GEANT4_MODE 1

#include "globals.hh"

/*
 * G4INCLAntinucleiAtrestEntryChannel.cc
 *
 * 
 * \author Olivier Lourgo
 */
#include "G4INCLAntinucleiAtrestEntryChannel.hh"
#include "G4INCLNbarAtrestEntryChannel.hh"
#include "G4INCLPbarAtrestEntryChannel.hh"
#include "G4INCLPhaseSpaceGenerator.hh"


namespace G4INCL{

	AntinucleiAtrestEntryChannel::AntinucleiAtrestEntryChannel(Nucleus *n, Cluster *ac, ThreeVector pos1, ThreeVector pos2)
	  :theNucleus(n), theantiComposite(ac), Posnbar(pos1), Pospbar(pos2){}
	
	AntinucleiAtrestEntryChannel::AntinucleiAtrestEntryChannel(Nucleus *n, Particle *p):theNucleus(n),Meson(p){}

	AntinucleiAtrestEntryChannel::~AntinucleiAtrestEntryChannel(){}

	ThreeVector AntinucleiAtrestEntryChannel::getAnnihilationPosition(ThreeVector nbarPos, ThreeVector pbarPos){ //Choose between the pbar or nbar annihilation position
	    if((nbarPos - pbarPos).mag2() <= ParticleTable::getLargestNuclearRadius(-theantiComposite->getA(),- theantiComposite->getZ())){
	    //If the annihilions positions are close (the radius of a deuteron) then we have 2 sources for the meson star
	    	return ThreeVector(999.,999.,999.);
	    }
	    if(nbarPos.mag2() <= pbarPos.mag2())
	    	return nbarPos;
	    else
	    	return pbarPos;
	}

	ParticleList AntinucleiAtrestEntryChannel::makeMesonStar(){
		ParticleList Antiparticles = theantiComposite->getParticles();
		Particle *nbar=nullptr;
		Particle *pbar=nullptr;
        for(ParticleIter p =Antiparticles.begin(), e=Antiparticles.end(); p!=e; ++p){
        	if((*p)->getType()==antiProton)
        		pbar = *p;
        	else if((*p)->getType()==antiNeutron)
        		nbar = *p;
        	else
        		INCL_ERROR("ERROR : something else than antiNeutron or antiProton in antiComposite");
        }
        PbarAtrestEntryChannel *pbarChannel = new PbarAtrestEntryChannel(theNucleus, pbar);
        NbarAtrestEntryChannel *nbarChannel = new NbarAtrestEntryChannel(theNucleus, nbar);
        
        ParticleList TotalStarList = pbarChannel->makeMesonStar(); //pbar in first because polarisation of the dbar (Coulomb)
        ParticleList nbarMesonStar = nbarChannel->makeMesonStar();
        pbarListSize = (G4int)TotalStarList.size();

        for(ParticleIter p=nbarMesonStar.begin(), e=nbarMesonStar.end(); p!=e; ++p){
        	TotalStarList.push_back(*p);
        }
        Pospbar = pbarChannel->getAnnihilationPosition();
        Posnbar = nbarChannel->getAnnihilationPosition();
        
        G4double EnergyofFinalMesonStar = 0;
        G4int a=theNucleus->getA();
        G4int z=theNucleus->getZ();
        G4int stra=theNucleus->getS();
        if(theNucleus->getAnnihilationType()==DNbarPPbarPType){
        	EnergyofFinalMesonStar = theantiComposite->getMass() + (ParticleTable::getTableMass(a+2,z+2,stra)- ParticleTable::getTableMass(a,z+1,stra));
        }
        else if(theNucleus->getAnnihilationType()==DNbarPPbarNType || theNucleus->getAnnihilationType()==DNbarNPbarPType){
        	//Correction for all but they all cancel out !  
        	EnergyofFinalMesonStar = theantiComposite->getMass() + (ParticleTable::getTableMass(a+2,z+1,stra)- ParticleTable::getTableMass(a,z,stra));

        }
        else if(theNucleus->getAnnihilationType()==DNbarNPbarNType){
        	EnergyofFinalMesonStar = theantiComposite->getMass() + (ParticleTable::getTableMass(a+2,z,stra) - ParticleTable::getTableMass(a,z-1,stra));
        }
        PhaseSpaceGenerator::generate(EnergyofFinalMesonStar, TotalStarList);
        
        return TotalStarList; 
	}

	IAvatarList AntinucleiAtrestEntryChannel::bringMesonStar(ParticleList const &pL, Nucleus * const n){
		ThreeVector ann_position = getAnnihilationPosition(Posnbar,Pospbar);
		IAvatarList theAvatarList;
		G4int cnt=1;
		if (ann_position.getX() == 999. && ann_position.getY() == 999. && ann_position.getZ() == 999.){
			INCL_DEBUG("Particle are close to each other : 2 sources of annihilation "<< '\n');
			for(ParticleIter p=pL.begin(), e=pL.end(); p!=e; ++p){
			    if(cnt <= pbarListSize){
			    	(*p)->setPosition(Pospbar);
			    }
			    else
			    	(*p)->setPosition(Posnbar);
			    theAvatarList.push_back(new ParticleEntryAvatar(0.0, n, *p, ADAR));
			    cnt++;
		    }
		} 
		else{ 
			for(ParticleIter p=pL.begin(), e=pL.end(); p!=e; ++p){
			    (*p)->setPosition(ann_position);
			    theAvatarList.push_back(new ParticleEntryAvatar(0.0, n, *p, ADAR));
			}
		}
		return theAvatarList;
	}

    void AntinucleiAtrestEntryChannel::fillFinalState(FinalState *fs){
    	const G4double energyBefore = Meson->getEnergy();
    	fs->addEnteringParticle(Meson);
    	INCL_DEBUG("Entering antiComposite annihilation product added " << '\n');
    	fs->setTotalEnergyBeforeInteraction(energyBefore);
    }
}
   
