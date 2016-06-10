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
//
// $Id: G4VRestProcess.hh 71231 2013-06-12 13:06:28Z gcosmo $
//
// 
// ------------------------------------------------------------
//	GEANT 4 class header file 
//
//	History: first implementation, based on object model of
//	2nd December 1995, G.Cosmo
// 
// Class Description 
//  Abstract class which defines the public behavior of
//  physics interactions at rest.
//
// ------------------------------------------------------------
//   New Physics scheme           18 Dec. 1996  H.Kurahige
// ------------------------------------------------------------
//   modified                     25 Feb. 1997  H.Kurahige
//   modified                      8 Mar. 1997  H.Kurahige
//   modified                     26 Mar. 1997  H.Kurahige
//   modified                     16 Apr. 1997  L.Urban    
//   modified                     17 Dec. 1997  L.Urban    
//   modified for new ParticleChange 12 Mar. 1998  H.Kurashige


#ifndef G4VRestProcess_h
#define G4VRestProcess_h 1

#include "globals.hh"
#include "G4ios.hh"

#include "G4VProcess.hh"


class G4VRestProcess : public G4VProcess 
{
  //  Abstract class which defines the public behavior of
  //  physics interactions at rest.

  public:
      G4VRestProcess(const G4String&  ,
		     G4ProcessType   aType = fNotDefined );
      G4VRestProcess(G4VRestProcess& );

      virtual ~G4VRestProcess();

  public:   //  with description
      virtual G4double AtRestGetPhysicalInteractionLength(
                             const G4Track& track,
			     G4ForceCondition* condition
			    );

      virtual G4VParticleChange* AtRestDoIt(
			     const G4Track& ,
			     const G4Step&
			    );

     //  no operation in  PostStepDoIt and  AlongStepDoIt
      virtual G4double AlongStepGetPhysicalInteractionLength(
                             const G4Track&,
			     G4double  ,
			     G4double  ,
			     G4double& ,
	                     G4GPILSelection*
			   ){ return -1.0; };

      virtual G4double PostStepGetPhysicalInteractionLength(
                             const G4Track& ,
			     G4double   ,
			     G4ForceCondition* 
			    ) { return -1.0; };

     //  no operation in  PostStepDoIt and  AlongStepDoIt
      virtual G4VParticleChange* PostStepDoIt(
			     const G4Track& ,
			     const G4Step&
			    ) {return 0;};

      virtual G4VParticleChange* AlongStepDoIt(
			     const G4Track& ,
			     const G4Step& 
			    ) {return 0;};
 
  protected: //  with description

      virtual G4double GetMeanLifeTime(const G4Track& aTrack,G4ForceCondition* condition)=0;
      //  Calculates the mean life-time (i.e. for decays) of the
      //  particle at rest due to the occurence of the given process,
      //  or converts the probability of interaction (i.e. for
      //  annihilation) into the life-time of the particle for the
      //  occurence of the given process.

 private:
  // hide default constructor and assignment operator as private 
      G4VRestProcess();
      G4VRestProcess & operator=(const G4VRestProcess &right);
};



#endif



