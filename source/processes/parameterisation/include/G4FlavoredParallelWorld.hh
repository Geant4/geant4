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
//
// $Id: G4FlavoredParallelWorld.hh,v 1.5 2001-07-11 10:08:24 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//  
//---------------------------------------------------------------
//
//  G4FlavoredParallelWorld.hh
//
//  Description:
//    Internal class to keep ParticleType X Parallel world 
//    relationship.
//
//  History:
//    June 98: Verderi && MoraDeFreitas - "G4ParallelWorld" becomes
//             "G4FlavoredParallelWorld".
//    Mars 98: Verderi && MoraDeFreitas - First Implementation.
//
//---------------------------------------------------------------

#ifndef  G4FlavoredParallelWorld_hh
#define  G4FlavoredParallelWorld_hh

#include "G4VFlavoredParallelWorld.hh"

class G4ParticleDefinition;
class G4VPhysicalVolume;

class G4FlavoredParallelWorld : public G4VFlavoredParallelWorld
{
public:
  //
  // Constructor
  G4FlavoredParallelWorld(G4ParticleDefinition *pParticleType,
			  G4VPhysicalVolume *pWorld) {
    ParticleType=pParticleType;
    World=pWorld;
  }
  //
  // Destructor
  ~G4FlavoredParallelWorld() {
    delete World;
    World = 0;
  }
  // Get/Set
  inline G4ParticleDefinition* GetTheParticleType() const {
    return ParticleType;
  }
  
  inline G4VPhysicalVolume* GetThePhysicalVolumeWorld() const {
    return World;
  }
  
  // operator == 
  inline G4bool 
  operator == (const G4FlavoredParallelWorld& pw) const {
    return (this==&pw) ? true : false;
  }
  
  
private:
  G4ParticleDefinition *ParticleType;
  G4VPhysicalVolume *World;
};

#endif 
// end of #ifndef G4FlavoredParallelWorld_hh
