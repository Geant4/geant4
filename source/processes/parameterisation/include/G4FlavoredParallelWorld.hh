// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4FlavoredParallelWorld.hh,v 1.2 1999-04-14 14:25:28 mora Exp $
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

#include "globals.hh"
#include "G4ParticleDefinition.hh"
#include "G4VPhysicalVolume.hh"

class G4FlavoredParallelWorld 
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
    delete [] World;
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
