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
// $Id: G4FlavoredParallelWorld.hh,v 1.6 2006/06/29 21:09:20 gunter Exp $
// GEANT4 tag $Name: geant4-08-01 $
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
