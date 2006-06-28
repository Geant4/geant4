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
// File name:     RadmonPhysicsInfo.hh
// Creation date: Nov 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonPhysicsInfo.hh,v 1.2 2006-06-28 13:54:45 gunter Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   Provides informations about a process
//

#ifndef   RADMONPHYSICSINFO_HH
 #define  RADMONPHYSICSINFO_HH
 
 // Include files
 #include "G4String.hh"
 #include "G4ParticleDefinition.hh"
 #include <iostream>
 
 class RadmonPhysicsInfo
 {
  public:
   inline                                       RadmonPhysicsInfo();
                                                RadmonPhysicsInfo(const RadmonPhysicsInfo & copy);
   inline                                      ~RadmonPhysicsInfo();
   
   RadmonPhysicsInfo &                          operator=(const RadmonPhysicsInfo & copy);
   
   G4bool                                       CollidesWith(const RadmonPhysicsInfo & other) const;
   
   inline const G4String &                      GetProcessName(void) const;
   inline const G4ParticleDefinition *          GetParticleDefinition(void) const;
   inline G4double                              GetMinEnergy(void) const;
   inline G4double                              GetMaxEnergy(void) const;

   inline void                                  SetProcessName(const G4String & processName);
   inline void                                  SetParticleDefinition(const G4ParticleDefinition * particleDefinition);
   inline void                                  SetMinEnergy(G4double energy);
   inline void                                  SetMaxEnergy(G4double energy);
   
  private:
  // Private attributes
   G4double                                     minEnergy;
   G4double                                     maxEnergy;
   G4String                                     name;
   const G4ParticleDefinition *                 particle;
 };
 
 std::ostream &                                 operator<<(std::ostream & out, const RadmonPhysicsInfo & info);
 
 // Inline implementations
 #include "RadmonPhysicsInfo.icc"
#endif /* RADMONPHYSICSINFO_HH */
