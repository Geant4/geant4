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
// File name:     RadmonPhysicsInfo.hh
// Creation date: Nov 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonPhysicsInfo.hh,v 1.3 2006/06/29 16:17:19 gunter Exp $
// Tag:           $Name: geant4-09-00 $
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
