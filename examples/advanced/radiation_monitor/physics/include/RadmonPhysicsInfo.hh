//
// File name:     RadmonPhysicsInfo.hh
// Creation date: Nov 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonPhysicsInfo.hh,v 1.1 2005-11-10 08:14:10 capra Exp $
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
