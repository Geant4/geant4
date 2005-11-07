//
// File name:     RadmonApplicationPhysicsSetup.hh
// Creation date: Oct 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonApplicationPhysicsSetup.hh,v 1.1 2005-11-07 17:52:36 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   Radmon application generators algoritms setup
//

#ifndef   RADMONAPPLICATIONPHYSICSSETUP_HH
 #define  RADMONAPPLICATIONPHYSICSSETUP_HH

 // Inglude files
 #include "globals.hh"

 // Forward declarations
 class RadmonApplicationOptions;
 class RadmonSubPhysicsListWithLabelFactory;

 class RadmonApplicationPhysicsSetup
 {
  public:
   inline                                       RadmonApplicationPhysicsSetup(const RadmonApplicationOptions & options);
   inline                                      ~RadmonApplicationPhysicsSetup();

  protected:
   G4bool                                       CreateSubPhysicsList(RadmonSubPhysicsListWithLabelFactory * factory);

  // Hidden constructors and operators
                                                RadmonApplicationPhysicsSetup();
                                                RadmonApplicationPhysicsSetup(const RadmonApplicationPhysicsSetup & copy);
   RadmonApplicationPhysicsSetup &              operator=(const RadmonApplicationPhysicsSetup & copy);
   
  // Private attributes
   const RadmonApplicationOptions &             currentOptions;
 };
 
 #include "RadmonApplicationPhysicsSetup.icc"
#endif /* RADMONAPPLICATIONPHYSICSSETUP_HH */
