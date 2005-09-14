//
// File name:     RadmonApplication.hh
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonApplication.hh,v 1.2 2005-09-14 12:30:57 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   Radmon application
//

#ifndef   RADMONAPPLICATION_HH
 #define  RADMONAPPLICATION_HH

 // Include files
 #include "globals.hh"
 
 // Forward declarations
 class RadmonApplicationOptions;
 class G4RunManager;
 class RadmonDetectorLayout;
 class RadmonVDetectorEntitiesConstructorsFactory;
 class G4VisManager;
 class G4UImanager;
 class RadmonDetectorMessenger;
 class G4UIsession;
 class G4UIdirectory;

 class RadmonApplication
 {
  public:
                                                RadmonApplication(const RadmonApplicationOptions & options);
                                               ~RadmonApplication();

   inline G4bool                                Valid(void) const;

  private:
   G4bool                                       CreateEntityConstructors(void);
   G4bool                                       RunMacro(const RadmonApplicationOptions & options, const char * fileName);
  
  // Hidden constructors and operators
                                                RadmonApplication(const RadmonApplication & copy);
   RadmonApplication &                          operator=(const RadmonApplication & copy);
   
  // Private attributes
   G4bool                                       valid;

   G4RunManager *                               runManager;
   RadmonDetectorLayout *                       layout;
   RadmonVDetectorEntitiesConstructorsFactory * factory;
   #ifdef    G4VIS_USE
    G4VisManager *                              visManager;
   #endif /* G4VIS_USE */
   G4UImanager *                                uiManager;
   RadmonDetectorMessenger *                    messenger;
   G4UIsession *                                session;
   G4UIdirectory *                              directory;
 };
 
 #include "RadmonApplication.icc"
#endif /* RADMONAPPLICATION_HH */
