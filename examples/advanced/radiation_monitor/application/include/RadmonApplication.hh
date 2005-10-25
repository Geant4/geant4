//
// File name:     RadmonApplication.hh
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonApplication.hh,v 1.4 2005-10-25 16:39:12 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   Radmon application
//

#ifndef   RADMONAPPLICATION_HH
 #define  RADMONAPPLICATION_HH

 // Include files
 #include "globals.hh"
 #include "RadmonApplicationDetectorSetup.hh"
 #include "RadmonApplicationGeneratorSetup.hh"

 // Forward declarations
 class RadmonApplicationOptions;
 class G4RunManager;
 class RadmonDetectorLayout;
 class RadmonGeneratorLayout;
 class RadmonDetectorLabelledEntitiesConstructorsFactory;
 class RadmonGeneratorsWithLabelFactory;
 class G4VisManager;
 class G4UImanager;
 class RadmonDetectorMessenger;
 class RadmonGeneratorMessenger;
 class G4UIsession;
 class G4UIdirectory;

 class RadmonApplication : public RadmonApplicationDetectorSetup, RadmonApplicationGeneratorSetup
 {
  public:
                                                RadmonApplication(const RadmonApplicationOptions & options);
                                               ~RadmonApplication();

   inline G4bool                                Valid(void) const;

  private:
   G4bool                                       RunMacro(const RadmonApplicationOptions & options, const char * fileName);
  
  // Hidden constructors and operators
                                                RadmonApplication();
                                                RadmonApplication(const RadmonApplication & copy);
   RadmonApplication &                          operator=(const RadmonApplication & copy);
   
  // Private attributes
   G4bool                                       valid;

   G4RunManager *                               runManager;
   RadmonDetectorLayout *                       detectorLayout;
   RadmonGeneratorLayout *                      generatorLayout;
   RadmonDetectorLabelledEntitiesConstructorsFactory * detectorsFactory;
   RadmonGeneratorsWithLabelFactory *           generatorsFactory;
   #ifdef    G4VIS_USE
    G4VisManager *                              visManager;
   #endif /* G4VIS_USE */
   G4UImanager *                                uiManager;
   RadmonDetectorMessenger *                    detectorMessenger;
   RadmonGeneratorMessenger *                   generatorMessenger;
   G4UIsession *                                session;
   G4UIdirectory *                              directory;
 };
 
 #include "RadmonApplication.icc"
#endif /* RADMONAPPLICATION_HH */
