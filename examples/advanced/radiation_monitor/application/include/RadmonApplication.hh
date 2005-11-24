//
// File name:     RadmonApplication.hh
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonApplication.hh,v 1.6 2005-11-24 02:33:32 capra Exp $
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
 #include "RadmonApplicationPhysicsSetup.hh"
 #ifdef    G4ANALYSIS_USE
  #include "RadmonApplicationAnalysisSetup.hh"
 #endif /* G4ANALYSIS_USE */

 // Forward declarations
 class RadmonApplicationOptions;
 class G4RunManager;
 class RadmonDetectorLayout;
 class RadmonGeneratorLayout;
 class RadmonPhysicsLayout;
 class RadmonAnalysisLayout;
 class RadmonDetectorLabelledEntitiesConstructorsFactory;
 class RadmonGeneratorsWithLabelFactory;
 class RadmonSubPhysicsListWithLabelFactory;
 class RadmonDataAnalysisWithLabelFactory;
 class RadmonAnalysis;
 class G4VisManager;
 class G4UImanager;
 class RadmonDetectorMessenger;
 class RadmonGeneratorMessenger;
 class RadmonPhysicsMessenger;
 class RadmonAnalysisMessenger;
 class RadmonApplicationMessenger;
 class G4UIsession;
 class G4UIdirectory;

 class RadmonApplication : public RadmonApplicationDetectorSetup, public RadmonApplicationGeneratorSetup, public RadmonApplicationPhysicsSetup
 #ifdef    G4ANALYSIS_USE
                         , public RadmonApplicationAnalysisSetup
 #endif /* G4ANALYSIS_USE */
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
   RadmonPhysicsLayout *                        physicsLayout;
   #ifdef    G4ANALYSIS_USE
    RadmonAnalysisLayout *                      analysisLayout;
   #endif /* G4ANALYSIS_USE */
   RadmonDetectorLabelledEntitiesConstructorsFactory * detectorsFactory;
   RadmonGeneratorsWithLabelFactory *           generatorsFactory;
   RadmonSubPhysicsListWithLabelFactory *       physicsFactory;
   #ifdef    G4ANALYSIS_USE
    RadmonDataAnalysisWithLabelFactory *        analysisFactory;
    RadmonAnalysis *                            analysis;
   #endif /* G4ANALYSIS_USE */
   #ifdef    G4VIS_USE
    G4VisManager *                              visManager;
   #endif /* G4VIS_USE */
   G4UImanager *                                uiManager;
   RadmonDetectorMessenger *                    detectorMessenger;
   RadmonGeneratorMessenger *                   generatorMessenger;
   RadmonPhysicsMessenger *                     physicsMessenger;
   #ifdef    G4ANALYSIS_USE
    RadmonAnalysisMessenger *                   analysisMessenger;
   #endif /* G4ANALYSIS_USE */
   RadmonApplicationMessenger *                 applicationMessenger;
   G4UIsession *                                session;
   G4UIdirectory *                              directory;
 };
 
 #include "RadmonApplication.icc"
#endif /* RADMONAPPLICATION_HH */
