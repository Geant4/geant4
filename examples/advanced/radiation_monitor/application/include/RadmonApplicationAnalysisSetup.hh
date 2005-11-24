//
// File name:     RadmonApplicationAnalysisSetup.hh
// Creation date: Oct 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonApplicationAnalysisSetup.hh,v 1.1 2005-11-24 02:33:32 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   Radmon application analysis setup
//

#ifndef   RADMONAPPLICATIONANALYSISSETUP_HH
 #define  RADMONAPPLICATIONANALYSISSETUP_HH

 // Inglude files
 #include "globals.hh"

 // Forward declarations
 class RadmonApplicationOptions;
 class RadmonDataAnalysisWithLabelFactory;

 class RadmonApplicationAnalysisSetup
 {
  public:
   inline                                       RadmonApplicationAnalysisSetup(const RadmonApplicationOptions & options);
   inline                                      ~RadmonApplicationAnalysisSetup();

  protected:
   G4bool                                       CreateDataAnalysis(RadmonDataAnalysisWithLabelFactory * factory);

  // Hidden constructors and operators
                                                RadmonApplicationAnalysisSetup();
                                                RadmonApplicationAnalysisSetup(const RadmonApplicationAnalysisSetup & copy);
   RadmonApplicationAnalysisSetup &             operator=(const RadmonApplicationAnalysisSetup & copy);
   
  // Private attributes
   const RadmonApplicationOptions &             currentOptions;
 };
 
 #include "RadmonApplicationAnalysisSetup.icc"
#endif /* RADMONAPPLICATIONANALYSISSETUP_HH */
