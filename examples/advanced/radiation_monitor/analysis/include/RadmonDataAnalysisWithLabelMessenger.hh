//
// File name:     RadmonDataAnalysisWithLabelMessenger.hh
// Creation date: Nov 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonDataAnalysisWithLabelMessenger.hh,v 1.1 2005-11-24 02:33:32 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   UI commands for dumping the available constructors
//

#ifndef   RADMONDATAANALYSISWITHLABELMESSENGER_HH
 #define  RADMONDATAANALYSISWITHLABELMESSENGER_HH

 // Include files
 #include "RadmonMessenger.hh"
 #include "G4String.hh"
 #include <set>

 class RadmonDataAnalysisWithLabelMessenger : public RadmonMessenger
 {
  public:
   static RadmonDataAnalysisWithLabelMessenger * Instance(void);
  
   void                                         AddAvailableDataAnalysis(const G4String & name);
   void                                         RemoveAvailableDataAnalysis(const G4String & name);
  
   virtual G4String                             GetCurrentValue(G4UIcommand * command);
   virtual void                                 SetNewValue(G4UIcommand * command, G4String newValue);

  private:
  // Hidden constructors and operators
                                                RadmonDataAnalysisWithLabelMessenger();
                                                RadmonDataAnalysisWithLabelMessenger(const RadmonDataAnalysisWithLabelMessenger & copy);
                                               ~RadmonDataAnalysisWithLabelMessenger();
   RadmonDataAnalysisWithLabelMessenger &       operator=(const RadmonDataAnalysisWithLabelMessenger & copy);

  // Private Data Types
   typedef std::set<G4String>                   AvailableDataAnalyses;
   
  // Private variables
   AvailableDataAnalyses                        availableDataAnalyses;
   
   static RadmonDataAnalysisWithLabelMessenger * instance;
   
  // Commands
   RADMON_DECLARE_COMMAND(Dump);
 };
#endif /* RADMONDATAANALYSISWITHLABELMESSENGER_HH */
