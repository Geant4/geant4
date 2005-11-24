//
// File name:     RadmonAnalysis.hh
// Creation date: Nov 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonAnalysis.hh,v 1.1 2005-11-24 02:33:32 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   Analysis main class
//

#ifndef   RADMONANALYSIS_HH
 #define  RADMONANALYSIS_HH

 // Include files
 #include "globals.hh"
 #include "RadmonEventActionObserver.hh"
 #include "RadmonVLayoutObserver.hh"
 #include "AIDA/IAnalysisFactory.h"
 #include "AIDA/ITreeFactory.h"
 #include "AIDA/ITree.h"
 #include <list>
 
 // Forward declaration
 class RadmonVDataAnalysisFactory;
 class RadmonVDataAnalysis;
 class RadmonVAnalysisLayout;

 class RadmonAnalysis : public RadmonEventActionObserver, public RadmonVLayoutObserver
 {
  public:
                                                RadmonAnalysis(RadmonVAnalysisLayout * layout, RadmonVDataAnalysisFactory * factory, AIDA::IAnalysisFactory * analysis);
   virtual                                     ~RadmonAnalysis();

   virtual void                                 OnLayoutChange(void);

   virtual void                                 OnBeginOfEvent(const G4Event * event);
   virtual void                                 OnEndOfEvent(const G4Event * event);

  private:
  // Private methods
   void                                         Destruct(void);

  // Hidden constructors and operators
                                                RadmonAnalysis();
                                                RadmonAnalysis(const RadmonAnalysis & copy);
   RadmonAnalysis &                             operator=(const RadmonAnalysis & copy);

  // Private data types
   typedef std::list<RadmonVDataAnalysis *>     DataAnalysesList;

  // Private attributes
   RadmonVAnalysisLayout *                      analysisLayout;
   RadmonVDataAnalysisFactory *                 dataFactory;
   
   AIDA::IAnalysisFactory *                     analysisFactory;
   AIDA::ITreeFactory *                         treeFactory;
   AIDA::ITree *                                tree;
  
   G4bool                                       changed;
  
   DataAnalysesList                             dataAnalysesList;
 };
#endif /* RADMONANALYSIS_HH */
