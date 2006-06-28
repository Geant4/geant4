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
// File name:     RadmonAnalysis.hh
// Creation date: Nov 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonAnalysis.hh,v 1.3 2006-06-28 13:43:47 gunter Exp $
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
 #include "AIDA/ITupleFactory.h"
 #include "AIDA/ITuple.h"
 #include <utility>
 #include <list>
 
 // Forward declaration
 class RadmonVDataAnalysisFactory;
 class RadmonVDataAnalysis;
 class RadmonVAnalysisLayout;
 class RadmonSensitiveDetector;

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
   G4bool                                       InitializeSensitiveDetectorsList(std::list<G4String> & tupleLabels, G4String & columns);
   G4bool                                       OpenFile(void);
   G4bool                                       InitializeTuple(const std::list<G4String> & tupleLabels, const G4String & columns);

   void                                         Destruct(void);
   
  // Hidden constructors and operators
                                                RadmonAnalysis();
                                                RadmonAnalysis(const RadmonAnalysis & copy);
   RadmonAnalysis &                             operator=(const RadmonAnalysis & copy);

  // Private data types
   typedef std::list<RadmonVDataAnalysis *>     DataAnalysesList;
   typedef std::pair<RadmonSensitiveDetector *, DataAnalysesList *> SensitiveDetectorPair;
   typedef std::list<SensitiveDetectorPair>     SensitiveDetectorsList;

  // Private attributes
   RadmonVAnalysisLayout *                      analysisLayout;
   RadmonVDataAnalysisFactory *                 dataFactory;
   
   AIDA::IAnalysisFactory *                     analysisFactory;
   AIDA::ITreeFactory *                         treeFactory;
   AIDA::ITree *                                tree;
   AIDA::ITupleFactory *                        tupleFactory;
   AIDA::ITuple *                               tuple;
   
   G4int                                        indexRunId;
   G4int                                        indexEventId;
  
   G4bool                                       changed;
  
   SensitiveDetectorsList                       sensitiveDetectorsList;
 };
#endif /* RADMONANALYSIS_HH */
