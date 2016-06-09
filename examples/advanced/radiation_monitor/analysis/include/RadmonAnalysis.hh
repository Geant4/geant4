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
// File name:     RadmonAnalysis.hh
// Creation date: Nov 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonAnalysis.hh,v 1.4 2006/06/29 16:06:41 gunter Exp $
// Tag:           $Name: geant4-09-02 $
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
