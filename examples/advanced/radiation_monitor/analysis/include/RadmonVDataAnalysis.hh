//
// File name:     RadmonVDataAnalysis.hh
// Creation date: Nov 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonVDataAnalysis.hh,v 1.2 2005-11-25 01:53:30 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   Abstract class of an analysis piece
//

#ifndef   RADMONVDATAANALYSIS_HH
 #define  RADMONVDATAANALYSIS_HH

 // Include files
 #include "RadmonSensitiveDetectorDataStorer.hh"
 #include "RadmonHit.hh"
 #include "G4String.hh"
 #include "AIDA/ITuple.h" 
 
 class RadmonVDataAnalysis : public RadmonSensitiveDetectorDataStorer
 {
  public:
   inline virtual                              ~RadmonVDataAnalysis();
    
   virtual void                                 SetDataAnalysisAttribute(const G4String & attributeName, const G4String &value) = 0;

   virtual G4String                             ObtainColumnsDeclaration(const G4String & prefix) = 0;
   virtual void                                 InitializeFromTuple(const G4String & prefix, const AIDA::ITuple * tuple) = 0;
   virtual void                                 StoreIntoTuple(RadmonHitsCollection * hitsCollection, AIDA::ITuple * tuple) = 0;
   
  protected:
   inline                                       RadmonVDataAnalysis();

  private:
  // Hidden constructors and operators
                                                RadmonVDataAnalysis(const RadmonVDataAnalysis & copy);
   RadmonVDataAnalysis &                        operator=(const RadmonVDataAnalysis & copy);
 };
 
 // Inline implementations
 #include "RadmonVDataAnalysis.icc"
#endif /* RADMONVDATAANALYSIS_HH */
