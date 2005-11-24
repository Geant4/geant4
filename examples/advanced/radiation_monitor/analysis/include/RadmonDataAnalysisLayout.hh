//
// File name:     RadmonDataAnalysisLayout.hh
// Creation date: Nov 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonDataAnalysisLayout.hh,v 1.1 2005-11-24 02:33:32 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   Internal class that describes layer informations
//

#ifndef   RADMONDATAANALYSISLAYOUT_HH
 #define  RADMONDATAANALYSISLAYOUT_HH

 // Include files
 #include "RadmonLayoutEntityWithAttributes.hh"

 class RadmonDataAnalysisLayout : public RadmonLayoutEntityWithAttributes
 {
  public:
   inline                                       RadmonDataAnalysisLayout();
   inline                                       RadmonDataAnalysisLayout(const RadmonDataAnalysisLayout & copy);
   inline                                      ~RadmonDataAnalysisLayout();

   inline RadmonDataAnalysisLayout &            operator=(const RadmonDataAnalysisLayout & copy);

   inline const G4String &                      GetLabel(void) const;
   inline void                                  SetLabel(const G4String & label);

   inline const G4String &                      GetType(void) const;
   inline void                                  SetType(const G4String & type);

   void                                         DumpLayout(std::ostream & out, const G4String &indent=G4String()) const;

  private:
  // Private attributes
   G4String                                     dataAnalysisLabel;
   G4String                                     dataAnalysisType;
 };

 // Inline implementations
 #include "RadmonDataAnalysisLayout.icc"
#endif /* RADMONDATAANALYSISLAYOUT_HH */
