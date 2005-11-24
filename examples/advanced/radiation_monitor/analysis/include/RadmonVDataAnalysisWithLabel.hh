//
// File name:     RadmonVDataAnalysisWithLabel.hh
// Creation date: Nov 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonVDataAnalysisWithLabel.hh,v 1.1 2005-11-24 02:33:32 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   Abstract class of an analysis piece with label
//                and attributes
//

#ifndef   RADMONVDATAANALYSISWITHLABEL_HH
 #define  RADMONVDATAANALYSISWITHLABEL_HH
 
 // Include files
 #include "RadmonLayoutEntityWithAttributes.hh"
 #include "RadmonVDataAnalysis.hh"
 
 class RadmonVDataAnalysisWithLabel : public RadmonVDataAnalysis, public RadmonLayoutEntityWithAttributes
 {
  public:
   inline virtual                              ~RadmonVDataAnalysisWithLabel();

   inline const G4String &                      GetLabel(void) const;
   inline virtual void                          SetDataAnalysisAttribute(const G4String & attributeName, const G4String & value);

   virtual RadmonVDataAnalysisWithLabel *       New(void) const = 0;

  protected:
   inline                                       RadmonVDataAnalysisWithLabel(const G4String & label);
   
  private:
  // Hidden constructors and operators
                                                RadmonVDataAnalysisWithLabel();
                                                RadmonVDataAnalysisWithLabel(const RadmonVDataAnalysisWithLabel & copy);
   RadmonVDataAnalysisWithLabel &               operator=(const RadmonVDataAnalysisWithLabel & copy);

  // Private attributes
   G4String                                     physiscListLabel;
 };
 
 #include "RadmonVDataAnalysisWithLabel.icc"
#endif /* RADMONVDATAANALYSISWITHLABEL_HH */
