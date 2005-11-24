//
// File name:     RadmonVDataAnalysisFactory.hh
// Creation date: Nov 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonVDataAnalysisFactory.hh,v 1.1 2005-11-24 02:33:32 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   Abstract class of a factory of analysis items
//

#ifndef   RADMONVDATAANALYSISFACTORY_HH
 #define  RADMONVDATAANALYSISFACTORY_HH

 // Forward declaration
 class RadmonVDataAnalysis;
 class G4String;

 class RadmonVDataAnalysisFactory
 {
  public:
   inline virtual                              ~RadmonVDataAnalysisFactory();

   virtual RadmonVDataAnalysis *                CreateDataAnalysis(const G4String & dataAnalysis) = 0;

  protected:
   inline                                       RadmonVDataAnalysisFactory();

  private:
  // Hidden constructors and operators
                                                RadmonVDataAnalysisFactory(const RadmonVDataAnalysisFactory & copy);
   RadmonVDataAnalysisFactory &                 operator=(const RadmonVDataAnalysisFactory & copy);
 };
 
 // Inline implementations
 #include "RadmonVDataAnalysisFactory.icc"
#endif /* RADMONVDATAANALYSISFACTORY_HH */
