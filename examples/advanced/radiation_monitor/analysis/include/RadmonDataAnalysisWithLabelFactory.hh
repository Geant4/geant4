//
// File name:     RadmonDataAnalysisWithLabelFactory.hh
// Creation date: Nov 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonDataAnalysisWithLabelFactory.hh,v 1.1 2005-11-24 02:33:32 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   Concrete factory that constructs 
//                RadmonSubPhysiscsListWithLabel objects
//

#ifndef   RADMONDATAANALYSISWITHLABELFACTORY_HH
 #define  RADMONDATAANALYSISWITHLABELFACTORY_HH

 // Include files
 #include "RadmonVDataAnalysisFactory.hh"
 #include <list>

 // Forward declarations
 class RadmonVDataAnalysis;
 class RadmonVDataAnalysisWithLabel;

 class RadmonDataAnalysisWithLabelFactory : public RadmonVDataAnalysisFactory
 {
  public:
   inline                                       RadmonDataAnalysisWithLabelFactory();
   virtual                                     ~RadmonDataAnalysisWithLabelFactory();

   virtual RadmonVDataAnalysis *                CreateDataAnalysis(const G4String & dataAnalysisName);

   void                                         AppendDataAnalysisWithLabel(RadmonVDataAnalysisWithLabel * physicsList);

  private:
  // Hidden constructors and operators
                                                RadmonDataAnalysisWithLabelFactory(const RadmonDataAnalysisWithLabelFactory & copy);
   RadmonDataAnalysisWithLabelFactory &         operator=(const RadmonDataAnalysisWithLabelFactory & copy);

  // Private attributes
   typedef std::list<RadmonVDataAnalysisWithLabel *> DataAnalyses;
   DataAnalyses                                 dataAnalyses;
 };
 
 // Inline implementations
 #include "RadmonDataAnalysisWithLabelFactory.icc"
#endif /* RADMONDATAANALYSISWITHLABELFACTORY_HH */
