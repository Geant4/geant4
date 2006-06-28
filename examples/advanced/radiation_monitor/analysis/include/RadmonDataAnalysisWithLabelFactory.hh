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
// File name:     RadmonDataAnalysisWithLabelFactory.hh
// Creation date: Nov 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonDataAnalysisWithLabelFactory.hh,v 1.2 2006-06-28 13:44:05 gunter Exp $
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
