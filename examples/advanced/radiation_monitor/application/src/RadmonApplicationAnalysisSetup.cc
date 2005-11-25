//
// File name:     RadmonApplicationAnalysisSetup.cc
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonApplicationAnalysisSetup.cc,v 1.2 2005-11-25 01:56:26 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//

// Include files
#include "RadmonApplicationAnalysisSetup.hh"
#include "RadmonApplicationOptions.hh"

#include "RadmonDataAnalysisDepositedEnergy.hh"

#include "RadmonDataAnalysisWithLabelFactory.hh"


#define DECLARE_ANALYSIS_CONSTRUCTOR(name)      constructor=new name();                                                                  \
                                                if (constructor==0)                                                                      \
                                                {                                                                                        \
                                                 G4cerr << currentOptions.ApplicationName() << ": Cannot allocate " #name "." << G4endl; \
                                                 return false;                                                                           \
                                                }                                                                                        \
                                                factory->AppendDataAnalysisWithLabel(constructor)

G4bool                                          RadmonApplicationAnalysisSetup :: CreateDataAnalysis(RadmonDataAnalysisWithLabelFactory * factory)
{
 RadmonVDataAnalysisWithLabel * constructor;

 DECLARE_ANALYSIS_CONSTRUCTOR(RadmonDataAnalysisDepositedEnergy);

 return true;
}
