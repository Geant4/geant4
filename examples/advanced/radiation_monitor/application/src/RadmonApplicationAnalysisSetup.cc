//
// File name:     RadmonApplicationAnalysisSetup.cc
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonApplicationAnalysisSetup.cc,v 1.1 2005-11-24 02:34:48 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//

// Include files
#include "RadmonApplicationAnalysisSetup.hh"
#include "RadmonApplicationOptions.hh"

#include "RadmonDataAnalysisWithLabelFactory.hh"


#define DECLARE_ANALYSIS_CONSTRUCTOR(name)      constructor=new name();                                                                  \
                                                if (constructor==0)                                                                      \
                                                {                                                                                        \
                                                 G4cerr << currentOptions.ApplicationName() << ": Cannot allocate " #name "." << G4endl; \
                                                 return false;                                                                           \
                                                }                                                                                        \
                                                factory->AppendAnalysis(constructor)

G4bool                                          RadmonApplicationAnalysisSetup :: CreateDataAnalysis(RadmonDataAnalysisWithLabelFactory * /* factory */)
{
// RadmonVAnalysisWithLabel * constructor;

 return true;
}
