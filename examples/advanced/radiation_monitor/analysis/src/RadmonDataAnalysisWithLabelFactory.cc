//
// File name:     RadmonDataAnalysisWithLabelFactory.cc
// Creation date: Nov 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonDataAnalysisWithLabelFactory.cc,v 1.1 2005-11-24 02:33:32 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//

// Include files
#include "RadmonDataAnalysisWithLabelFactory.hh"
#include "RadmonDataAnalysisWithLabelMessenger.hh"
#include "RadmonVDataAnalysisWithLabel.hh"


                                                RadmonDataAnalysisWithLabelFactory :: ~RadmonDataAnalysisWithLabelFactory()
{
 DataAnalyses::iterator i(dataAnalyses.begin());
 DataAnalyses::iterator end(dataAnalyses.end());
 
 RadmonDataAnalysisWithLabelMessenger * messenger(RadmonDataAnalysisWithLabelMessenger::Instance());
 
 while (i!=end)
 {
  messenger->RemoveAvailableDataAnalysis((*i)->GetLabel());

  delete (*i);
  i++;
 }
}





RadmonVDataAnalysis *                         RadmonDataAnalysisWithLabelFactory :: CreateDataAnalysis(const G4String & dataAnalysisName)
{
 DataAnalyses::iterator i(dataAnalyses.begin());
 DataAnalyses::iterator end(dataAnalyses.end());
 
 while (i!=end)
 {
  if ((*i)->GetLabel()==dataAnalysisName)
   return (*i)->New();

  i++;
 }
 
 return 0;
}



void                                            RadmonDataAnalysisWithLabelFactory :: AppendDataAnalysisWithLabel(RadmonVDataAnalysisWithLabel * dataAnalysis)
{
 dataAnalyses.push_back(dataAnalysis);
 
 RadmonDataAnalysisWithLabelMessenger * messenger(RadmonDataAnalysisWithLabelMessenger::Instance());
 messenger->AddAvailableDataAnalysis(dataAnalysis->GetLabel());
}
