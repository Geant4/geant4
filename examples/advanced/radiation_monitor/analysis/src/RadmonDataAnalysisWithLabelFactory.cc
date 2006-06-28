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
// File name:     RadmonDataAnalysisWithLabelFactory.cc
// Creation date: Nov 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonDataAnalysisWithLabelFactory.cc,v 1.2 2006-06-28 13:44:39 gunter Exp $
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
