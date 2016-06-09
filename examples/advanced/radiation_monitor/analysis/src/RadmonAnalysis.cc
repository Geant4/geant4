//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// File name:     RadmonAnalysis.cc
// Creation date: Nov 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonAnalysis.cc,v 1.2.2.2.4.1 2009/08/11 14:20:35 gcosmo Exp $
// Tag:           $Name: geant4-09-02-patch-03 $
//

// Include files
#include "RadmonAnalysis.hh"

#include "RadmonVDataAnalysisFactory.hh"
#include "RadmonVDataAnalysis.hh"
#include "RadmonVAnalysisLayout.hh"
#include "RadmonSensitiveDetector.hh"
#include "G4UIcommand.hh"
#include "G4Event.hh"
#include "G4SDManager.hh"
#include "G4RunManager.hh"
#include "G4Run.hh"
#include <fstream>
#include <exception>

using namespace std;

RadmonAnalysis :: RadmonAnalysis(RadmonVAnalysisLayout * layout, 
				 RadmonVDataAnalysisFactory * factory, 
				 AIDA::IAnalysisFactory * analysis)
  :
  analysisLayout(layout),
  dataFactory(factory),
  analysisFactory(analysis),
  treeFactory(0),
  tree(0),
  changed(true)
{
  if (analysisLayout==0)
    G4Exception("RadmonAnalysis::RadmonAnalysis: layout==0.");

  if (dataFactory==0)
    G4Exception("RadmonAnalysis::RadmonAnalysis: factory==0.");
 
  if (analysisFactory==0)
    G4Exception("RadmonAnalysis::RadmonAnalysis: analysis==0.");
 
  analysisLayout->AttachObserver(this);
 
  treeFactory=analysisFactory->createTreeFactory();
}



RadmonAnalysis :: ~RadmonAnalysis()
{
  analysisLayout->DetachObserver(this);
 
  Destruct();
 
  delete treeFactory;
  delete analysisFactory;
  delete dataFactory;
}

void RadmonAnalysis :: OnLayoutChange(void)
{
  changed=true;
}

void RadmonAnalysis :: OnBeginOfEvent(const G4Event * /* event */)
{
  if (!changed)
    return;

  Destruct();
 
  list<G4String> tupleLabels;
  G4String columns("int runId; int eventId");
 
  if (!InitializeSensitiveDetectorsList(tupleLabels, columns))
    {
      changed=true;
      return;
    }

  if (!OpenFile())
    {
      changed=true;
      return;
    }
 
  if (!InitializeTuple(tupleLabels, columns))
    {
      changed=false;
      return;
    }
 
  indexRunId=tuple->findColumn("runId");
  indexEventId=tuple->findColumn("eventId");
 
  changed=false;
}

void RadmonAnalysis :: OnEndOfEvent(const G4Event * event)
{
  if (!tuple)
    return;
  
  if (sensitiveDetectorsList.empty())
    return;

  SensitiveDetectorsList::iterator sd(sensitiveDetectorsList.begin());
  const SensitiveDetectorsList::iterator sdEnd(sensitiveDetectorsList.end());

  G4SDManager * manager(G4SDManager::GetSDMpointer());
 
  if (!manager)
    return;
  
  tuple->fill(indexRunId, G4RunManager::GetRunManager()->GetCurrentRun()->GetRunID());
  tuple->fill(indexEventId, event->GetEventID());
  
  while (sd!=sdEnd)
    {
      RadmonSensitiveDetector * radmonSensitiveDetector(sd->first);
      DataAnalysesList * dataAnalysesList(sd->second);

      RadmonHitsCollection * collection(radmonSensitiveDetector->GetDetectorCollection());

      DataAnalysesList::iterator i(dataAnalysesList->begin());
      DataAnalysesList::iterator end(dataAnalysesList->end());

      while (i!=end)
	{
	  (*i)->StoreIntoTuple(collection, tuple);
	  i++;
	}
  
      sd++;
    }
 
  tuple->addRow();
}

G4bool RadmonAnalysis :: InitializeSensitiveDetectorsList(list<G4String> & tupleLabels, G4String & columns)
{
  G4SDManager * manager(G4SDManager::GetSDMpointer());
 
  if (!manager)
    return false;

  RadmonVAnalysisLayout const * const constLayout(analysisLayout);
  G4int n(constLayout->GetNSensitiveDetectors());
 
  tupleLabels.clear();
 
  while (n>0)
    {
      n--;
  
      const G4String & sensitiveDetectorLabel(constLayout->GetSensitiveDetectorLabel(n));
      G4VSensitiveDetector * sensitiveDetector(manager->FindSensitiveDetector(sensitiveDetectorLabel, false));
  
      if (!sensitiveDetector)
	{
	  G4cout << "RadmonAnalysis::InitializeSensitiveDetectorsList: \"" << sensitiveDetectorLabel << "\" is not attached to any volume." << G4endl;
	  continue;
	}
   
      RadmonSensitiveDetector * radmonSensitiveDetector;
      try
	{
	  radmonSensitiveDetector=dynamic_cast<RadmonSensitiveDetector *>(sensitiveDetector);
	}
      catch (exception e)
	{
	  G4cout << "RadmonAnalysis::InitializeSensitiveDetectorsList: \"" << sensitiveDetectorLabel << "\" is not a RadmonSensitiveDetector object." << G4endl;
	  radmonSensitiveDetector=0;
	}
  
      if (!radmonSensitiveDetector)
	continue;
  
      const G4String & sensitiveDetectorType(constLayout->GetSensitiveDetectorType(sensitiveDetectorLabel));
  
      G4int m(constLayout->GetNDataAnalyses(sensitiveDetectorType));
  
      if (m==0)
	continue;
  
      DataAnalysesList * dataAnalysesList(new DataAnalysesList);
  
      G4String sensitiveDetectorColumns;
      while (m>0)
	{
	  m--;
   
	  const G4String & dataAnalysisLabel(constLayout->GetDataAnalysisLabel(sensitiveDetectorType, m));
	  const G4String & dataAnalysisType(constLayout->GetDataAnalysisType(sensitiveDetectorType, dataAnalysisLabel));
   
	  RadmonVDataAnalysis * dataAnalysis(dataFactory->CreateDataAnalysis(dataAnalysisType));
   
	  if (dataAnalysis==0)
	    {
	      G4cout << "RadmonAnalysis::InitializeSensitiveDetectorsList: \"" << dataAnalysisLabel << "\" defined in \"" << sensitiveDetectorType << "\" not found in the analysis factory." << G4endl;
	      continue;
	    }
   
	  G4String tupleLabel(sensitiveDetectorLabel);
	  tupleLabel+='_';
	  tupleLabel+=dataAnalysisLabel;
    
	  G4String dataAnalysisColumns(dataAnalysis->ObtainColumnsDeclaration(tupleLabel));
   
	  if (dataAnalysisColumns.empty())
	    continue;
    
	  dataAnalysesList->push_back(dataAnalysis);
	  tupleLabels.push_back(tupleLabel);
  
	  columns+="; ";
	  columns+=dataAnalysisColumns;
	}

      if (dataAnalysesList->empty())
	delete dataAnalysesList;
      else
	sensitiveDetectorsList.push_back(SensitiveDetectorPair(radmonSensitiveDetector, dataAnalysesList));
    }
 
  return !sensitiveDetectorsList.empty();
}

G4bool RadmonAnalysis :: OpenFile(void)
{
  RadmonVAnalysisLayout const * const constLayout(analysisLayout);
  const G4String & original(constLayout->GetOutputFileName());
  G4String fileName(original);
 
  for(G4int index(1); ; index++)
    {
      ifstream in(fileName);
  
      if (!in.is_open())
	break;
   
      fileName=(original+'_')+G4UIcommand::ConvertToString(index);
    }
 
  if (fileName!=original)
    G4cout << "RadmonAnalysis::OpenFile: \"" << original << "\" just exists. \"" << fileName << "\" will be used." << G4endl;

  try
    {
      tree=treeFactory->create(fileName, constLayout->GetOutputFileFormat(), false, true);
    }
  catch (exception e)
    {
      G4cout << "RadmonAnalysis::OpenFile: " << e.what() << G4endl;
      tree=0;
    }
 
  if (!tree)
    return false;
 
  tupleFactory=analysisFactory->createTupleFactory(* tree);
 
  if (!tupleFactory)
    {
      delete tree;
      tree=0;
      return false;
    }
 
  return true;
}

G4bool RadmonAnalysis :: InitializeTuple(const list<G4String> & tupleLabels, const G4String & columns)
{
  try
    {
      tuple=tupleFactory->create("1", "Radmon", columns);
    }
  catch (exception e)
    {
      G4cout << "RadmonAnalysis::InitializeTuple: " << e.what() << G4endl;
      tuple=0;
    }
  
  SensitiveDetectorsList::iterator sd(sensitiveDetectorsList.begin());
  const SensitiveDetectorsList::iterator sdEnd(sensitiveDetectorsList.end());
  list<G4String>::const_iterator j(tupleLabels.begin());

  while (sd!=sdEnd)
    {
      RadmonSensitiveDetector * radmonSensitiveDetector(sd->first);
      DataAnalysesList * dataAnalysesList(sd->second);

      DataAnalysesList::iterator i(dataAnalysesList->begin());
      DataAnalysesList::iterator end(dataAnalysesList->end());

      if (tuple==0)
	{
	  while (i!=end)
	    {
	      delete *i;
	      i=dataAnalysesList->erase(i);
	    }
	}
  
      if (dataAnalysesList->empty())
	{
	  delete dataAnalysesList;
	  sd=sensitiveDetectorsList.erase(sd);
	}
      else
	{
	  radmonSensitiveDetector->Activate(true);
   
	  i=dataAnalysesList->begin();
	  end=dataAnalysesList->end();
   
	  while (i!=end)
	    {
	      radmonSensitiveDetector->AttachDataStorer(*i);
	      (*i)->InitializeFromTuple(*j, tuple);
	      j++;
	      i++;
	    }
   
	  sd++;
	}
    }
 
  return !sensitiveDetectorsList.empty();
}

void RadmonAnalysis :: Destruct(void)
{
  if (tree)
    {
      tree->commit();
      tree->close();
     
      delete tupleFactory;
      delete tree;
      tree=0;
      tupleFactory=0;
    }
 
  while (!sensitiveDetectorsList.empty())
    {
      SensitiveDetectorPair sensitiveDetectorPair(sensitiveDetectorsList.back());
      sensitiveDetectorsList.pop_back();
      sensitiveDetectorPair.first->Activate(false);
      sensitiveDetectorPair.first->ClearDataStorersList();
  
      DataAnalysesList::iterator i(sensitiveDetectorPair.second->begin());  
      const DataAnalysesList::iterator end(sensitiveDetectorPair.second->end());  
  
      while (i!=end)
	{
	  delete (*i);
	  i++;
	}

      delete sensitiveDetectorPair.second;  
    }
}
