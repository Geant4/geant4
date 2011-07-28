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
// $Id$

// Author: Ivana Hrivnacova, 15/06/2011  (ivana@ipno.in2p3.fr)

#include "G4CsvAnalysisManager.hh"
#include "G4UnitsTable.hh"

#include "tools/waxml/begend"
#include "tools/waxml/histos"

#include <iostream>

//_____________________________________________________________________________
G4CsvAnalysisManager::G4CsvAnalysisManager()
 : fFile(0),
   fNtupleName(),
   fNtupleTitle(),
   fNtuple(0),
   fNtupleIColumnMap(),
   fNtupleFColumnMap(),
   fNtupleDColumnMap()
{
}

//_____________________________________________________________________________
G4CsvAnalysisManager::~G4CsvAnalysisManager()
{  
  delete fNtuple;
  delete fFile;
}

// 
// private methods
//

//_____________________________________________________________________________
tools::wcsv::ntuple::column<int>*    
G4CsvAnalysisManager::GetNtupleIColumn(G4int id) const
{
  std::map<G4int, tools::wcsv::ntuple::column<int>* >::const_iterator it
    = fNtupleIColumnMap.find(id);
  if ( it == fNtupleIColumnMap.end() ) {
    G4cerr << "---> Warning from G4CsvAnalysisManager::GetNtupleIColumn():"
           << " column " << id << " does not exist. " << G4endl;
    return 0;
  }
  
  return it->second;
}  
    
//_____________________________________________________________________________
tools::wcsv::ntuple::column<float>*  
G4CsvAnalysisManager::GetNtupleFColumn(G4int id) const
{
  std::map<G4int, tools::wcsv::ntuple::column<float>* >::const_iterator it
    = fNtupleFColumnMap.find(id);
  if ( it == fNtupleFColumnMap.end() ) {
    G4cerr << "---> Warning from G4CsvAnalysisManager::GetNtupleFColumn():"
           << " column " << id << " does not exist. " << G4endl;
    return 0;
  }
  
  return it->second;
}  


//_____________________________________________________________________________
tools::wcsv::ntuple::column<double>* 
G4CsvAnalysisManager::GetNtupleDColumn(G4int id) const
{
  std::map<G4int, tools::wcsv::ntuple::column<double>* >::const_iterator it
    = fNtupleDColumnMap.find(id);
  if ( it == fNtupleDColumnMap.end() ) {
    G4cerr << "---> Warning from G4CsvAnalysisManager::GetNtupleDColumn():"
           << " column " << id << " does not exist. " << G4endl;
    return 0;
  }
  
  return it->second;
}  
 
// 
// public methods
//

//_____________________________________________________________________________
G4bool G4CsvAnalysisManager::OpenFile(const G4String& fileName)
{
  // Add file extension .root if no extension is given
  G4String name(fileName);
  if ( name.find(".") == std::string::npos ) name.append(".csv");

  fFile = new std::ofstream(name);
  if ( fFile->fail() ) {
    G4cerr << "---> Warning from G4CsvAnalysisManager::OpenFile(): "
           << "Cannot open file " << name << std::endl;
    return false;
  }

  return true;
}  
  
//_____________________________________________________________________________
G4bool G4CsvAnalysisManager::Write() 
{
  // nothing to be done here
  return true;
}

//_____________________________________________________________________________
G4bool G4CsvAnalysisManager::CloseFile()
{
  fFile->close(); 
  return true; 
} 
   
//_____________________________________________________________________________
G4int G4CsvAnalysisManager::CreateH1(const G4String& /*name*/, 
                               const G4String& /*title*/, G4int /*nbins*/, 
                               G4double /*xmin*/, G4double /*xmax*/)
{
  G4cerr << "---> Warning from  G4CsvAnalysisManager::CreateH1():  "
         << " histograms are not supported." << G4endl;
  return 0;
}                                         

//_____________________________________________________________________________
void G4CsvAnalysisManager::CreateNtuple(const G4String& name, 
                                        const G4String& title)
{
  if ( fNtuple ) {
    G4cerr << "---> Warning from G4CsvAnalysisManager::CreateNtuple():"
           << " Ntuple already exists. (Only one ntuple is currently supported.)" 
           << G4endl;
    return;       
  }

  fNtuple = new tools::wcsv::ntuple(*fFile);
  fNtupleName = name;
  fNtupleTitle = title;
}                                         

//_____________________________________________________________________________
G4int G4CsvAnalysisManager::CreateNtupleIColumn(const G4String& name)
{
  G4int index = fNtuple->columns().size();
  tools::wcsv::ntuple::column<int>* column = fNtuple->create_column<int>(name);  
  fNtupleIColumnMap[index] = column;
  return index + fFirstNtupleId;       
}                                         

//_____________________________________________________________________________
G4int G4CsvAnalysisManager::CreateNtupleFColumn(const G4String& name)
{
  G4int index = fNtuple->columns().size();
  tools::wcsv::ntuple::column<float>* column = fNtuple->create_column<float>(name);  
  fNtupleFColumnMap[index] = column;
  return index + fFirstNtupleId;       
}                                         

//_____________________________________________________________________________
G4int G4CsvAnalysisManager::CreateNtupleDColumn(const G4String& name)   
{
  G4int index = fNtuple->columns().size();
  tools::wcsv::ntuple::column<double>* column = fNtuple->create_column<double>(name);  
  fNtupleDColumnMap[index] = column;
  return index + fFirstNtupleId;       
}                                         

//_____________________________________________________________________________
void G4CsvAnalysisManager::FinishNtuple()
{ 
  // nothing to be done here
}
   
  
//_____________________________________________________________________________
G4bool G4CsvAnalysisManager::FillH1(G4int /*id*/, 
                                    G4double /*value*/, G4double /*weight*/)
{
  G4cerr << "---> Warning from  G4CsvAnalysisManager::FillH1():  "
         << " histograms are not supported." << G4endl;
  return false;
}

//_____________________________________________________________________________
G4bool G4CsvAnalysisManager::FillNtupleIColumn(G4int id, G4int value)
{
  tools::wcsv::ntuple::column<int>* column = GetNtupleIColumn(id);
  if ( ! column ) {
    G4cout << "---> Warning from G4CsvAnalysisManager::FillNtupleIColumn():"
           << " column " << id << " does not exist. " << G4endl;
    return false;
  }  
  
  column->fill(value);
  return true;       
}                                         
//_____________________________________________________________________________
G4bool G4CsvAnalysisManager::FillNtupleFColumn(G4int id, G4float value)
{
  tools::wcsv::ntuple::column<float>* column = GetNtupleFColumn(id);
  if ( ! column ) {
    G4cout << "---> Warning from G4CsvAnalysisManager::FillNtupleFColumn():"
           << " column " << id << " does not exist. " << G4endl;
    return false;
  }  
  
  column->fill(value);
  return true;       
}                                         

//_____________________________________________________________________________
G4bool G4CsvAnalysisManager::FillNtupleDColumn(G4int id, G4double value)
{
  tools::wcsv::ntuple::column<double>* column = GetNtupleDColumn(id);
  if ( ! column ) {
    G4cout << "---> Warning from G4CsvAnalysisManager::FillNtupleDColumn():"
           << " column " << id << " does not exist. " << G4endl;
    return false;
  }  
  
  column->fill(value);
  return true;       
}                                         

//_____________________________________________________________________________
G4bool G4CsvAnalysisManager::AddNtupleRow()
{ 
  if ( ! fNtuple ) {
    G4cout << "---> Warning from G4CsvAnalysisManager::AddNtupleRow(): " 
           << " ntuple does not exist. " << G4endl;
    return false;
  }  
  
  fNtuple->add_row();
  return true;
}
