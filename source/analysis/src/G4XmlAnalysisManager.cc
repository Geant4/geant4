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

#include "G4XmlAnalysisManager.hh"
#include "G4UnitsTable.hh"

#include "tools/waxml/begend"
#include "tools/waxml/histos"

#include <iostream>

//_____________________________________________________________________________
G4XmlAnalysisManager::G4XmlAnalysisManager()
 : fFile(0),
   fH1Vector(),   
   fH1MapByName(),
   fNtupleName(),
   fNtupleTitle(),
   fNtuple(0),
   fNtupleIColumnMap(),
   fNtupleFColumnMap(),
   fNtupleDColumnMap()
{
}

//_____________________________________________________________________________
G4XmlAnalysisManager::~G4XmlAnalysisManager()
{  
  std::vector<tools::histo::h1d*>::iterator it;
  for ( it = fH1Vector.begin(); it != fH1Vector.end(); it++ ) {
    delete *it;
  }  
  delete fNtuple;
  delete fFile;  
}

// 
// private methods
//

//_____________________________________________________________________________
tools::waxml::ntuple::column<int>*    
G4XmlAnalysisManager::GetNtupleIColumn(G4int id) const
{
  std::map<G4int, tools::waxml::ntuple::column<int>* >::const_iterator it
    = fNtupleIColumnMap.find(id);
  if ( it == fNtupleIColumnMap.end() ) {
    G4cerr << "---> Warning from G4XmlAnalysisManager::GetNtupleIColumn():"
           << " column " << id << " does not exist. " << G4endl;
    return 0;
  }
  
  return it->second;
}  
    
//_____________________________________________________________________________
tools::waxml::ntuple::column<float>*  
G4XmlAnalysisManager::GetNtupleFColumn(G4int id) const
{
  std::map<G4int, tools::waxml::ntuple::column<float>* >::const_iterator it
    = fNtupleFColumnMap.find(id);
  if ( it == fNtupleFColumnMap.end() ) {
    G4cerr << "---> Warning from G4XmlAnalysisManager::GetNtupleFColumn():"
           << " column " << id << " does not exist. " << G4endl;
    return 0;
  }
  
  return it->second;
}  


//_____________________________________________________________________________
tools::waxml::ntuple::column<double>* 
G4XmlAnalysisManager::GetNtupleDColumn(G4int id) const
{
  std::map<G4int, tools::waxml::ntuple::column<double>* >::const_iterator it
    = fNtupleDColumnMap.find(id);
  if ( it == fNtupleDColumnMap.end() ) {
    G4cerr << "---> Warning from G4XmlAnalysisManager::GetNtupleDColumn():"
           << " column " << id << " does not exist. " << G4endl;
    return 0;
  }
  
  return it->second;
}  
 
// 
// public methods
//

//_____________________________________________________________________________
G4bool G4XmlAnalysisManager::OpenFile(const G4String& fileName)
{
  // Add file extension .root if no extension is given
  G4String name(fileName);
  if ( name.find(".") == std::string::npos ) name.append(".xml");

  fFile = new std::ofstream(name);
  if ( fFile->fail() ) {
    G4cerr << "---> Warning from G4XmlAnalysisManager::OpenFile(): "
           << "Cannot open file " << name << std::endl;
    return false;
  }

  tools::waxml::begin(*fFile);
  return true;
}  
  
//_____________________________________________________________________________
G4bool G4XmlAnalysisManager::Write() 
{
  // ntuple 
  if ( fNtuple ) fNtuple->write_trailer();

  // histograms
  std::map<G4String, tools::histo::h1d*>::iterator it;
  for ( it = fH1MapByName.begin(); it != fH1MapByName.end(); it++ ) {
    G4bool result
      = tools::waxml::write(*fFile, *(it->second), fHistoDirectoryName, it->first);
    if ( ! result ) {
      G4cerr << "---> Warning from  G4XmlAnalysisManager::Write():"
             << " saving histo " << it->first << " failed" << G4endl;
      return false;       
    } 
  }

  return true;
}

//_____________________________________________________________________________
G4bool G4XmlAnalysisManager::CloseFile()
{
  tools::waxml::end(*fFile);
  fFile->close(); 
  return true; 
} 
   
//_____________________________________________________________________________
G4int G4XmlAnalysisManager::CreateH1(const G4String& name, const G4String& title, 
                               G4int nbins, G4double xmin, G4double xmax)
{
  G4int index = fH1Vector.size();
  tools::histo::h1d* h1 = new tools::histo::h1d(title, nbins, xmin, xmax);
  fH1Vector.push_back(h1);
  fH1MapByName[name] = h1;
  return index + fFirstHistoId;
}                                         

//_____________________________________________________________________________
void G4XmlAnalysisManager::CreateNtuple(const G4String& name, 
                                        const G4String& title)
{
  if ( fNtuple ) {
    G4cerr << "---> Warning from G4XmlAnalysisManager::CreateNtuple():"
           << " Ntuple already exists. (Only one ntuple is currently supported.)" 
           << G4endl;
    return;       
  }

  fNtuple = new tools::waxml::ntuple(*fFile);
  fNtupleName = name;
  fNtupleTitle = title;
}                                         

//_____________________________________________________________________________
G4int G4XmlAnalysisManager::CreateNtupleIColumn(const G4String& name)
{
  G4int index = fNtuple->columns().size();
  tools::waxml::ntuple::column<int>* column = fNtuple->create_column<int>(name);  
  fNtupleIColumnMap[index] = column;
  return index + fFirstNtupleId;       
}                                         

//_____________________________________________________________________________
G4int G4XmlAnalysisManager::CreateNtupleFColumn(const G4String& name)
{
  G4int index = fNtuple->columns().size();
  tools::waxml::ntuple::column<float>* column = fNtuple->create_column<float>(name);  
  fNtupleFColumnMap[index] = column;
  return index + fFirstNtupleId;       
}                                         

//_____________________________________________________________________________
G4int G4XmlAnalysisManager::CreateNtupleDColumn(const G4String& name)   
{
  G4int index = fNtuple->columns().size();
  tools::waxml::ntuple::column<double>* column = fNtuple->create_column<double>(name);  
  fNtupleDColumnMap[index] = column;
  return index + fFirstNtupleId;       
}                                         

//_____________________________________________________________________________
void G4XmlAnalysisManager::FinishNtuple()
{ 
  G4String path = "/";
  path.append(fNtupleDirectoryName);
  fNtuple->write_header(path, fNtupleName, fNtupleTitle);  
}
   
  
//_____________________________________________________________________________
G4bool G4XmlAnalysisManager::FillH1(G4int id, G4double value, G4double weight)
{
  tools::histo::h1d* h1d = GetH1(id, false);
  if ( ! h1d ) {
    G4cerr << "---> Warning from G4XmlAnalysisManager::FillH1():"
           << " histo " << id << " does not exist. " << G4endl;
    return false;
  }  

  h1d->fill(value, weight);
  return true;
}

//_____________________________________________________________________________
G4bool G4XmlAnalysisManager::FillNtupleIColumn(G4int id, G4int value)
{
  tools::waxml::ntuple::column<int>* column = GetNtupleIColumn(id);
  if ( ! column ) {
    G4cerr << "---> Warning from G4XmlAnalysisManager::FillNtupleIColumn():"
           << " column " << id << " does not exist. " << G4endl;
    return false;
  }  
  
  column->fill(value);
  return true;       
}                                         
//_____________________________________________________________________________
G4bool G4XmlAnalysisManager::FillNtupleFColumn(G4int id, G4float value)
{
  tools::waxml::ntuple::column<float>* column = GetNtupleFColumn(id);
  if ( ! column ) {
    G4cerr << "---> Warning from G4XmlAnalysisManager::FillNtupleFColumn():"
           << " column " << id << " does not exist. " << G4endl;
    return false;
  }  
  
  column->fill(value);
  return true;       
}                                         

//_____________________________________________________________________________
G4bool G4XmlAnalysisManager::FillNtupleDColumn(G4int id, G4double value)
{
  tools::waxml::ntuple::column<double>* column = GetNtupleDColumn(id);
  if ( ! column ) {
    G4cerr << "---> Warning from G4XmlAnalysisManager::FillNtupleDColumn():"
           << " column " << id << " does not exist. " << G4endl;
    return false;
  }  
  
  column->fill(value);
  return true;       
}                                         

//_____________________________________________________________________________
G4bool G4XmlAnalysisManager::AddNtupleRow()
{ 
  if ( ! fNtuple ) {
    G4cerr << "---> Warning from G4XmlAnalysisManager::AddNtupleRow(): " 
           << " ntuple does not exist. " << G4endl;
    return false;
  }  
  
  fNtuple->add_row();
  return true;
}
 
//_____________________________________________________________________________
tools::histo::h1d*  G4XmlAnalysisManager::GetH1(G4int id, G4bool warn) const 
{
  G4int index = id - fFirstHistoId;
  if ( index < 0 || index >= G4int(fH1Vector.size()) ) {
    if ( warn) {
      G4cerr << "---> Warning from G4XmlAnalysisManager::GetH1():"
             << " histo " << id << " does not exist. " << G4endl;
    }
    return 0;         
  }
  return fH1Vector[index];
}
