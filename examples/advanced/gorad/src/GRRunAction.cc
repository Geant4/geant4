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
//  Gorad (Geant4 Open-source Radiation Analysis and Design)
//
//  Author : Makoto Asai (SLAC National Accelerator Laboratory)
//
//  Development of Gorad is funded by NASA Johnson Space Center (JSC)
//  under the contract NNJ15HK11B.
//
// ********************************************************************
//
// GRRunAction.cc
//   Gorad Run Action class that takes care of defining and handling
//   histograms and n-tuple.
//   Filling histograms is taken care by GRRun class.
//
// History
//   September 8th, 2020 : first implementation
//
// ********************************************************************

#include "GRRunAction.hh"
#include "GRRunActionMessenger.hh"

#include "G4AnalysisManager.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

GRRunAction::GRRunAction()
:fileName("GoradOut"), fileOpen(false), verbose(0), ifCarry(false),
 id_offset(100), id_factor(100)
{ 
  messenger = new GRRunActionMessenger(this);
  auto analysisManager = G4AnalysisManager::Instance();
  analysisManager->SetDefaultFileType("csv");
  G4cout << "Using " << analysisManager->GetType() << G4endl;

  analysisManager->SetVerboseLevel(verbose);
  //analysisManager->SetNtupleMerging(true);
}

GRRunAction::~GRRunAction()
{
  delete messenger;
}

void GRRunAction::BeginOfRunAction(const G4Run* /*run*/)
{ 
  // Open an output file
  //
  OpenFile();

  // Define nTuple column if needed
  //
  DefineNTColumn();
}

void GRRunAction::OpenFile()
{
  if(!fileOpen)
  {
    auto analysisManager = G4AnalysisManager::Instance();
    analysisManager->OpenFile(fileName);
    if(verbose>0) G4cout << "GRRunAction::BeginOfRunAction ### <" << fileName << "> is opened." << G4endl;
    fileOpen = true;
  }
}

void GRRunAction::EndOfRunAction(const G4Run* /*run*/)
{
  // print histogram statistics
  //
  if(!ifCarry) Flush();
}

void GRRunAction::Flush()
{
  auto analysisManager = G4AnalysisManager::Instance();

  analysisManager->Write();
  analysisManager->CloseFile();
  if(verbose>0) G4cout << "GRRunAction::Flush ### <" << fileName << "> is closed." << G4endl;

  fileOpen = false;

  if(IsMaster()) MergeNtuple();
}

void GRRunAction::SetVerbose(G4int val)
{
  verbose = val;
  auto analysisManager = G4AnalysisManager::Instance();
  analysisManager->SetVerboseLevel(verbose);
}

void GRRunAction::ListHistograms()
{
  G4cout << "################## registered histograms/plots" << G4endl;
  G4cout << "id\thistID\thistType\tdetName-X\tpsName-X\tcollID-X\tcopyNo-X\tdetName-Y\tpsName-Y\tcollID-Y\tcopyNo-Y" << G4endl;
  for(auto itr : IDMap)
  {
    G4cout << itr.first << "\t" << itr.second->histID << "\t";
    if(itr.second->histType==1) // 1D histogram
    { G4cout << "1-D hist\t" << itr.second->meshName << "\t" << itr.second->primName << "\t" << itr.second->collID << "\t" << itr.second->idx; }
    else if(itr.second->histType==2) // 1D profile
    { G4cout << "1-D prof\t" << itr.second->meshName << "\t" << itr.second->primName << "\t" << itr.second->collID; }
    G4cout << G4endl;
  }
}

G4bool GRRunAction::Open(G4int id)
{
  auto hItr = IDMap.find(id);
  return (hItr!=IDMap.end());
}

#include "G4SDManager.hh"
using namespace G4Analysis;

G4bool GRRunAction::SetAllPlotting(G4bool val)
{
  G4bool valid = true;
  for(auto hItr : IDMap)
  { 
    valid = SetPlotting(hItr.first,val);
    if(!valid) break;
  }
  return valid;
}

G4bool GRRunAction::SetPlotting(G4int id,G4bool val)
{
  auto hItr = IDMap.find(id);
  if(hItr==IDMap.end()) return false;
  auto ht = hItr->second;
  auto hTyp = ht->histType;
  auto analysisManager = G4AnalysisManager::Instance();
  if(hTyp==1) // 1D histogram
  { analysisManager->SetH1Plotting(ht->histID,val); }
  else if(hTyp==2) // 1D profile
  { analysisManager->SetP1Plotting(ht->histID,val); }
  else
  { return false; }
  return true;
}

// ------------- 1D histogram

G4int GRRunAction::Create1D(G4String& mName,G4String& pName,G4int cn)
{
  G4String collName = mName;
  collName += "/";
  collName += pName;
  auto cID = G4SDManager::GetSDMpointer()->GetCollectionID(collName);
  if(cID<0) return cID;

  G4int id = (cID+id_offset)*id_factor+cn+1;
  auto histoTypeItr = IDMap.find(id);
  if(histoTypeItr!=IDMap.end()) return false;
  if(verbose) G4cout << "GRRunAction::Create1D for <" << collName
                     << ", copyNo=" << cn << "> is registered for hitCollectionID "
                     << cID << G4endl;
  
  auto histTyp = new GRHistoType;
  histTyp->collID = cID;
  histTyp->histType = 1; // 1D histogram
  histTyp->meshName = mName;
  histTyp->primName = pName;
  histTyp->idx = cn;
  IDMap[id] = histTyp;
  return id;
}

G4int GRRunAction::Create1DForPrimary(G4String& mName,G4bool wgt)
{
  G4int cn = wgt ? 1 : 0;

  G4int id = 99999 - cn;
  auto histoTypeItr = IDMap.find(id);
  if(histoTypeItr!=IDMap.end()) return false;
  if(verbose) G4cout << "GRRunAction::Create1D for <" << mName
                     << "(weighted : " << cn << ")> is registered " << G4endl;

  auto histTyp = new GRHistoType;
  histTyp->collID = -999;
  histTyp->histType = 1; // 1D histogram
  histTyp->meshName = "PrimPEnergy";
  histTyp->primName = mName;
  histTyp->biasf = cn;
  IDMap[id] = histTyp;
  return id;
}

#include "G4SDManager.hh"
#include "G4ScoringManager.hh"
#include "G4VScoringMesh.hh"
#include "G4VPrimitivePlotter.hh"
G4int GRRunAction::Create1DForPlotter(G4String& mName,G4String& pName,G4bool /*wgt*/)
{
  using MeshShape = G4VScoringMesh::MeshShape;

  G4String collName = mName;
  collName += "/";
  collName += pName;
  auto cID = G4SDManager::GetSDMpointer()->GetCollectionID(collName);
  if(cID<0) return cID;

  auto sm = G4ScoringManager::GetScoringManagerIfExist();
  assert(sm!=nullptr);
  auto mesh = sm->FindMesh(mName);
  if(mesh==nullptr) 
  { return -2; }
  auto shape = mesh->GetShape();
  if(shape!=MeshShape::realWorldLogVol && shape!=MeshShape::probe)
  { return -3; }
  G4int nBin[3];
  mesh->GetNumberOfSegments(nBin);

  auto prim = mesh->GetPrimitiveScorer(pName);
  if(prim==nullptr)
  { return -3; }
  auto pp = dynamic_cast<G4VPrimitivePlotter*>(prim);
  if(pp==nullptr)
  { return -4; }
  
  G4int id0 = (cID+id_offset)*id_factor+1;
  for(G4int cn=0; cn<nBin[0]; cn++)
  {
    G4int id = id0+cn;
    auto histoTypeItr = IDMap.find(id);
    if(histoTypeItr!=IDMap.end())
    { return -5; }
    if(verbose) G4cout << "GRRunAction::Create1D for <" << collName
                     << ", copyNo=" << cn << "> is registered for hitCollectionID "
                     << cID << G4endl;
  
    auto histTyp = new GRHistoType;
    histTyp->collID = cID;
    histTyp->histType = 1; // 1D histogram
    histTyp->histDup = nBin[0];
    histTyp->meshName = mName;
    histTyp->primName = pName;
    histTyp->idx = cn;
    histTyp->pplotter = pp;
    IDMap[id] = histTyp;
  }
  return id0;
}  

#include "G4UIcommand.hh"
G4bool GRRunAction::Set1D(G4int id0,G4int nBin,G4double valMin,G4double valMax,G4String& unit,
                          G4String& schem, G4bool logVal)
{
  OpenFile();

  auto hIt = IDMap.find(id0);
  if(hIt==IDMap.end()) return false;

  auto analysisManager = G4AnalysisManager::Instance();
  auto dup = (hIt->second)->histDup;
  for(G4int ii=0;ii<dup;ii++)
  {
    G4int id = id0 + ii;
    auto hItr = IDMap.find(id);
    auto ht = hItr->second;
    G4String mNam = ht->primName;
    G4String nam = ht->meshName + "_" + ht->primName;
    if(ht->idx>-1)
    { 
      mNam += "_";
      mNam += G4UIcommand::ConvertToString(ht->idx);
      nam += "_";
      nam += G4UIcommand::ConvertToString(ht->idx);
    }
    G4int hid = -1;
    if(schem=="linear")
    { hid = analysisManager->CreateH1(mNam,nam,nBin,valMin,valMax,unit,"none","linear"); }
    else
    {
      if(logVal)
      { hid = analysisManager->CreateH1(mNam,nam,nBin,valMin,valMax,unit,"log10","linear"); }
      else
      {
        hid = analysisManager->CreateH1(mNam,nam,nBin,valMin,valMax,unit,"none","log");
        analysisManager->SetH1XAxisIsLog(hid,true);
      }
    }

    if(verbose) G4cout << "GRRunAction::Set1D for " << mNam << " / " << nam
                       << " has the histogram ID " << hid << G4endl;
    ht->histID = hid;
    auto pp = ht->pplotter;
    if(pp!=nullptr) pp->Plot(ht->idx,hid);
  }
  return true;
}

G4bool GRRunAction::Set1DTitle(G4int id,G4String& title,G4String& x_axis,G4String&y_axis)
{
  auto hItr = IDMap.find(id);
  if(hItr==IDMap.end()) return false;

  auto analysisManager = G4AnalysisManager::Instance();
  auto hid = hItr->second->histID;
  analysisManager->SetH1Title(hid,title);
  analysisManager->SetH1XAxisTitle(hid,x_axis);
  analysisManager->SetH1YAxisTitle(hid,y_axis);
  return true;
}

G4bool GRRunAction::Set1DYAxisLog(G4int id0,G4bool val)
{
  auto hIt = IDMap.find(id0);
  if(hIt==IDMap.end()) return false;
  auto analysisManager = G4AnalysisManager::Instance();
  auto dup = (hIt->second)->histDup;
  for(G4int ii=0;ii<dup;ii++)
  {
    G4int id = id0 + ii;
    auto hItr = IDMap.find(id);
    analysisManager->SetH1YAxisIsLog(hItr->second->histID,val);
  }
  return true;
}
  
// ------------- 1D profile

G4int GRRunAction::Create1P(G4String& mName,G4String& pName,G4int cn)
{
  G4String collName = mName;
  collName += "/";
  collName += pName;
  auto cID = G4SDManager::GetSDMpointer()->GetCollectionID(collName);
  if(cID<0) return cID;

  G4int id = (cID+2*id_offset)*id_factor;
  auto histoTypeItr = IDMap.find(id);
  if(histoTypeItr!=IDMap.end()) return false;
  if(verbose) G4cout << "GRRunAction::Create1P for <" << collName
                     << "> is registered for hitCollectionID "
                     << cID << G4endl;

  auto histTyp = new GRHistoType;
  histTyp->collID = cID;
  histTyp->histType = 2; // 1D profile
  histTyp->meshName = mName;
  histTyp->primName = pName;
  histTyp->idx = cn;
  IDMap[id] = histTyp;
  return id;
}

G4bool GRRunAction::Set1P(G4int id,G4double valYMin,G4double valYMax,G4String& unit,
        G4String& funcX,G4String& funcY,G4String& schem)
{
  OpenFile();

  if(verbose) G4cout << "GRRunAction::Set1P for id = " << id << G4endl;
  auto hItr = IDMap.find(id);
  if(hItr==IDMap.end()) return false;

  auto ht = hItr->second;
  if(verbose) G4cout << "GRRunAction::Set1P for " << ht->meshName << " / " << ht->primName << G4endl;
  auto analysisManager = G4AnalysisManager::Instance();
  auto nBin = ht->idx;
  G4double valMin = -0.5;
  G4double valMax = G4double(nBin) - 0.5;
  G4String nam = ht->meshName + "_" + ht->primName;
  auto hid = analysisManager->CreateP1(nam,ht->primName,nBin,
              valMin,valMax,valYMin,valYMax,"none",unit,funcX,funcY,schem);

  if(verbose) G4cout << "GRRunAction::Set1P for " << ht->meshName << " / " << ht->primName
                     << " has the histogram ID " << hid << G4endl;
  ht->histID = hid;
  return true;
}

G4bool GRRunAction::Set1PTitle(G4int id,G4String& title,G4String& x_axis,G4String&y_axis)
{
  auto hItr = IDMap.find(id);
  if(hItr==IDMap.end()) return false;

  auto analysisManager = G4AnalysisManager::Instance();
  auto hid = hItr->second->histID;
  analysisManager->SetP1Title(hid,title);
  analysisManager->SetP1XAxisTitle(hid,x_axis);
  analysisManager->SetP1YAxisTitle(hid,y_axis);
  return true;
}

// ------------- Ntuple

G4int GRRunAction::NtupleColumn(G4String& mName,G4String& pName,G4String& unit,G4int cn)
{
  G4String collName = mName;
  collName += "/";
  collName += pName;
  auto cID = G4SDManager::GetSDMpointer()->GetCollectionID(collName);
  if(cID<0) return cID;

  G4int id = NTMap.size();
  if(verbose) G4cout << "GRRunAction::NtupleColumn : <" << collName
                     << ", copyNo=" << cn << "> is registered for nTuple column "
                     << id << G4endl;

  auto histTyp = new GRHistoType;
  histTyp->collID = cID;
  histTyp->meshName = mName;
  histTyp->primName = pName;
  histTyp->meshName2 = unit;
  if(unit!="none")
  { histTyp->fuct = 1./(G4UnitDefinition::GetValueOf(unit)); }
  histTyp->idx = cn;
  NTMap[id] = histTyp;
  return id;
}

#include "G4UIcommand.hh"

void GRRunAction::DefineNTColumn()
{
  if(NTMap.size()==0) return;

  auto analysisManager = G4AnalysisManager::Instance();
  analysisManager->CreateNtuple("GRimNtuple","Scores for each event");
  
  for(auto itr : NTMap)
  {
    G4String colNam = itr.second->meshName;
    colNam += "_";
    colNam += itr.second->primName;
    if(itr.second->idx != -1)
    { colNam += "_"; colNam += G4UIcommand::ConvertToString(itr.second->idx); }
    if(itr.second->meshName2 != "none")
    { colNam += "["; colNam += itr.second->meshName2; colNam += "]"; }
    analysisManager->CreateNtupleDColumn(colNam);
  }

  analysisManager->FinishNtuple();
}

#include <fstream>
#include "G4Threading.hh"
#include "G4UImanager.hh"

void GRRunAction::MergeNtuple()
{
  if(NTMap.size()==0) return;
  if(!(G4Threading::IsMultithreadedApplication())) return;

  auto analysisManager = G4AnalysisManager::Instance();

  // This MergeNtuple() method is valid only for CSV file format
  if(analysisManager->GetType()!="Csv") return;

  std::fstream target;
  G4String targetFN = "GRimOut_nt_GRimNtuple_total.csv";
  target.open(targetFN,std::ios::out);

  enum { BUFSIZE = 4096 };
  char* line = new char[BUFSIZE];

  G4String titleFN = "GRimOut_nt_GRimNtuple.csv";
  std::ifstream title;
  title.open(titleFN,std::ios::in);
  while(1)
  {
    title.getline(line,BUFSIZE);
    if(title.eof()) break;
    G4cout << line << G4endl;
    target << line << G4endl;
  }
  title.close();

  auto nWorker = G4Threading::GetNumberOfRunningWorkerThreads();
  G4String sourceFNBase = "GRimOut_nt_GRimNtuple_t";
  for(G4int i = 0; i < nWorker; i++)
  {
    G4String sourceFN = sourceFNBase;
    sourceFN += G4UIcommand::ConvertToString(i);
    sourceFN += ".csv";
    std::ifstream source;
    source.open(sourceFN,std::ios::in);
    if(!source)
    {
      G4ExceptionDescription ed; ed << "Source file <" << sourceFN << "> is not found.";
      G4Exception("GRRunAction::MergeNtuple()","GRim12345",FatalException,ed);
    }
    while(1)
    {
      source.getline(line,BUFSIZE);
      if(line[0]=='#') continue;
      if(source.eof()) break;
      target << line << G4endl;
    }
    source.close();
    G4String scmd = "rm -f ";
    scmd += sourceFN;
    auto rc = system(scmd);
    if(rc<0)
    {
      G4ExceptionDescription ed; 
      ed << "File <" << sourceFN << "> could not be deleted, thought it is merged.";
      G4Exception("GRRunAction::MergeNtuple()","GRim12345",JustWarning,ed);
    }
  }

  target.close();

  G4String cmd = "mv ";
  cmd += targetFN;
  cmd += " ";
  cmd += titleFN;
  auto rcd = system(cmd);
  if(rcd<0)
  {
    G4ExceptionDescription ed; 
    ed << "File <" << targetFN << "> could not be renamed.";
    G4Exception("GRRunAction::MergeNtuple()","GRim12345",JustWarning,ed);
  }
}


