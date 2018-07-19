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
// $Id: G4VAnalysisManager.cc 70604 2013-06-03 11:27:06Z ihrivnac $

// Author: Ivana Hrivnacova, 09/07/2013  (ivana@ipno.in2p3.fr)

#include "G4VAnalysisManager.hh"
#include "G4AnalysisMessenger.hh"
#include "G4AnalysisUtilities.hh"
#include "G4HnManager.hh"
#include "G4VH1Manager.hh"
#include "G4VH2Manager.hh"
#include "G4VH3Manager.hh"
#include "G4VP1Manager.hh"
#include "G4VP2Manager.hh"
#include "G4VNtupleManager.hh"
#include "G4VFileManager.hh"

#include <iostream>

using namespace G4Analysis;

//_____________________________________________________________________________
G4VAnalysisManager::G4VAnalysisManager(const G4String& type, G4bool isMaster)
 : fState(type, isMaster),
   fVFileManager(nullptr),
   fMessenger(G4Analysis::make_unique<G4AnalysisMessenger>(this)),
   fH1HnManager(nullptr),
   fH2HnManager(nullptr),
   fH3HnManager(nullptr),
   fP1HnManager(nullptr),
   fP2HnManager(nullptr),
   fVH1Manager(nullptr),
   fVH2Manager(nullptr),
   fVH3Manager(nullptr),
   fVP1Manager(nullptr),
   fVP2Manager(nullptr),
   fVNtupleManager(nullptr)
{
  //fMessenger = G4Analysis::make_unique<G4AnalysisMessenger>(this);
}

//_____________________________________________________________________________
G4VAnalysisManager::~G4VAnalysisManager()
{
  delete fVNtupleManager;
}

// 
// protected methods
//

//_____________________________________________________________________________
void G4VAnalysisManager::SetH1Manager(G4VH1Manager* h1Manager)
{
  fVH1Manager.reset(h1Manager);
  fH1HnManager = h1Manager->GetHnManager();
  fMessenger->SetH1HnManager(*fH1HnManager);
} 

//_____________________________________________________________________________
void G4VAnalysisManager::SetH2Manager(G4VH2Manager* h2Manager)
{
  fVH2Manager.reset(h2Manager);
  fH2HnManager = h2Manager->GetHnManager();
  fMessenger->SetH2HnManager(*fH2HnManager);
}  

//_____________________________________________________________________________
void G4VAnalysisManager::SetH3Manager(G4VH3Manager* h3Manager)
{
  fVH3Manager.reset(h3Manager);
  fH3HnManager = h3Manager->GetHnManager();
  fMessenger->SetH3HnManager(*fH3HnManager);
}  

//_____________________________________________________________________________
void G4VAnalysisManager::SetP1Manager(G4VP1Manager* p1Manager)
{
  fVP1Manager.reset(p1Manager);
  fP1HnManager = p1Manager->GetHnManager();
  fMessenger->SetP1HnManager(*fP1HnManager);
} 

//_____________________________________________________________________________
void G4VAnalysisManager::SetP2Manager(G4VP2Manager* p2Manager)
{
  fVP2Manager.reset(p2Manager);
  fP2HnManager = p2Manager->GetHnManager();
  fMessenger->SetP2HnManager(*fP2HnManager);
} 

//_____________________________________________________________________________
void G4VAnalysisManager::SetNtupleManager(G4VNtupleManager* ntupleManager)
{
  // fVNtupleManager.reset(ntupleManager);
  fVNtupleManager = ntupleManager;
}  

//_____________________________________________________________________________
void G4VAnalysisManager::SetFileManager(std::shared_ptr<G4VFileManager> fileManager)
{
  fVFileManager = fileManager;
}  

//_____________________________________________________________________________
G4bool G4VAnalysisManager::WriteAscii(const G4String& fileName)
{
  G4bool finalResult = true;

  // Replace or add file extension .ascii
  G4String name(fileName);
  if ( name.find(".") != std::string::npos ) { 
    name.erase(name.find("."), name.length()); 
  }
  name.append(".ascii");

#ifdef G4VERBOSE
  if ( fState.GetVerboseL3() ) 
    fState.GetVerboseL3()->Message("write ASCII", "file", name);
#endif
     
  std::ofstream output(name, std::ios::out);
  if ( ! output ) {
    G4ExceptionDescription description;
    description 
      << "Cannot open file. File name is not defined.";
    G4Exception("G4VAnalysisManager::WriteAscii()",
                "Analysis_W001", JustWarning, description);
    return false;
  }
  output.setf( std::ios::scientific, std::ios::floatfield );

  G4bool result = fVH1Manager->WriteOnAscii(output);
  finalResult = finalResult && result;
  
  result = fVH2Manager->WriteOnAscii(output);
  finalResult = finalResult && result;  

  result = fVH3Manager->WriteOnAscii(output);
  finalResult = finalResult && result;  

  //result = fVP1Manager->WriteOnAscii(output);
  //finalResult = finalResult && result;  

  //result = fVP2Manager->WriteOnAscii(output);
  //finalResult = finalResult && result;  

#ifdef G4VERBOSE
    if ( fState.GetVerboseL1() ) 
      fState.GetVerboseL1()->Message("write ASCII", "file",  name, result);
#endif
  
  return finalResult;
}     

// 
// public methods
//

//_____________________________________________________________________________
G4bool G4VAnalysisManager::OpenFile(const G4String& fileName)
{
  if ( fileName != "" ) {
    return OpenFileImpl(fileName);
  }
  else {  
    if ( fVFileManager->GetFileName() == "" ) {
      G4ExceptionDescription description;
      description 
        << "Cannot open file. File name is not defined.";
      G4Exception("G4VFileManager::OpenFile()",
                  "Analysis_W001", JustWarning, description);
      return false;
    }           
    return OpenFileImpl(fVFileManager->GetFileName());
  }  
} 

//_____________________________________________________________________________
G4bool G4VAnalysisManager::Write()
{
  G4bool finalResult = true;

  G4bool result = WriteImpl();
  finalResult = finalResult && result;
 
  if ( IsPlotting() ) {
    result = PlotImpl();
    finalResult = finalResult && result;
  }

  return finalResult;
}  

//_____________________________________________________________________________
G4bool G4VAnalysisManager::CloseFile()
{
  return CloseFileImpl();
}  

//_____________________________________________________________________________
G4bool G4VAnalysisManager::Merge(tools::histo::hmpi* hmpi)
{
  return MergeImpl(hmpi);
}  

//_____________________________________________________________________________
G4bool G4VAnalysisManager::Plot()
{
  return PlotImpl();
}  

//_____________________________________________________________________________
G4bool G4VAnalysisManager::IsOpenFile() const
{
  return IsOpenFileImpl();
}  

//_____________________________________________________________________________
G4bool G4VAnalysisManager::SetFileName(const G4String& fileName)
{ 
  return fVFileManager->SetFileName(fileName); 
}

//_____________________________________________________________________________
G4bool G4VAnalysisManager::SetHistoDirectoryName(const G4String& dirName)
{ 
  return fVFileManager->SetHistoDirectoryName(dirName); 
}

//_____________________________________________________________________________
G4bool G4VAnalysisManager::SetNtupleDirectoryName(const G4String& dirName)
{ 
  return fVFileManager->SetNtupleDirectoryName(dirName); 
}

//_____________________________________________________________________________
void G4VAnalysisManager::SetCompressionLevel(G4int level)
{
  fState.SetCompressionLevel(level);
}

//_____________________________________________________________________________
G4String G4VAnalysisManager::GetFileName() const 
{  
  return fVFileManager->GetFileName(); 
}

//_____________________________________________________________________________
G4String G4VAnalysisManager::GetHistoDirectoryName() const 
{  
  return fVFileManager->GetHistoDirectoryName(); 
}
 
//_____________________________________________________________________________
G4String G4VAnalysisManager::GetNtupleDirectoryName() const
{
  return fVFileManager->GetNtupleDirectoryName(); 
}

//_____________________________________________________________________________
G4int G4VAnalysisManager::GetCompressionLevel() const
{
  return fState.GetCompressionLevel();
}

//_____________________________________________________________________________
G4int G4VAnalysisManager::CreateH1(const G4String& name,  const G4String& title,
                               G4int nbins, G4double xmin, G4double xmax,
                               const G4String& unitName, const G4String& fcnName,
                               const G4String& binSchemeName)
{
  if ( ! CheckName(name, "H1") ) return kInvalidId;
  if ( ! CheckNbins(nbins) ) return kInvalidId;
  if ( ! CheckMinMax(xmin, xmax, fcnName, binSchemeName) ) return kInvalidId;

  return fVH1Manager->CreateH1(name, title, nbins, xmin, xmax, 
                               unitName, fcnName, binSchemeName);
}                                         

//_____________________________________________________________________________
G4int G4VAnalysisManager::CreateH1(const G4String& name,  const G4String& title,
                               const std::vector<G4double>& edges,
                               const G4String& unitName, const G4String& fcnName)
{
  if ( ! CheckName(name, "H1") ) return kInvalidId;
  if ( ! CheckEdges(edges) ) return kInvalidId;

  return fVH1Manager->CreateH1(name, title, edges, unitName, fcnName);
}                                         

//_____________________________________________________________________________
G4int G4VAnalysisManager::CreateH2(const G4String& name,  const G4String& title,
                               G4int nxbins, G4double xmin, G4double xmax,
                               G4int nybins, G4double ymin, G4double ymax,
                               const G4String& xunitName, const G4String& yunitName,
                               const G4String& xfcnName, const G4String& yfcnName,
                               const G4String& xbinSchemeName, 
                               const G4String& ybinSchemeName)
                               
{
  if ( ! CheckName(name, "H2") ) return kInvalidId;
  
  if ( ! CheckNbins(nxbins) ) return kInvalidId;
  if ( ! CheckMinMax(xmin, xmax, xfcnName, xbinSchemeName) ) return kInvalidId;

  if ( ! CheckNbins(nybins) ) return kInvalidId;
  if ( ! CheckMinMax(ymin, ymax, yfcnName, ybinSchemeName) ) return kInvalidId;

  return fVH2Manager->CreateH2(name, title, 
                               nxbins, xmin, xmax, nybins, ymin, ymax, 
                               xunitName, yunitName, xfcnName, yfcnName, 
                               xbinSchemeName, ybinSchemeName);
}                                         

//_____________________________________________________________________________
G4int G4VAnalysisManager::CreateH2(const G4String& name,  const G4String& title,
                               const std::vector<G4double>& xedges,
                               const std::vector<G4double>& yedges,
                               const G4String& xunitName, const G4String& yunitName,
                               const G4String& xfcnName, const G4String& yfcnName)
                               
{
  if ( ! CheckName(name, "H2") ) return kInvalidId;
  
  if ( ! CheckEdges(xedges) ) return kInvalidId;
  if ( ! CheckEdges(yedges) ) return kInvalidId;

  return fVH2Manager->CreateH2(name, title, 
                               xedges, yedges,
                               xunitName, yunitName, xfcnName, yfcnName);
}                                         

//_____________________________________________________________________________
G4int G4VAnalysisManager::CreateH3(const G4String& name,  const G4String& title,
                               G4int nxbins, G4double xmin, G4double xmax,
                               G4int nybins, G4double ymin, G4double ymax,
                               G4int nzbins, G4double zmin, G4double zmax,
                               const G4String& xunitName, const G4String& yunitName,
                               const G4String& zunitName,
                               const G4String& xfcnName, const G4String& yfcnName, 
                               const G4String& zfcnName,
                               const G4String& xbinSchemeName, 
                               const G4String& ybinSchemeName,
                               const G4String& zbinSchemeName)
                               
{
  if ( ! CheckName(name, "H3") ) return kInvalidId;
  
  if ( ! CheckNbins(nxbins) ) return kInvalidId;
  if ( ! CheckMinMax(xmin, xmax, xfcnName, xbinSchemeName) ) return kInvalidId;

  if ( ! CheckNbins(nybins) ) return kInvalidId;
  if ( ! CheckMinMax(ymin, ymax, yfcnName, ybinSchemeName) ) return kInvalidId;

  if ( ! CheckNbins(nzbins) ) return kInvalidId;
  if ( ! CheckMinMax(zmin, zmax, zfcnName, zbinSchemeName) ) return kInvalidId;

  return fVH3Manager->CreateH3(name, title, 
                               nxbins, xmin, xmax, nybins, ymin, ymax, 
                               nzbins, zmin, zmax,
                               xunitName, yunitName, zunitName, 
                               xfcnName, yfcnName, zfcnName,
                               xbinSchemeName, ybinSchemeName, zbinSchemeName);
}                                         

//_____________________________________________________________________________
G4int G4VAnalysisManager::CreateH3(const G4String& name,  const G4String& title,
                               const std::vector<G4double>& xedges,
                               const std::vector<G4double>& yedges,
                               const std::vector<G4double>& zedges,
                               const G4String& xunitName, const G4String& yunitName,
                               const G4String& zunitName,
                               const G4String& xfcnName, const G4String& yfcnName, 
                               const G4String& zfcnName)
                               
{
  if ( ! CheckName(name, "H3") ) return kInvalidId;
  
  if ( ! CheckEdges(xedges) ) return kInvalidId;
  if ( ! CheckEdges(yedges) ) return kInvalidId;
  if ( ! CheckEdges(zedges) ) return kInvalidId;

  return fVH3Manager->CreateH3(name, title, 
                               xedges, yedges, zedges, 
                               xunitName, yunitName, zunitName,
                               xfcnName, yfcnName, zfcnName);
}                                         

//_____________________________________________________________________________
G4bool G4VAnalysisManager::SetH1(G4int id,
                                G4int nbins, G4double xmin, G4double xmax,
                                const G4String& unitName, const G4String& fcnName,
                                const G4String& binSchemeName)
{                                
  if ( ! CheckNbins(nbins) ) return kInvalidId;
  if ( ! CheckMinMax(xmin, xmax, fcnName, binSchemeName) ) return kInvalidId;

  return fVH1Manager->SetH1(id, nbins, xmin, xmax, unitName, fcnName, binSchemeName); 
}
  
//_____________________________________________________________________________
G4bool G4VAnalysisManager::SetH1(G4int id,
                                const std::vector<G4double>& edges,
                                const G4String& unitName, const G4String& fcnName)
{                                
  if ( ! CheckEdges(edges) ) return kInvalidId;

  return fVH1Manager->SetH1(id, edges, unitName, fcnName); 
}
  
//_____________________________________________________________________________
G4bool G4VAnalysisManager::SetH2(G4int id,
                                G4int nxbins, G4double xmin, G4double xmax, 
                                G4int nybins, G4double ymin, G4double ymax,
                                const G4String& xunitName, const G4String& yunitName,
                                const G4String& xfcnName, const G4String& yfcnName,
                                const G4String& xbinSchemeName, 
                                const G4String& ybinSchemeName)
{                                
  if ( ! CheckNbins(nxbins) ) return kInvalidId;
  if ( ! CheckMinMax(xmin, xmax, xfcnName, xbinSchemeName) ) return kInvalidId;

  if ( ! CheckNbins(nybins) ) return kInvalidId;
  if ( ! CheckMinMax(ymin, ymax, yfcnName, ybinSchemeName) ) return kInvalidId;

  return fVH2Manager->SetH2(id, nxbins, xmin, xmax, nybins, ymin, ymax, 
                            xunitName, yunitName, xfcnName, yfcnName,
                            xbinSchemeName, ybinSchemeName);
}
                                  
//_____________________________________________________________________________
G4bool G4VAnalysisManager::SetH2(G4int id,
                                const std::vector<G4double>& xedges,
                                const std::vector<G4double>& yedges,
                                const G4String& xunitName, const G4String& yunitName,
                                const G4String& xfcnName, const G4String& yfcnName)
{                                
  if ( ! CheckEdges(xedges) ) return kInvalidId;
  if ( ! CheckEdges(yedges) ) return kInvalidId;

  return fVH2Manager->SetH2(id, xedges, yedges, 
                            xunitName, yunitName, xfcnName, yfcnName);
}
                                  
//_____________________________________________________________________________
G4bool G4VAnalysisManager::SetH3(G4int id,
                                G4int nxbins, G4double xmin, G4double xmax, 
                                G4int nybins, G4double ymin, G4double ymax,
                                G4int nzbins, G4double zmin, G4double zmax,
                                const G4String& xunitName, const G4String& yunitName,
                                const G4String& zunitName,
                                const G4String& xfcnName, const G4String& yfcnName, 
                                const G4String& zfcnName,
                                const G4String& xbinSchemeName, 
                                const G4String& ybinSchemeName,
                                const G4String& zbinSchemeName)
{                                
  if ( ! CheckNbins(nxbins) ) return kInvalidId;
  if ( ! CheckMinMax(xmin, xmax, xfcnName, xbinSchemeName) ) return kInvalidId;

  if ( ! CheckNbins(nybins) ) return kInvalidId;
  if ( ! CheckMinMax(ymin, ymax, yfcnName, ybinSchemeName) ) return kInvalidId;

  if ( ! CheckNbins(nzbins) ) return kInvalidId;
  if ( ! CheckMinMax(zmin, zmax, zfcnName, zbinSchemeName) ) return kInvalidId;

  return fVH3Manager->SetH3(id, 
                            nxbins, xmin, xmax, nybins, ymin, ymax, 
                            nzbins, zmin, zmax,
                            xunitName, yunitName, zunitName, 
                            xfcnName, yfcnName, zfcnName,
                            xbinSchemeName, ybinSchemeName, zbinSchemeName);
}
                                  
//_____________________________________________________________________________
G4bool G4VAnalysisManager::SetH3(G4int id,
                                const std::vector<G4double>& xedges,
                                const std::vector<G4double>& yedges,
                                const std::vector<G4double>& zedges,
                                const G4String& xunitName, const G4String& yunitName,
                                const G4String& zunitName,
                                const G4String& xfcnName, const G4String& yfcnName, 
                                const G4String& zfcnName)
{                                
  if ( ! CheckEdges(xedges) ) return kInvalidId;
  if ( ! CheckEdges(yedges) ) return kInvalidId;
  if ( ! CheckEdges(zedges) ) return kInvalidId;

  return fVH3Manager->SetH3(id, xedges, yedges, zedges,
                            xunitName, yunitName, zunitName, 
                            xfcnName, yfcnName, zfcnName);
}
                                  
//_____________________________________________________________________________
G4bool G4VAnalysisManager::ScaleH1(G4int id, G4double factor)
{
  return fVH1Manager->ScaleH1(id, factor);
}  

//_____________________________________________________________________________
G4bool G4VAnalysisManager::ScaleH2(G4int id, G4double factor)
{
  return fVH2Manager->ScaleH2(id, factor);
}  
                           
//_____________________________________________________________________________
G4bool G4VAnalysisManager::ScaleH3(G4int id, G4double factor)
{
  return fVH3Manager->ScaleH3(id, factor);
}  
                           
//_____________________________________________________________________________
G4int G4VAnalysisManager::CreateP1(const G4String& name,  const G4String& title,
                               G4int nbins, G4double xmin, G4double xmax,
                               G4double ymin, G4double ymax,
                               const G4String& xunitName, const G4String& yunitName,
                               const G4String& xfcnName, const G4String& yfcnName,
                               const G4String& xbinSchemeName)
{
  if ( ! CheckName(name, "P1") ) return kInvalidId;
  if ( ! CheckNbins(nbins) ) return kInvalidId;
  if ( ! CheckMinMax(xmin, xmax, xfcnName, xbinSchemeName) ) return kInvalidId;
  if ( ymin != 0. || ymax != 0. ) {
    // Do not check  default values
    if ( ! CheckMinMax(ymin, ymax) ) return kInvalidId;
  }

  return fVP1Manager->CreateP1(name, title, nbins, xmin, xmax, ymin, ymax,
                               xunitName, yunitName, xfcnName, yfcnName, 
                               xbinSchemeName);
}                                         

//_____________________________________________________________________________
G4int G4VAnalysisManager::CreateP1(const G4String& name,  const G4String& title,
                               const std::vector<G4double>& edges,
                               G4double ymin, G4double ymax,
                               const G4String& xunitName, const G4String& yunitName,
                               const G4String& xfcnName, const G4String& yfcnName)
{
  if ( ! CheckName(name, "P1") ) return kInvalidId;
  if ( ! CheckEdges(edges) ) return kInvalidId;
  if ( ymin != 0. || ymax != 0. ) {
    // Do not check  default values
    if ( ! CheckMinMax(ymin, ymax) ) return kInvalidId;
  }

  return fVP1Manager->CreateP1(name, title, edges, ymin, ymax, 
                               xunitName, yunitName, xfcnName, yfcnName);
}                                         

//_____________________________________________________________________________
G4int G4VAnalysisManager::CreateP2(const G4String& name, const G4String& title,
                              G4int nxbins, G4double xmin, G4double xmax,
                              G4int nybins, G4double ymin, G4double ymax, 
                              G4double zmin, G4double zmax,
                              const G4String& xunitName, const G4String& yunitName,
                              const G4String& zunitName,
                              const G4String& xfcnName, const G4String& yfcnName,
                              const G4String& zfcnName,
                              const G4String& xbinSchemeName, 
                              const G4String& ybinSchemeName)
{
  if ( ! CheckName(name, "P2") ) return kInvalidId;
  if ( ! CheckNbins(nxbins) ) return kInvalidId;
  if ( ! CheckMinMax(xmin, xmax, xfcnName, xbinSchemeName) ) return kInvalidId;
  if ( ! CheckMinMax(ymin, ymax, yfcnName, xbinSchemeName) ) return kInvalidId;
  if ( zmin != 0. || zmax != 0. ) {
    // Do not check  default values
    if ( ! CheckMinMax(zmin, zmax) ) return kInvalidId;
  }

  return fVP2Manager->CreateP2(name, title, 
                               nxbins, xmin, xmax, nybins, ymin, ymax,
                               zmin, zmax,
                               xunitName, yunitName, zunitName,
                               xfcnName, yfcnName, zfcnName,
                               xbinSchemeName, ybinSchemeName);
}                               

//_____________________________________________________________________________
G4int G4VAnalysisManager::CreateP2(const G4String& name, const G4String& title,
                              const std::vector<G4double>& xedges,
                              const std::vector<G4double>& yedges,
                              G4double zmin, G4double zmax,
                              const G4String& xunitName, const G4String& yunitName,
                              const G4String& zunitName,
                              const G4String& xfcnName, const G4String& yfcnName,
                              const G4String& zfcnName)
{
  if ( ! CheckName(name, "P2") ) return kInvalidId;
  if ( ! CheckEdges(xedges) ) return kInvalidId;
  if ( ! CheckEdges(yedges) ) return kInvalidId;
  if ( zmin != 0. || zmax != 0. ) {
    // Do not check  default values
    if ( ! CheckMinMax(zmin, zmax) ) return kInvalidId;
  }

  return fVP2Manager->CreateP2(name, title, xedges, yedges, zmin, zmax, 
                               xunitName, yunitName, zunitName,
                               xfcnName, yfcnName, zfcnName);
}                                         

//_____________________________________________________________________________
G4bool G4VAnalysisManager::SetP1(G4int id,
                                G4int nbins, G4double xmin, G4double xmax,
                                G4double ymin, G4double ymax,
                                const G4String& xunitName, const G4String& yunitName,
                                const G4String& xfcnName, const G4String& yfcnName,
                                const G4String& xbinSchemeName)
{                                
  if ( ! CheckNbins(nbins) ) return kInvalidId;
  if ( ! CheckMinMax(xmin, xmax, xfcnName, xbinSchemeName) ) return kInvalidId;
  if ( ymin != 0. || ymax != 0. ) {
    // Do not check  default values
    if ( ! CheckMinMax(ymin, ymax) ) return kInvalidId;
  }

  return fVP1Manager->SetP1(id, nbins, xmin, xmax, ymin, ymax, 
                            xunitName, yunitName, xfcnName, yfcnName, 
                            xbinSchemeName);
}
  
//_____________________________________________________________________________
G4bool G4VAnalysisManager::SetP1(G4int id,
                                const std::vector<G4double>& edges,
                                G4double ymin, G4double ymax,
                                const G4String& xunitName, const G4String& yunitName,
                                const G4String& xfcnName, const G4String& yfcnName)
{                                
  if ( ! CheckEdges(edges) ) return kInvalidId;
  if ( ymin != 0. || ymax != 0. ) {
    // Do not check  default values
    if ( ! CheckMinMax(ymin, ymax) ) return kInvalidId;
  }

  return fVP1Manager->SetP1(id, edges, ymin, ymax, 
                            xunitName, yunitName, xfcnName, yfcnName); 
}
  
//_____________________________________________________________________________
G4bool G4VAnalysisManager::SetP2(G4int id,
                              G4int nxbins, G4double xmin, G4double xmax, 
                              G4int nybins, G4double ymin, G4double ymax, 
                              G4double zmin, G4double zmax,
                              const G4String& xunitName, const G4String& yunitName,
                              const G4String& zunitName,
                              const G4String& xfcnName, const G4String& yfcnName,
                              const G4String& zfcnName,
                              const G4String& xbinSchemeName, 
                              const G4String& ybinSchemeName)
{
  if ( ! CheckNbins(nxbins) ) return kInvalidId;
  if ( ! CheckNbins(nybins) ) return kInvalidId;
  if ( ! CheckMinMax(xmin, xmax, xfcnName, xbinSchemeName) ) return kInvalidId;
  if ( ! CheckMinMax(ymin, ymax, yfcnName, ybinSchemeName) ) return kInvalidId;
  if ( zmin != 0. || zmax != 0. ) {
    // Do not check  default values
    if ( ! CheckMinMax(zmin, zmax) ) return kInvalidId;
  }

  return fVP2Manager->SetP2(id, nxbins, xmin, xmax, nybins, ymin, ymax, 
                            zmin, zmax,
                            xunitName, yunitName, zunitName,
                            xfcnName, yfcnName, zfcnName,
                            xbinSchemeName, ybinSchemeName);
}
  
//_____________________________________________________________________________
G4bool G4VAnalysisManager::SetP2(G4int id,
                              const std::vector<G4double>& xedges,
                              const std::vector<G4double>& yedges,
                              G4double zmin, G4double zmax,
                              const G4String& xunitName, 
                              const G4String& yunitName,
                              const G4String& zunitName,
                              const G4String& xfcnName, 
                              const G4String& yfcnName,
                              const G4String& zfcnName)
{
  if ( ! CheckEdges(xedges) ) return kInvalidId;
  if ( ! CheckEdges(yedges) ) return kInvalidId;
  if ( zmin != 0. || zmax != 0. ) {
    // Do not check  default values
    if ( ! CheckMinMax(zmin, zmax) ) return kInvalidId;
  }

  return fVP2Manager->SetP2(id, xedges, yedges, zmin, zmax, 
                            xunitName, yunitName, zunitName,
                            xfcnName, yfcnName, zfcnName);
}

//_____________________________________________________________________________
G4bool G4VAnalysisManager::ScaleP1(G4int id, G4double factor)
{
  return fVP1Manager->ScaleP1(id, factor);
}  

//_____________________________________________________________________________
G4bool G4VAnalysisManager::ScaleP2(G4int id, G4double factor)
{
  return fVP2Manager->ScaleP2(id, factor);
}  

//_____________________________________________________________________________
G4int G4VAnalysisManager::CreateNtuple(const G4String& name, 
                                          const G4String& title)
{
  if ( ! CheckName(name, "Ntuple") ) return kInvalidId;

  return fVNtupleManager->CreateNtuple(name, title);
}                                         

//_____________________________________________________________________________
G4int G4VAnalysisManager::CreateNtupleIColumn(const G4String& name)
{
  if ( ! CheckName(name, "NtupleIColumn") ) return kInvalidId;

  return fVNtupleManager->CreateNtupleIColumn(name, 0);
}  

//_____________________________________________________________________________
G4int G4VAnalysisManager::CreateNtupleFColumn(const G4String& name)
{
  if ( ! CheckName(name, "NtupleFColumn") ) return kInvalidId;

  return fVNtupleManager->CreateNtupleFColumn(name, 0);
}  

//_____________________________________________________________________________
G4int G4VAnalysisManager::CreateNtupleDColumn(const G4String& name)
{
  if ( ! CheckName(name, "NtupleDColumn") ) return kInvalidId;

  return fVNtupleManager->CreateNtupleDColumn(name, 0);
}  

//_____________________________________________________________________________
G4int G4VAnalysisManager::CreateNtupleSColumn(const G4String& name)
{
  if ( ! CheckName(name, "NtupleSColumn") ) return kInvalidId;

  return fVNtupleManager->CreateNtupleSColumn(name);
}  

//_____________________________________________________________________________
G4int G4VAnalysisManager::CreateNtupleIColumn(const G4String& name, 
                                              std::vector<int>& vector)
{
  if ( ! CheckName(name, "NtupleIColumn") ) return kInvalidId;

  return fVNtupleManager->CreateNtupleIColumn(name, &vector);
}  

//_____________________________________________________________________________
G4int G4VAnalysisManager::CreateNtupleFColumn(const G4String& name, 
                                              std::vector<float>& vector)
{
  if ( ! CheckName(name, "NtupleFColumn") ) return kInvalidId;

  return fVNtupleManager->CreateNtupleFColumn(name, &vector);
}  

//_____________________________________________________________________________
G4int G4VAnalysisManager::CreateNtupleDColumn(const G4String& name, 
                                              std::vector<double>& vector)
{
  if ( ! CheckName(name, "NtupleDColumn") ) return kInvalidId;

  return fVNtupleManager->CreateNtupleDColumn(name, &vector);
}  

//_____________________________________________________________________________
void G4VAnalysisManager::FinishNtuple()
{ 
  return fVNtupleManager->FinishNtuple();
}
   
//_____________________________________________________________________________
G4int G4VAnalysisManager::CreateNtupleIColumn(G4int ntupleId, 
                                              const G4String& name)
{
  if ( ! CheckName(name, "NtupleIColumn") ) return kInvalidId;

  return fVNtupleManager->CreateNtupleIColumn(ntupleId, name, 0);
}                                         

//_____________________________________________________________________________
G4int G4VAnalysisManager::CreateNtupleFColumn(G4int ntupleId, 
                                              const G4String& name)
{
  if ( ! CheckName(name, "NtupleFColumn") ) return kInvalidId;

  return fVNtupleManager->CreateNtupleFColumn(ntupleId, name, 0);
}                                         


//_____________________________________________________________________________
G4int G4VAnalysisManager::CreateNtupleDColumn(G4int ntupleId, 
                                              const G4String& name)   
{
  if ( ! CheckName(name, "NtupleDColumn") ) return kInvalidId;

  return fVNtupleManager->CreateNtupleDColumn(ntupleId, name, 0);
}                                         

//_____________________________________________________________________________
G4int G4VAnalysisManager::CreateNtupleSColumn(G4int ntupleId, 
                                              const G4String& name)   
{
  if ( ! CheckName(name, "NtupleSColumn") ) return kInvalidId;

  return fVNtupleManager->CreateNtupleSColumn(ntupleId, name);
}                                         

//_____________________________________________________________________________
G4int G4VAnalysisManager::CreateNtupleIColumn(G4int ntupleId, 
                                              const G4String& name, 
                                              std::vector<int>& vector)
{
  if ( ! CheckName(name, "NtupleIColumn") ) return kInvalidId;

  return fVNtupleManager->CreateNtupleIColumn(ntupleId, name, &vector);
}  

//_____________________________________________________________________________
G4int G4VAnalysisManager::CreateNtupleFColumn(G4int ntupleId, 
                                              const G4String& name, 
                                              std::vector<float>& vector)
{
  if ( ! CheckName(name, "NtupleFColumn") ) return kInvalidId;

  return fVNtupleManager->CreateNtupleFColumn(ntupleId, name, &vector);
}  

//_____________________________________________________________________________
G4int G4VAnalysisManager::CreateNtupleDColumn(G4int ntupleId, 
                                              const G4String& name, 
                                              std::vector<double>& vector)
{
  if ( ! CheckName(name, "NtupleDColumn") ) return kInvalidId;

  return fVNtupleManager->CreateNtupleDColumn(ntupleId, name, &vector);
}  

//_____________________________________________________________________________
void G4VAnalysisManager::FinishNtuple(G4int ntupleId)
{ 
  return fVNtupleManager->FinishNtuple(ntupleId);
}
   
//_____________________________________________________________________________
G4bool G4VAnalysisManager::SetFirstHistoId(G4int firstId) 
{
  G4bool finalResult = true;

  G4bool result = SetFirstH1Id(firstId);
  finalResult = finalResult && result;
  
  result = SetFirstH2Id(firstId);
  finalResult = finalResult && result;

  result = SetFirstH3Id(firstId);
  finalResult = finalResult && result;
   
  return finalResult; 
}  

//_____________________________________________________________________________
G4bool G4VAnalysisManager::SetFirstH1Id(G4int firstId) 
{
  return fH1HnManager->SetFirstId(firstId);
}  

//_____________________________________________________________________________
G4bool G4VAnalysisManager::SetFirstH2Id(G4int firstId) 
{
  return fH2HnManager->SetFirstId(firstId);
}  

//_____________________________________________________________________________
G4bool G4VAnalysisManager::SetFirstH3Id(G4int firstId) 
{
  return fH3HnManager->SetFirstId(firstId);
}  

//_____________________________________________________________________________
G4bool G4VAnalysisManager::SetFirstProfileId(G4int firstId) 
{
  G4bool finalResult = true;

  G4bool result = SetFirstP1Id(firstId);
  finalResult = finalResult && result;
  
  result = SetFirstP2Id(firstId);
  finalResult = finalResult && result;
   
  return finalResult; 
}  

//_____________________________________________________________________________
G4bool G4VAnalysisManager::SetFirstP1Id(G4int firstId) 
{
  return fP1HnManager->SetFirstId(firstId);
}  

//_____________________________________________________________________________
G4bool G4VAnalysisManager::SetFirstP2Id(G4int firstId) 
{
  return fP2HnManager->SetFirstId(firstId);
}  

//_____________________________________________________________________________
G4bool G4VAnalysisManager::SetFirstNtupleId(G4int firstId) 
{
  return fVNtupleManager->SetFirstId(firstId);
}  

//_____________________________________________________________________________
G4bool G4VAnalysisManager::SetFirstNtupleColumnId(G4int firstId) 
{
  return fVNtupleManager->SetFirstNtupleColumnId(firstId);
}

// Fill methods in .icc

//_____________________________________________________________________________
void  G4VAnalysisManager::SetActivation(G4bool activation) 
{
  fState.SetIsActivation(activation);
}

// GetActivation() in .icc

//_____________________________________________________________________________
G4bool G4VAnalysisManager::IsActive() const
{
// Return true if activation option is selected and any of managers has 
// an activated object.

  return fState.GetIsActivation() && 
         ( fH1HnManager->IsActive() || 
           fH2HnManager->IsActive() || 
           fH3HnManager->IsActive() || 
           fP1HnManager->IsActive() || 
           fP2HnManager->IsActive() );
}  

//_____________________________________________________________________________
G4bool G4VAnalysisManager::IsAscii() const
{
// Return true any of managers has an object with activated ASCII option.

  return ( fH1HnManager->IsAscii() || 
           fH2HnManager->IsAscii() ||
           fH3HnManager->IsAscii() ||
           fP1HnManager->IsAscii() ||
           fP2HnManager->IsAscii() );
}  

//_____________________________________________________________________________
G4bool G4VAnalysisManager::IsPlotting() const
{
// Return true any of managers has an object with activated plotting option.

  return ( fH1HnManager->IsPlotting() || 
           fH2HnManager->IsPlotting() ||
           fH3HnManager->IsPlotting() ||
           fP1HnManager->IsPlotting() ||
           fP2HnManager->IsPlotting() );
}  

//_____________________________________________________________________________
G4int G4VAnalysisManager::GetFirstH1Id() const
{
// Return first H1 id

  return fH1HnManager->GetFirstId();
}  

//_____________________________________________________________________________
G4int G4VAnalysisManager::GetFirstH2Id() const
{
// Return first H2 id

  return fH2HnManager->GetFirstId();
}  

//_____________________________________________________________________________
G4int G4VAnalysisManager::GetFirstH3Id() const
{
// Return first H3 id

  return fH3HnManager->GetFirstId();
}  

//_____________________________________________________________________________
G4int G4VAnalysisManager::GetFirstP1Id() const
{
// Return first P1 id

  return fP1HnManager->GetFirstId();
}  

//_____________________________________________________________________________
G4int G4VAnalysisManager::GetFirstP2Id() const
{
// Return first P2 id

  return fP2HnManager->GetFirstId();
}  

//_____________________________________________________________________________
G4int G4VAnalysisManager::GetFirstNtupleId() const
{
// Return first Ntuple id

  return fVNtupleManager->GetFirstId();
}  

//_____________________________________________________________________________
G4int G4VAnalysisManager::GetFirstNtupleColumnId() const
{
// Return first Ntuple column id

  return fVNtupleManager->GetFirstNtupleColumnId();
}  

//_____________________________________________________________________________
G4int G4VAnalysisManager::GetNofH1s() const
{
  return fH1HnManager->GetNofHns();
}  

//_____________________________________________________________________________
G4int G4VAnalysisManager::GetNofH2s() const
{
  return fH2HnManager->GetNofHns();
}  

//_____________________________________________________________________________
G4int G4VAnalysisManager::GetNofH3s() const
{
  return fH3HnManager->GetNofHns();
}  

//_____________________________________________________________________________
G4int G4VAnalysisManager::GetNofP1s() const
{
  return fP1HnManager->GetNofHns();
}  

//_____________________________________________________________________________
G4int G4VAnalysisManager::GetNofP2s() const
{
  return fP2HnManager->GetNofHns();
}  

//_____________________________________________________________________________
G4int G4VAnalysisManager::GetNofNtuples() const
{
  return fVNtupleManager->GetNofNtuples();
}  

// GetH1Id(), GetH2Id in .icc 

//_____________________________________________________________________________
void  G4VAnalysisManager::SetH1Activation(G4int id, G4bool activation)
{
// Set activation to a given H1 object

  fH1HnManager->SetActivation(id, activation);
}    

//_____________________________________________________________________________
void  G4VAnalysisManager::SetH1Activation(G4bool activation)
{
// Set activation to all H1 objects

  fH1HnManager->SetActivation(activation);
}    

//_____________________________________________________________________________
void  G4VAnalysisManager::SetH1Ascii(G4int id, G4bool ascii)
{
  fH1HnManager->SetAscii(id, ascii);
}    

//_____________________________________________________________________________
void  G4VAnalysisManager::SetH1Plotting(G4int id, G4bool plotting)
{
  fH1HnManager->SetPlotting(id, plotting);
}    

//_____________________________________________________________________________
void  G4VAnalysisManager::SetH2Activation(G4int id, G4bool activation)
{
// Set activation to a given H2 object

  fH2HnManager->SetActivation(id, activation);
}    

//_____________________________________________________________________________
void  G4VAnalysisManager::SetH2Activation(G4bool activation)
{
// Set activation to all H2 objects

  fH2HnManager->SetActivation(activation);
}    

//_____________________________________________________________________________
void  G4VAnalysisManager::SetH2Ascii(G4int id, G4bool ascii)
{
  fH2HnManager->SetAscii(id, ascii);
}

//_____________________________________________________________________________
void  G4VAnalysisManager::SetH2Plotting(G4int id, G4bool plotting)
{
  fH2HnManager->SetPlotting(id, plotting);
}    

//_____________________________________________________________________________
void  G4VAnalysisManager::SetH3Activation(G4int id, G4bool activation)
{
// Set activation to a given H3 object

  fH3HnManager->SetActivation(id, activation);
}    

//_____________________________________________________________________________
void  G4VAnalysisManager::SetH3Activation(G4bool activation)
{
// Set activation to all H3 objects

  fH3HnManager->SetActivation(activation);
}    

//_____________________________________________________________________________
void  G4VAnalysisManager::SetH3Ascii(G4int id, G4bool ascii)
{
  fH3HnManager->SetAscii(id, ascii);
}

//_____________________________________________________________________________
void  G4VAnalysisManager::SetH3Plotting(G4int id, G4bool plotting)
{
  fH3HnManager->SetPlotting(id, plotting);
}

//_____________________________________________________________________________
void  G4VAnalysisManager::SetP1Activation(G4int id, G4bool activation)
{
// Set activation to a given P1 object

  fP1HnManager->SetActivation(id, activation);
}    

//_____________________________________________________________________________
void  G4VAnalysisManager::SetP1Activation(G4bool activation)
{
// Set activation to all P1 objects

  fP1HnManager->SetActivation(activation);
}    

//_____________________________________________________________________________
void  G4VAnalysisManager::SetP1Ascii(G4int id, G4bool ascii)
{
  fP1HnManager->SetAscii(id, ascii);
}

//_____________________________________________________________________________
void  G4VAnalysisManager::SetP1Plotting(G4int id, G4bool plotting)
{
  fP1HnManager->SetPlotting(id, plotting);
}  

//_____________________________________________________________________________
void  G4VAnalysisManager::SetP2Activation(G4int id, G4bool activation)
{
// Set activation to a given P2 object

  fP2HnManager->SetActivation(id, activation);
}    

//_____________________________________________________________________________
void  G4VAnalysisManager::SetP2Activation(G4bool activation)
{
// Set activation to all P2 objects

  fP2HnManager->SetActivation(activation);
}    

//_____________________________________________________________________________
void  G4VAnalysisManager::SetP2Ascii(G4int id, G4bool ascii)
{
  fP2HnManager->SetAscii(id, ascii);
}

//_____________________________________________________________________________
void  G4VAnalysisManager::SetP2Plotting(G4int id, G4bool plotting)
{
  fP2HnManager->SetPlotting(id, plotting);
} 

//_____________________________________________________________________________
void  G4VAnalysisManager::SetNtupleActivation(G4int id, G4bool activation)
{
// Set activation to a given P2 object

  fVNtupleManager->SetActivation(id, activation);
}    

//_____________________________________________________________________________
void  G4VAnalysisManager::SetNtupleActivation(G4bool activation)
{
// Set activation to all P2 objects

  fVNtupleManager->SetActivation(activation);
}    

// Access methods in .icc

//_____________________________________________________________________________
void G4VAnalysisManager::SetVerboseLevel(G4int verboseLevel) 
{
  fState.SetVerboseLevel(verboseLevel);
} 

// GetVerboseLevel() in .icc 

