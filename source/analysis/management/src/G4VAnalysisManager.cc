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

// Author: Ivana Hrivnacova, 09/07/2013  (ivana@ipno.in2p3.fr)

#include "G4VAnalysisManager.hh"
#include "G4AnalysisMessenger.hh"
#include "G4AnalysisUtilities.hh"
#include "G4HnManager.hh"
#include "G4VNtupleManager.hh"
#include "G4VNtupleFileManager.hh"
#include "G4VFileManager.hh"
#include "G4NtupleBookingManager.hh"
#include "G4Threading.hh"
#include "G4AutoLock.hh"

using namespace G4Analysis;

namespace {
  //Mutex to lock master manager when merging histograms
  G4Mutex registerWorkerMutex = G4MUTEX_INITIALIZER;
}

namespace {

//_____________________________________________________________________________
void NtupleMergingWarning(std::string_view className,
                          std::string_view functionName,
                          const G4String& outputType)
{
  G4Analysis::Warn(
    "Ntuple merging is not available with " + outputType + " output.\n" +
    "Setting is ignored.", className, functionName);
}

}

//
// ctor/dtor
//

//_____________________________________________________________________________
G4VAnalysisManager::G4VAnalysisManager(const G4String& type)
 : fState(type, ! G4Threading::IsWorkerThread())
{
  fMessenger = std::make_unique<G4AnalysisMessenger>(this);
  fNtupleBookingManager = std::make_shared<G4NtupleBookingManager>(fState);

  // Set master/worker instances
  // used only in "FromUI" functions
  if ( ! G4Threading::IsWorkerThread() ) {
    fgMasterInstance = this;
  }
  else {
    if (fgMasterInstance != nullptr) {
      G4AutoLock lock(&registerWorkerMutex);
      fgMasterInstance->fWorkerManagers.push_back(this);
      lock.unlock();
    }
  }
}

//_____________________________________________________________________________
G4VAnalysisManager::~G4VAnalysisManager() = default;

//
// private methods
//

//_____________________________________________________________________________
G4bool G4VAnalysisManager::WriteFromUI()
{  
// Write is performed on workers first, then on master

  if (! fState.GetIsMaster() ) return true;

  auto result = true;

  // Process first all workers
  for (auto workerManger : fWorkerManagers) {
    // Update G4Threading 
    auto g4ThreadingValue = G4Threading::G4GetThreadId();
    G4Threading::G4SetThreadId(workerManger->fState.GetThreadId());

    result &= workerManger->Write();

    // Set G4Threading back
    G4Threading::G4SetThreadId(g4ThreadingValue);
  }

  // Process write on master
  result &= Write();

  return result;
}

//_____________________________________________________________________________
G4bool G4VAnalysisManager::CloseFileFromUI(G4bool reset)
{
// Close file is performed on workers first, then on master

  if (! fState.GetIsMaster() ) return true;

  auto result = true;

  // Process first all workers
  for (auto workerManger : fWorkerManagers) {
    // Update G4Threading 
    auto g4ThreadingValue = G4Threading::G4GetThreadId();
    G4Threading::G4SetThreadId(workerManger->fState.GetThreadId());

    result &= workerManger->CloseFile(reset);
 
    // Set G4Threading back
    G4Threading::G4SetThreadId(g4ThreadingValue);
  }

  // Process write on master
  result &= CloseFile(reset);

  return result;
}

//_____________________________________________________________________________
G4bool G4VAnalysisManager::ResetFromUI()
{  
// Reset file is performed on workers first, then on master

  if (! fState.GetIsMaster() ) return true;

  auto result = true;

  // Process first all workers
  for (auto workerManger : fWorkerManagers) {
    // Update G4Threading 
    auto g4ThreadingValue = G4Threading::G4GetThreadId();
    G4Threading::G4SetThreadId(workerManger->fState.GetThreadId());

    result &= workerManger->Reset();
 
    // Set G4Threading back
    G4Threading::G4SetThreadId(g4ThreadingValue);
  }

  // Process write on master
  result &= Reset();

  return result;
}

//
// protected methods
//

//_____________________________________________________________________________
void G4VAnalysisManager::SetH1Manager(G4VTBaseHnManager<kDim1>* h1Manager)
{
  fVH1Manager.reset(h1Manager);
  fH1HnManager = h1Manager->GetHnManager();
  if (fVFileManager != nullptr ) fH1HnManager->SetFileManager(fVFileManager);
}

//_____________________________________________________________________________
void G4VAnalysisManager::SetH2Manager(G4VTBaseHnManager<kDim2>* h2Manager)
{
  fVH2Manager.reset(h2Manager);
  fH2HnManager = h2Manager->GetHnManager();
  if (fVFileManager != nullptr ) fH2HnManager->SetFileManager(fVFileManager);
}

//_____________________________________________________________________________
void G4VAnalysisManager::SetH3Manager(G4VTBaseHnManager<kDim3>* h3Manager)
{
  fVH3Manager.reset(h3Manager);
  fH3HnManager = h3Manager->GetHnManager();
  if (fVFileManager != nullptr ) fH3HnManager->SetFileManager(fVFileManager);
}

//_____________________________________________________________________________
void G4VAnalysisManager::SetP1Manager(G4VTBaseHnManager<kDim2>* p1Manager)
{
  fVP1Manager.reset(p1Manager);
  fP1HnManager = p1Manager->GetHnManager();
  if (fVFileManager != nullptr ) fP1HnManager->SetFileManager(fVFileManager);
}

//_____________________________________________________________________________
void G4VAnalysisManager::SetP2Manager(G4VTBaseHnManager<kDim3>* p2Manager)
{
  fVP2Manager.reset(p2Manager);
  fP2HnManager = p2Manager->GetHnManager();
  if (fVFileManager != nullptr ) fP2HnManager->SetFileManager(fVFileManager);
}

//_____________________________________________________________________________
void G4VAnalysisManager::SetNtupleManager(std::shared_ptr<G4VNtupleManager> ntupleManager)
{
  fVNtupleManager = std::move(ntupleManager);
  fVNtupleManager->SetFirstId(fNtupleBookingManager->GetFirstId());
  fVNtupleManager->SetFirstNtupleColumnId(fNtupleBookingManager->GetFirstNtupleColumnId());
}

//_____________________________________________________________________________
void G4VAnalysisManager::SetNtupleFileManager(
  std::shared_ptr<G4VNtupleFileManager> ntupleFileManager)
{
  fVNtupleFileManager = std::move(ntupleFileManager);
}

//_____________________________________________________________________________
void G4VAnalysisManager::SetFileManager(std::shared_ptr<G4VFileManager> fileManager)
{
  fVFileManager = fileManager;

  if ( fH1HnManager != nullptr ) fH1HnManager->SetFileManager(fileManager);
  if ( fH2HnManager != nullptr ) fH2HnManager->SetFileManager(fileManager);
  if ( fH3HnManager != nullptr ) fH3HnManager->SetFileManager(fileManager);
  if ( fP1HnManager != nullptr ) fP1HnManager->SetFileManager(fileManager);
  if ( fP2HnManager != nullptr ) fP2HnManager->SetFileManager(std::move(fileManager));
}

//_____________________________________________________________________________
G4bool G4VAnalysisManager::WriteAscii(const G4String& fileName)
{
  // Do not write on workers
  if ( ! fState.GetIsMaster() ) return true;

  auto result = true;

  // Replace or add file extension .ascii
  G4String name(fileName);
  if (name.find('.') != std::string::npos) {
    name.erase(name.find('.'), name.length());
  }
  name.append(".ascii");

  Message(kVL3, "write ASCII", "file", name);

  std::ofstream output(name, std::ios::out);
  if ( ! output ) {
    Warn("Cannot open file. File name is not defined.",
      fkClass, "WriteAscii");
    return false;
  }
  output.setf( std::ios::scientific, std::ios::floatfield );

  result &= fVH1Manager->WriteOnAscii(output);
  result &= fVH2Manager->WriteOnAscii(output);
  result &= fVH3Manager->WriteOnAscii(output);
  result &= fVP1Manager->WriteOnAscii(output);
  result &= fVP2Manager->WriteOnAscii(output);

  Message(kVL1, "write ASCII", "file", name, result);

  return result;
}

//_____________________________________________________________________________
std::shared_ptr<G4VFileManager>
G4VAnalysisManager::GetFileManager(const G4String& fileName)
{
  // Check if file type corresponds the manager output type
  G4String extension = GetExtension(fileName);
  if ((extension.size() != 0u) && extension != GetFileType()) {
    Warn(
      "The file extension differs from " + GetFileType() + " output type.\n" +
      GetFileType() + " output type will be used.",
      fkClass, "GetFileManager");
  }

  return fVFileManager;
}

//
// public methods
//

//_____________________________________________________________________________
G4bool G4VAnalysisManager::OpenFile(const G4String& fileName)
{
  // Protection against opening file twice
  // (Seems to happen when opening file via UI command after the first run)
  if (IsOpenFile()) {
    // G4cout << "Skipping OpenFile. File is already open" << G4endl;
    return true;
  }

  if ( fileName != "" ) {
    return OpenFileImpl(fileName);
  }
  if (fVFileManager->GetFileName() == "") {
    Warn("Cannot open file. File name is not defined.", fkClass, "OpenFile");
    return false;
  }
  return OpenFileImpl(fVFileManager->GetFileName());
}

//_____________________________________________________________________________
G4bool G4VAnalysisManager::Write()
{
  auto result = true;

  result &= WriteImpl();
  if ( IsPlotting() ) {
    result &= PlotImpl();
  }

  // Increment cycle number
  fState.IncrementCycle();

  return result;
}

//_____________________________________________________________________________
G4bool G4VAnalysisManager::CloseFile(G4bool reset)
{
  auto result = CloseFileImpl(reset);

  // Notify about new cycle
  fState.ResetCycle();
  if (fVNtupleManager != nullptr) {
    fVNtupleManager->SetNewCycle(false);
  }

  return result;
}

//_____________________________________________________________________________
G4bool G4VAnalysisManager::Reset()
{
  // Notify about new cycle
  // (as reset causes deleting ntuples)
  if (fVNtupleManager != nullptr) {
    fVNtupleManager->SetNewCycle(true);
  }

  return ResetImpl();
}

//_____________________________________________________________________________
void G4VAnalysisManager::Clear()
{
  Message(kVL4, "clear", "all data");

  // clear hntools objects
  ClearImpl();

  // clear remaining data
  fNtupleBookingManager->ClearData();
  if ( fVNtupleManager != nullptr ) fVNtupleManager->Clear();
  if ( fVFileManager != nullptr ) fVFileManager->Clear();

  Message(kVL1, "clear", "all data");
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
  std::array<G4HnDimension, kDim1> bins = {
    G4HnDimension(nbins, xmin, xmax)};
  std::array<G4HnDimensionInformation, kDim1> info = {
    G4HnDimensionInformation(unitName, fcnName, binSchemeName) };

  return fVH1Manager->Create(name, title, bins, info);
}

//_____________________________________________________________________________
G4int G4VAnalysisManager::CreateH1(const G4String& name,  const G4String& title,
                               const std::vector<G4double>& edges,
                               const G4String& unitName, const G4String& fcnName)
{
  std::array<G4HnDimension, kDim1> bins = {
    G4HnDimension(edges)};
  std::array<G4HnDimensionInformation, kDim1> info = {
    G4HnDimensionInformation(unitName, fcnName, "user")};

  return fVH1Manager->Create(name, title, bins, info);
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
  std::array<G4HnDimension, kDim2> bins = {
    G4HnDimension(nxbins, xmin, xmax),
    G4HnDimension(nybins, ymin, ymax) };
  std::array<G4HnDimensionInformation, kDim2> info = {
    G4HnDimensionInformation(xunitName, xfcnName, xbinSchemeName),
    G4HnDimensionInformation(yunitName, yfcnName, ybinSchemeName)};

  return fVH2Manager->Create(name, title, bins, info);
}

//_____________________________________________________________________________
G4int G4VAnalysisManager::CreateH2(const G4String& name,  const G4String& title,
                               const std::vector<G4double>& xedges,
                               const std::vector<G4double>& yedges,
                               const G4String& xunitName, const G4String& yunitName,
                               const G4String& xfcnName, const G4String& yfcnName)

{
  std::array<G4HnDimension, kDim2> bins = {
    G4HnDimension(xedges), G4HnDimension(yedges)};
  std::array<G4HnDimensionInformation, kDim2> info = {
    G4HnDimensionInformation(xunitName, xfcnName, "user"),
    G4HnDimensionInformation(yunitName, yfcnName, "user")};

  return fVH2Manager->Create(name, title, bins, info);
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
  std::array<G4HnDimension, kDim3> bins = {
    G4HnDimension(nxbins, xmin, xmax),
    G4HnDimension(nybins, ymin, ymax),
    G4HnDimension(nzbins, zmin, zmax)};
  std::array<G4HnDimensionInformation, kDim3> info = {
    G4HnDimensionInformation(xunitName, xfcnName, xbinSchemeName),
    G4HnDimensionInformation(yunitName, yfcnName, ybinSchemeName),
    G4HnDimensionInformation(zunitName, zfcnName, zbinSchemeName)};

  return fVH3Manager->Create(name, title, bins, info);
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
  std::array<G4HnDimension, kDim3> bins = {
    G4HnDimension(xedges),  G4HnDimension(yedges), G4HnDimension(zedges) };
  std::array<G4HnDimensionInformation, kDim3> info = {
    G4HnDimensionInformation(xunitName, xfcnName, "user"),
    G4HnDimensionInformation(yunitName, yfcnName, "user"),
    G4HnDimensionInformation(zunitName, zfcnName, "user")};

  return fVH3Manager->Create(name, title, bins, info);
}

//_____________________________________________________________________________
G4bool G4VAnalysisManager::SetH1(G4int id,
                                G4int nbins, G4double xmin, G4double xmax,
                                const G4String& unitName, const G4String& fcnName,
                                const G4String& binSchemeName)
{
  std::array<G4HnDimension, kDim1> bins = {
    G4HnDimension(nbins, xmin, xmax)};
  std::array<G4HnDimensionInformation, kDim1> info = {
    G4HnDimensionInformation(unitName, fcnName, binSchemeName) };

  return fVH1Manager->Set(id, bins, info);
}

//_____________________________________________________________________________
G4bool G4VAnalysisManager::SetH1(G4int id,
                                const std::vector<G4double>& edges,
                                const G4String& unitName, const G4String& fcnName)
{
  std::array<G4HnDimension, kDim1> bins = {
    G4HnDimension(edges)};
  std::array<G4HnDimensionInformation, kDim1> info = {
    G4HnDimensionInformation(unitName, fcnName, "user")};

  return fVH1Manager->Set(id, bins, info);
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
  std::array<G4HnDimension, kDim2> bins = {
    G4HnDimension(nxbins, xmin, xmax),
    G4HnDimension(nybins, ymin, ymax) };
  std::array<G4HnDimensionInformation, kDim2> info = {
    G4HnDimensionInformation(xunitName, xfcnName, xbinSchemeName),
    G4HnDimensionInformation(yunitName, yfcnName, ybinSchemeName)};

  return fVH2Manager->Set(id, bins, info);
}

//_____________________________________________________________________________
G4bool G4VAnalysisManager::SetH2(G4int id,
                                const std::vector<G4double>& xedges,
                                const std::vector<G4double>& yedges,
                                const G4String& xunitName, const G4String& yunitName,
                                const G4String& xfcnName, const G4String& yfcnName)
{
  std::array<G4HnDimension, kDim2> bins = {
    G4HnDimension(xedges), G4HnDimension(yedges)};
  std::array<G4HnDimensionInformation, kDim2> info = {
    G4HnDimensionInformation(xunitName, xfcnName, "user"),
    G4HnDimensionInformation(yunitName, yfcnName, "user")};

  return fVH2Manager->Set(id, bins, info);
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
  std::array<G4HnDimension, kDim3> bins = {
    G4HnDimension(nxbins, xmin, xmax),
    G4HnDimension(nybins, ymin, ymax),
    G4HnDimension(nzbins, zmin, zmax)};
  std::array<G4HnDimensionInformation, kDim3> info = {
    G4HnDimensionInformation(xunitName, xfcnName, xbinSchemeName),
    G4HnDimensionInformation(yunitName, yfcnName, ybinSchemeName),
    G4HnDimensionInformation(zunitName, zfcnName, zbinSchemeName)};

  return fVH3Manager->Set(id, bins, info);
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
  std::array<G4HnDimension, kDim3> bins = {
    G4HnDimension(xedges), G4HnDimension(yedges), G4HnDimension(zedges) };
  std::array<G4HnDimensionInformation, kDim3> info = {
    G4HnDimensionInformation(xunitName, xfcnName, "user"),
    G4HnDimensionInformation(yunitName, yfcnName, "user"),
    G4HnDimensionInformation(zunitName, zfcnName, "user")};

  return fVH3Manager->Set(id, bins, info);
}

//_____________________________________________________________________________
G4bool G4VAnalysisManager::ScaleH1(G4int id, G4double factor)
{
  return fVH1Manager->Scale(id, factor);
}

//_____________________________________________________________________________
G4bool G4VAnalysisManager::ScaleH2(G4int id, G4double factor)
{
  return fVH2Manager->Scale(id, factor);
}

//_____________________________________________________________________________
G4bool G4VAnalysisManager::ScaleH3(G4int id, G4double factor)
{
  return fVH3Manager->Scale(id, factor);
}

//_____________________________________________________________________________
G4int G4VAnalysisManager::CreateP1(const G4String& name,  const G4String& title,
                               G4int nbins, G4double xmin, G4double xmax,
                               G4double ymin, G4double ymax,
                               const G4String& xunitName, const G4String& yunitName,
                               const G4String& xfcnName, const G4String& yfcnName,
                               const G4String& xbinSchemeName)
{
  std::array<G4HnDimension, kDim2> bins = {
    G4HnDimension(nbins, xmin, xmax),
    G4HnDimension(0, ymin, ymax) };
  std::array<G4HnDimensionInformation, kDim2> info = {
    G4HnDimensionInformation(xunitName, xfcnName, xbinSchemeName),
    G4HnDimensionInformation(yunitName, yfcnName)};

  return fVP1Manager->Create(name, title, bins, info);
}

//_____________________________________________________________________________
G4int G4VAnalysisManager::CreateP1(const G4String& name,  const G4String& title,
                               const std::vector<G4double>& edges,
                               G4double ymin, G4double ymax,
                               const G4String& xunitName, const G4String& yunitName,
                               const G4String& xfcnName, const G4String& yfcnName)
{
  std::array<G4HnDimension, kDim2> bins = {
    G4HnDimension(edges), G4HnDimension(0, ymin, ymax)};
  std::array<G4HnDimensionInformation, kDim2> info = {
    G4HnDimensionInformation(xunitName, xfcnName),
    G4HnDimensionInformation(yunitName, yfcnName)};

  return fVP1Manager->Create(name, title, bins, info);
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
  std::array<G4HnDimension, kDim3> bins = {
    G4HnDimension(nxbins, xmin, xmax),
    G4HnDimension(nybins, ymin, ymax),
    G4HnDimension(0, zmin, zmax)};
  std::array<G4HnDimensionInformation, kDim3> info = {
    G4HnDimensionInformation(xunitName, xfcnName, xbinSchemeName),
    G4HnDimensionInformation(yunitName, yfcnName, ybinSchemeName),
    G4HnDimensionInformation(zunitName, zfcnName)};

  return fVP2Manager->Create(name, title, bins, info);
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
  std::array<G4HnDimension, kDim3> bins = {
    G4HnDimension(xedges),  G4HnDimension(yedges), G4HnDimension(0, zmin, zmax)};
  std::array<G4HnDimensionInformation, kDim3> info = {
    G4HnDimensionInformation(xunitName, xfcnName),
    G4HnDimensionInformation(yunitName, yfcnName),
    G4HnDimensionInformation(zunitName, zfcnName)};

  return fVP2Manager->Create(name, title, bins, info);
}

//_____________________________________________________________________________
G4bool G4VAnalysisManager::SetP1(G4int id,
                                G4int nbins, G4double xmin, G4double xmax,
                                G4double ymin, G4double ymax,
                                const G4String& xunitName, const G4String& yunitName,
                                const G4String& xfcnName, const G4String& yfcnName,
                                const G4String& xbinSchemeName)
{
  std::array<G4HnDimension, kDim2> bins = {
    G4HnDimension(nbins, xmin, xmax),
    G4HnDimension(0, ymin, ymax) };
  std::array<G4HnDimensionInformation, kDim2> info = {
    G4HnDimensionInformation(xunitName, xfcnName, xbinSchemeName),
    G4HnDimensionInformation(yunitName, yfcnName)};

  return fVP1Manager->Set(id, bins, info);
}

//_____________________________________________________________________________
G4bool G4VAnalysisManager::SetP1(G4int id,
                                const std::vector<G4double>& edges,
                                G4double ymin, G4double ymax,
                                const G4String& xunitName, const G4String& yunitName,
                                const G4String& xfcnName, const G4String& yfcnName)
{
  std::array<G4HnDimension, kDim2> bins = {
    G4HnDimension(edges), G4HnDimension(0, ymin, ymax)};
  std::array<G4HnDimensionInformation, kDim2> info = {
    G4HnDimensionInformation(xunitName, xfcnName),
    G4HnDimensionInformation(yunitName, yfcnName)};

  return fVP1Manager->Set(id, bins, info);
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
  std::array<G4HnDimension, kDim3> bins = {
    G4HnDimension(nxbins, xmin, xmax),
    G4HnDimension(nybins, ymin, ymax),
    G4HnDimension(0, zmin, zmax)};
  std::array<G4HnDimensionInformation, kDim3> info = {
    G4HnDimensionInformation(xunitName, xfcnName, xbinSchemeName),
    G4HnDimensionInformation(yunitName, yfcnName, ybinSchemeName),
    G4HnDimensionInformation(zunitName, zfcnName)};

  return fVP2Manager->Set(id, bins, info);
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
  std::array<G4HnDimension, kDim3> bins = {
    G4HnDimension(xedges), G4HnDimension(yedges), G4HnDimension(0, zmin, zmax)};
  std::array<G4HnDimensionInformation, kDim3> info = {
    G4HnDimensionInformation(xunitName, xfcnName),
    G4HnDimensionInformation(yunitName, yfcnName),
    G4HnDimensionInformation(zunitName, zfcnName)};

  return fVP2Manager->Set(id, bins, info);
}

//_____________________________________________________________________________
G4bool G4VAnalysisManager::ScaleP1(G4int id, G4double factor)
{
  return fVP1Manager->Scale(id, factor);
}

//_____________________________________________________________________________
G4bool G4VAnalysisManager::ScaleP2(G4int id, G4double factor)
{
  return fVP2Manager->Scale(id, factor);
}

//_____________________________________________________________________________
G4int G4VAnalysisManager::CreateNtuple(const G4String& name,
                                          const G4String& title)
{
  return fNtupleBookingManager->CreateNtuple(name, title);
}

//_____________________________________________________________________________
G4int G4VAnalysisManager::CreateNtupleIColumn(const G4String& name)
{
  return fNtupleBookingManager->CreateNtupleIColumn(name, nullptr);
}

//_____________________________________________________________________________
G4int G4VAnalysisManager::CreateNtupleFColumn(const G4String& name)
{
  return fNtupleBookingManager->CreateNtupleFColumn(name, nullptr);
}

//_____________________________________________________________________________
G4int G4VAnalysisManager::CreateNtupleDColumn(const G4String& name)
{
  return fNtupleBookingManager->CreateNtupleDColumn(name, nullptr);
}

//_____________________________________________________________________________
G4int G4VAnalysisManager::CreateNtupleSColumn(const G4String& name)
{
  return fNtupleBookingManager->CreateNtupleSColumn(name, nullptr);
}

//_____________________________________________________________________________
G4int G4VAnalysisManager::CreateNtupleIColumn(const G4String& name,
                                              std::vector<int>& vector)
{
  return fNtupleBookingManager->CreateNtupleIColumn(name, &vector);
}

//_____________________________________________________________________________
G4int G4VAnalysisManager::CreateNtupleFColumn(const G4String& name,
                                              std::vector<float>& vector)
{
  return fNtupleBookingManager->CreateNtupleFColumn(name, &vector);
}

//_____________________________________________________________________________
G4int G4VAnalysisManager::CreateNtupleDColumn(const G4String& name,
                                              std::vector<double>& vector)
{
  return fNtupleBookingManager->CreateNtupleDColumn(name, &vector);
}

//_____________________________________________________________________________
G4int G4VAnalysisManager::CreateNtupleSColumn(const G4String& name,
                                              std::vector<std::string>& vector)
{
  return fNtupleBookingManager->CreateNtupleSColumn(name, &vector);
}

//_____________________________________________________________________________
void G4VAnalysisManager::FinishNtuple()
{
  auto ntupleBooking = fNtupleBookingManager->FinishNtuple();

  if ( fVNtupleManager ) {
    fVNtupleManager->CreateNtuple(ntupleBooking);
  }
}

//_____________________________________________________________________________
void G4VAnalysisManager::SetNtupleMerging(G4bool /*mergeNtuples*/,
                   G4int /*nofReducedNtupleFiles*/)
{
// The function is overridden in the managers which supports ntuple merging
// Here we give just a warning that the feature is not available.

  NtupleMergingWarning(fkClass, "SetNtupleMerging", GetType());
}

//_____________________________________________________________________________
void G4VAnalysisManager::SetNtupleRowWise(G4bool /*rowWise*/,
                                          G4bool /*rowMode*/)
{
// The function is overridden in the managers which supports ntuple merging
// Here we give just a warning that the feature is not available.

  NtupleMergingWarning(fkClass, "SetNtupleRowWise", GetType());
}

//_____________________________________________________________________________
void G4VAnalysisManager::SetBasketSize(unsigned int /*basketSize*/)
{
// The function is overridden in the managers which supports ntuple merging
// Here we give just a warning that the feature is not available.

  NtupleMergingWarning(fkClass, "SetBasketSize", GetType());
}

//_____________________________________________________________________________
void G4VAnalysisManager::SetBasketEntries(unsigned int /*basketEntries*/)
{
// The function is overridden in the managers which supports ntuple merging
// Here we give just a warning that the feature is not available.

  NtupleMergingWarning(fkClass, "SetBasketEntries", GetType());
}

//_____________________________________________________________________________
G4int G4VAnalysisManager::CreateNtupleIColumn(G4int ntupleId,
                                              const G4String& name)
{
  return fNtupleBookingManager->CreateNtupleIColumn(ntupleId, name, nullptr);
}

//_____________________________________________________________________________
G4int G4VAnalysisManager::CreateNtupleFColumn(G4int ntupleId,
                                              const G4String& name)
{
  return fNtupleBookingManager->CreateNtupleFColumn(ntupleId, name, nullptr);
}


//_____________________________________________________________________________
G4int G4VAnalysisManager::CreateNtupleDColumn(G4int ntupleId,
                                              const G4String& name)
{
  return fNtupleBookingManager->CreateNtupleDColumn(ntupleId, name, nullptr);
}

//_____________________________________________________________________________
G4int G4VAnalysisManager::CreateNtupleSColumn(G4int ntupleId,
                                              const G4String& name)
{
  return fNtupleBookingManager->CreateNtupleSColumn(ntupleId, name, nullptr);
}

//_____________________________________________________________________________
G4int G4VAnalysisManager::CreateNtupleIColumn(G4int ntupleId,
                                              const G4String& name,
                                              std::vector<int>& vector)
{
  return fNtupleBookingManager->CreateNtupleIColumn(ntupleId, name, &vector);
}

//_____________________________________________________________________________
G4int G4VAnalysisManager::CreateNtupleFColumn(G4int ntupleId,
                                              const G4String& name,
                                              std::vector<float>& vector)
{
  return fNtupleBookingManager->CreateNtupleFColumn(ntupleId, name, &vector);
}

//_____________________________________________________________________________
G4int G4VAnalysisManager::CreateNtupleDColumn(G4int ntupleId,
                                              const G4String& name,
                                              std::vector<double>& vector)
{
  return fNtupleBookingManager->CreateNtupleDColumn(ntupleId, name, &vector);
}

//_____________________________________________________________________________
G4int G4VAnalysisManager::CreateNtupleSColumn(G4int ntupleId,
                                              const G4String& name,
                                              std::vector<std::string>& vector)
{
  return fNtupleBookingManager->CreateNtupleSColumn(ntupleId, name, &vector);
}

//_____________________________________________________________________________
void G4VAnalysisManager::FinishNtuple(G4int ntupleId)
{
  auto ntupleBooking = fNtupleBookingManager->FinishNtuple(ntupleId);

  if ( fVNtupleManager ) {
    fVNtupleManager->CreateNtuple(ntupleBooking);
  }
}

//_____________________________________________________________________________
G4bool G4VAnalysisManager::SetFirstHistoId(G4int firstId)
{
  auto result = true;

  result &= SetFirstH1Id(firstId);
  result &= SetFirstH2Id(firstId);
  result &= SetFirstH3Id(firstId);

  return result;
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
  auto result = true;

  result &= SetFirstP1Id(firstId);
  result &= SetFirstP2Id(firstId);

  return result;
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
  auto result = true;

  result &= fNtupleBookingManager->SetFirstId(firstId);
  if ( fVNtupleManager ) {
    result &= fVNtupleManager->SetFirstId(firstId);
  }

  return result;
}

//_____________________________________________________________________________
G4bool G4VAnalysisManager::SetFirstNtupleColumnId(G4int firstId)
{
  auto result = true;

  result &= fNtupleBookingManager->SetFirstNtupleColumnId(firstId);
  if ( fVNtupleManager ) {
    result &= fVNtupleManager->SetFirstNtupleColumnId(firstId);
  }

  return result;
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

  return fNtupleBookingManager->GetFirstId();
}

//_____________________________________________________________________________
G4int G4VAnalysisManager::GetFirstNtupleColumnId() const
{
// Return first Ntuple column id

  return fNtupleBookingManager->GetFirstNtupleColumnId();
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
  if (fVNtupleManager != nullptr) {
    return fVNtupleManager->GetNofNtuples();
  }

  return 0;
}

// GetH1Id(), GetH2Id in .icc

//_____________________________________________________________________________
G4bool G4VAnalysisManager::ListH1(G4bool onlyIfActive) const
{
  return fVH1Manager->List(G4cout, onlyIfActive);
}

//_____________________________________________________________________________
G4bool G4VAnalysisManager::ListH2(G4bool onlyIfActive) const
{
  return fVH2Manager->List(G4cout, onlyIfActive);
}

//_____________________________________________________________________________
G4bool G4VAnalysisManager::ListH3(G4bool onlyIfActive) const
{
  return fVH3Manager->List(G4cout, onlyIfActive);
}

//_____________________________________________________________________________
G4bool G4VAnalysisManager::ListP1(G4bool onlyIfActive) const
{
  return fVP1Manager->List(G4cout, onlyIfActive);
}

//_____________________________________________________________________________
G4bool G4VAnalysisManager::ListP2(G4bool onlyIfActive) const
{
  return fVP2Manager->List(G4cout, onlyIfActive);
}

//_____________________________________________________________________________
G4bool G4VAnalysisManager::ListNtuple(G4bool onlyIfActive) const
{
  if (fVNtupleManager != nullptr) {
    return fVNtupleManager->List(G4cout, onlyIfActive);
  }
  return false;
}

//_____________________________________________________________________________
G4bool G4VAnalysisManager::List(G4bool onlyIfActive) const
{
  auto result = true;
  result &= ListH1(onlyIfActive);
  result &= ListH2(onlyIfActive);
  result &= ListH3(onlyIfActive);
  result &= ListP1(onlyIfActive);
  result &= ListP2(onlyIfActive);
  result &= ListNtuple(onlyIfActive);

  return result;
}

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
void  G4VAnalysisManager::SetH1FileName(G4int id, const G4String& fileName)
{
  fH1HnManager->SetFileName(id, fileName);
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
void  G4VAnalysisManager::SetH2FileName(G4int id, const G4String& fileName)
{
  fH2HnManager->SetFileName(id, fileName);
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
void  G4VAnalysisManager::SetH3FileName(G4int id, const G4String& fileName)
{
  fH3HnManager->SetFileName(id, fileName);
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
void  G4VAnalysisManager::SetP1FileName(G4int id, const G4String& fileName)
{
  fP1HnManager->SetFileName(id, fileName);
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
void  G4VAnalysisManager::SetP2FileName(G4int id, const G4String& fileName)
{
  fP2HnManager->SetFileName(id, fileName);
}

//_____________________________________________________________________________
void  G4VAnalysisManager::SetNtupleActivation(G4int id, G4bool activation)
{
// Set activation to a given ntuple object

  fNtupleBookingManager->SetActivation(id, activation);
  if ( fVNtupleManager ) {
    fVNtupleManager->SetActivation(id, activation);
  }
}

//_____________________________________________________________________________
void  G4VAnalysisManager::SetNtupleActivation(G4bool activation)
{
// Set activation to all ntuple objects

  fNtupleBookingManager->SetActivation(activation);
  if ( fVNtupleManager ) {
    fVNtupleManager->SetActivation(activation);
  }
}

//_____________________________________________________________________________
void  G4VAnalysisManager::SetNtupleFileName(G4int id, const G4String& fileName)
{
// Set activation to a given P2 object

  fNtupleBookingManager->SetFileName(id, fileName);
}

//_____________________________________________________________________________
void  G4VAnalysisManager::SetNtupleFileName(const G4String& fileName)
{
// Set activation to all P2 objects

  fNtupleBookingManager->SetFileName(fileName);
}

// Access methods in .icc

//_____________________________________________________________________________
void G4VAnalysisManager::SetVerboseLevel(G4int verboseLevel)
{
  fState.SetVerboseLevel(verboseLevel);
}

// GetVerboseLevel() in .icc

