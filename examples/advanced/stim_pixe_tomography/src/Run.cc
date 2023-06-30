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
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "Run.hh"

#include "PrimaryGeneratorAction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
Run::Run(G4bool isP, G4int nbProjs, G4int nbS, G4int nbP, G4int resumeProjIndex)
  : G4Run(),
    isPIXE(isP),
    fProjectionIndex(0),
    fSliceIndex(0),
    fPixelIndex(0),
    fNbProjections(nbProjs),
    fNbSlices(nbS),
    fNbPixels(nbP),
    fResumeProjIndex(resumeProjIndex)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Run::~Run() = default;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::Merge(const G4Run* run)
{
  G4cout << "--------------Merge is called---------------" << G4endl;
  const Run* localRun = static_cast<const Run*>(run);

  for (size_t i = 0; i != localRun->gammaAtExit.size(); i++) {
    gammaAtExit.push_back(localRun->gammaAtExit[i]);
  }
  for (size_t i = 0; i != localRun->gammaAtCreation.size(); i++) {
    gammaAtCreation.push_back(localRun->gammaAtCreation[i]);
  }

  for (size_t i = 0; i != localRun->protonAtExit.size(); i++) {
    protonAtExit.push_back(localRun->protonAtExit[i]);
  }

  G4Run::Merge(run);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::EndOfRun()
{
  // Write X-ray data or proton data into files
  //    G4int runID = GetRunID();
  fProjectionIndex = GetCurrentProjection();
  fSliceIndex = GetCurrentSlice();
  fPixelIndex = GetCurrentPixel();

  RunInfo runInfo((uint8_t)fProjectionIndex, (uint16_t)fSliceIndex, (uint16_t)fPixelIndex);

  if (isPIXE) {
    runInfo.nbParticle = (uint32_t)gammaAtCreation.size();
    WriteFile("GammaAtCreation.dat", runInfo, gammaAtCreation);

    runInfo.nbParticle = (uint32_t)gammaAtExit.size();
    WriteFile("GammaAtExit.dat", runInfo, gammaAtExit);
  }
  else {
    runInfo.nbParticle = (uint32_t)protonAtExit.size();
    WriteFile("ProtonAtExit.dat", runInfo, protonAtExit);
  }

  ClearVecs();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool Run::GetIsPIXE()
{
  return isPIXE;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::FillGammaAtExit(ParticleInfo gammaInfo)
{
  gammaAtExit.push_back(gammaInfo);
}

void Run::FillGammaAtCreation(ParticleInfo gammaInfo)
{
  gammaAtCreation.push_back(gammaInfo);
}

void Run::FillProtonAtExit(ParticleInfo protonInfo)
{
  protonAtExit.push_back(protonInfo);
}

void Run::ClearVecs()
{
  if (gammaAtExit.size()) {
    gammaAtExit.clear();
  }
  if (gammaAtCreation.size()) {
    gammaAtCreation.clear();
  }
  if (protonAtExit.size()) {
    protonAtExit.clear();
  }
}

void Run::WriteFile(const std::string fName, RunInfo& runInfo, std::vector<ParticleInfo>& vec)
{
  std::ofstream ofs;
  if (runID) {
    ofs.open(fName.c_str(),
             std::ios::out | std::ios::app | std::ios::binary);  // appendix
  }
  else {
    // if runId ==0, delete the existing file and create a new file
    ofs.open(fName.c_str(), std::ios::out | std::ios::trunc | std::ios::binary);
  }
  ofs.write((const char*)&runInfo, sizeof(RunInfo));
  if (vec.size()) {
    ofs.write((const char*)vec.data(), sizeof(ParticleInfo) * vec.size());
  }

  ofs.close();
}
G4int Run::GetCurrentProjection()
{
  if (fResumeProjIndex > fNbProjections - 1) {
    G4Exception("fResumeProjIndex", "Run::GetCurrentProjection()", FatalException,
                "To resume a simulation, the start of index of projection must be lower than the "
                "maximal index");
  }
  G4int projIndex = fResumeProjIndex + runID / (fNbSlices * fNbPixels);
  return projIndex;
}
G4int Run::GetCurrentSlice()
{
  G4int remain = runID % (fNbSlices * fNbPixels);
  return remain / fNbPixels;
}

G4int Run::GetCurrentPixel()
{
  G4int remain = runID % (fNbSlices * fNbPixels);
  return remain % fNbPixels;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
