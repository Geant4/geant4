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

#ifndef Run_h
#define Run_h 1

#include "G4Electron.hh"
#include "G4Gamma.hh"
#include "G4Positron.hh"
#include "G4Run.hh"
#include "globals.hh"

#include <vector>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//
//// struct for keeping particle information
// energy + momentum
//
struct ParticleInfo
{
  G4float energy;
  G4float mx;
  G4float my;
  G4float mz;
  ParticleInfo(G4float e, G4float vx, G4float vy, G4float vz) : energy(e), mx(vx), my(vy), mz(vz) {}
};

// struct ParticleInfo
// {
//   G4float energy;
//   G4float mx;
//   G4float my;
//   G4float mz;
//   G4float x;
//   G4float y;
//   G4float z;
//   ParticleInfo(G4float e, G4float vx, G4float vy, G4float vz,
//                         G4float x0, G4float y0, G4float z0)
//     : energy(e)
//     , mx(vx)
//     , my(vy)
//     , mz(vz)
//     , x(x0)
//     , y(y0)
//     , z(z0)
//   {}
// };

struct RunInfo
{
  uint8_t projectionIndex;  // 1 byte
  uint16_t sliceIndex;  //
  uint16_t pixelIndex;
  uint32_t nbParticle;  // 4 bytes
  RunInfo(uint8_t proI, uint16_t sI, uint16_t pI, uint32_t nb = 0)
    : projectionIndex(proI), sliceIndex(sI), pixelIndex(pI), nbParticle(nb)
  {}
};

class Run : public G4Run
{
public:
  Run(G4bool isP, G4int nbProjs, G4int nbSlices, G4int nbPixels, G4int resumeProjectionIndex);
  ~Run() override;

public:
  void Merge(const G4Run*) override;
  void EndOfRun();

  G4bool GetIsPIXE();
  G4int GetCurrentProjection();
  G4int GetCurrentSlice();
  G4int GetCurrentPixel();

  void FillGammaAtExit(ParticleInfo gammaInfo);
  void FillGammaAtCreation(ParticleInfo gammaInfo);
  void FillProtonAtExit(ParticleInfo protonInfo);

  void ClearVecs();
  void WriteFile(const std::string fName, RunInfo& runInfo, std::vector<ParticleInfo>& vec);

private:
  G4bool isPIXE;
  std::vector<ParticleInfo> gammaAtExit;  // gamma at exit for a run
  std::vector<ParticleInfo> gammaAtCreation;  // gamma at creation for a run
  std::vector<ParticleInfo> protonAtExit;  // proton at exit for a run

  G4int fProjectionIndex;
  G4int fSliceIndex;
  G4int fPixelIndex;

  G4int fNbProjections;
  G4int fNbSlices;
  G4int fNbPixels;
  G4int fResumeProjIndex;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
