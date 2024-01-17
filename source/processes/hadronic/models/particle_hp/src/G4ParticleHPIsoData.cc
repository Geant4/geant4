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
// particle_hp -- source file
// J.P. Wellisch, Nov-1996
// A prototype of the low energy neutron transport model.
//
// 080901 Avoiding troubles which caused by G4PhysicsVecotor of length 0 by T. Koi
//
// P. Arce, June-2014 Conversion neutron_hp to particle_hp
//
#include "G4ParticleHPIsoData.hh"

#include "G4Neutron.hh"
#include "G4ParticleHPDataUsed.hh"
#include "G4ParticleHPManager.hh"
#include "G4Alpha.hh"
#include "G4Deuteron.hh"
#include "G4He3.hh"
#include "G4Proton.hh"
#include "G4Triton.hh"
//#include <stdlib.h>
#include <fstream>

void G4ParticleHPIsoData::FillChannelData(G4ParticleHPVector* aBuffer)
{
  if (theChannelData != nullptr) {
    G4Exception("G4ParticleHPIsoData::FillChannelData","hadhp02",
                FatalException, "Inconsistency: the data uploaded next time");
  }
  theChannelData = new G4ParticleHPVector;
  for (G4int i = 0; i < aBuffer->GetVectorLength(); ++i) {
    theChannelData->SetPoint(i, aBuffer->GetPoint(i));
  }
  theChannelData->Hash();
}

G4bool G4ParticleHPIsoData::Init(G4int A, G4int Z, G4int M, G4double abun, G4String dirName,
                                 G4String aFSType)
{
  theChannelData = nullptr;
  G4double abundance = abun / 100.;
  G4String filename;
  G4bool result = true;
  G4ParticleHPDataUsed aFile = theNames.GetName(A, Z, M, dirName, aFSType, result);
  filename = aFile.GetName();

  std::istringstream theChannel(filename, std::ios::in);
  auto man = G4ParticleHPManager::GetInstance();
  man->GetDataStream(filename, theChannel);

#ifdef G4PHPDEBUG
  if (man->GetDEBUG())
    G4cout << "G4ParticleHPIsoData::Init = " << filename << " " << A << " " << Z << G4endl;
#endif

  if (Z == 1 && (aFile.GetZ() != Z || aFile.GetA() != A)) {
    if (man->GetDEBUG())
      G4cout << "Skipped = " << filename << " " << A << " " << Z << G4endl;
    // 080901 TKDB No more necessary below protection, cross sections set to 0 
  }
  if (!theChannel) { /*theChannel.close()*/
    return false;
  }
  // accommodating deficiencie of some compilers
  if (theChannel.eof()) { /*theChannel.close()*/
    return false;
  }
  G4int dummy;
  theChannel >> dummy >> dummy;
  theChannelData = new G4ParticleHPVector;
  G4int nData;
  theChannel >> nData;
  theChannelData->Init(theChannel, nData, CLHEP::eV, abundance * CLHEP::barn);
  return result;
}

void G4ParticleHPIsoData::Init(G4int A, G4int Z, G4int M, G4double abun,
                               G4ParticleDefinition* projectile, const char*)
// fill PhysicsVector for this Isotope
{
  auto man = G4ParticleHPManager::GetInstance();

  G4String baseName = man->GetParticleHPPath(projectile);

  G4String dirName;
  if (projectile == G4Neutron::Neutron()) {
    dirName = baseName + "/Fission";
    if (Z > 87)  // TK Modifed for ENDF VII.0
    {
      Init(A, Z, M, abun, dirName, "/CrossSection");
    }
    else {
      theChannelData = new G4ParticleHPVector;
    }
    theFissionData = theChannelData;
    theChannelData = nullptr;  // fast fix for double delete; revisit later. @@@@@@@

    dirName = baseName + "/Capture";
    Init(A, Z, M, abun, dirName, "/CrossSection");
    theCaptureData = theChannelData;
    theChannelData = nullptr;

    dirName = baseName + "/Elastic";
    Init(A, Z, M, abun, dirName, "/CrossSection");
    theElasticData = theChannelData;
    theChannelData = nullptr;
  }

  dirName = baseName + "/Inelastic";
  Init(A, Z, M, abun, dirName, "/CrossSection");
  theInelasticData = theChannelData;
  theChannelData = nullptr;
}

G4String G4ParticleHPIsoData::GetName(G4int A, G4int Z, G4String base, G4String rest)
{
  G4bool dbool;
  return (theNames.GetName(A, Z, base, rest, dbool)).GetName();
}
