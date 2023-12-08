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
// G4MicroElecCapture.cc,
//                 2011/08/29 A.Valentin, M. Raine are with CEA [a]
//                 2020/05/20 P. Caron, C. Inguimbert are with ONERA [b] 
//                            Q. Gibaru is with CEA [a], ONERA [b] and CNES [c]
//                            M. Raine and D. Lambert are with CEA [a]
//
// A part of this work has been funded by the French space agency(CNES[c])
// [a] CEA, DAM, DIF - 91297 ARPAJON, France
// [b] ONERA - DPHY, 2 avenue E.Belin, 31055 Toulouse, France
// [c] CNES, 18 av.E.Belin, 31401 Toulouse CEDEX, France
//
// Based on the following publications
// - A.Valentin, M. Raine, 
//   Inelastic cross-sections of low energy electrons in silicon
//   for the simulation of heavy ion tracks with the Geant4-DNA toolkit,
//   NSS Conf. Record 2010, pp. 80-85
//   https://doi.org/10.1109/NSSMIC.2010.5873720
//
// - A.Valentin, M. Raine, M.Gaillardin, P.Paillet
//   Geant4 physics processes for microdosimetry simulation:
//   very low energy electromagnetic models for electrons in Silicon,
//   https://doi.org/10.1016/j.nimb.2012.06.007
//   NIM B, vol. 288, pp. 66-73, 2012, part A
//   heavy ions in Si, NIM B, vol. 287, pp. 124-129, 2012, part B
//   https://doi.org/10.1016/j.nimb.2012.07.028
//
// - M. Raine, M. Gaillardin, P. Paillet
//   Geant4 physics processes for silicon microdosimetry simulation: 
//   Improvements and extension of the energy-range validity up to 10 GeV/nucleon
//   NIM B, vol. 325, pp. 97-100, 2014
//   https://doi.org/10.1016/j.nimb.2014.01.014
//
// - J. Pierron, C. Inguimbert, M. Belhaj, T. Gineste, J. Puech, M. Raine
//   Electron emission yield for low energy electrons: 
//   Monte Carlo simulation and experimental comparison for Al, Ag, and Si
//   Journal of Applied Physics 121 (2017) 215107. 
//   https://doi.org/10.1063/1.4984761
//
// - P. Caron,
//   Study of Electron-Induced Single-Event Upset in Integrated Memory Devices
//   PHD, 16th October 2019
//
// - Q.Gibaru, C.Inguimbert, P.Caron, M.Raine, D.Lambert, J.Puech, 
//   Geant4 physics processes for microdosimetry and secondary electron emission simulation : 
//   Extension of MicroElec to very low energies and new materials
//   NIM B, 2020, in review.
//
//----------------------------------------------------------------------------
//
// ClassName: G4MicroElecCapture derivated from G4ElectronCapture (V Ivanchenko)
//
// Description: The process to kill particles to save CPU
//
// Author: C. Inguimbert 31 january 2022 derivated from G4ElectronCapture (V.Ivanchenko 31 August 2010)
//
//----------------------------------------------------------------------------
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4MicroElecCapture.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleDefinition.hh"
#include "G4Step.hh"
#include "G4PhysicalConstants.hh"
#include "G4Track.hh"
#include "G4Region.hh"
#include "G4RegionStore.hh"
#include "G4Electron.hh"
#include "G4Pow.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4MicroElecCapture::G4MicroElecCapture(const G4String& regName, G4double ekinlim)
  : G4VDiscreteProcess("MicroElecCapture", fElectromagnetic), kinEnergyThreshold(ekinlim),
    regionName(regName), region(0)
{
  if(regName == "" || regName == "world")
  { 
    regionName = "DefaultRegionForTheWorld";
  }
  isInitialised = false;
  pParticleChange = &fParticleChange;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4MicroElecCapture::~G4MicroElecCapture() 
{
  for (auto pos = tableWF.cbegin(); pos != tableWF.cend(); ++pos)
  {
    G4MicroElecMaterialStructure* table = pos->second;
    delete table;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4MicroElecCapture::SetKinEnergyLimit(G4double val)
{
  kinEnergyThreshold = val;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4MicroElecCapture::BuildPhysicsTable(const G4ParticleDefinition&)
{
  region = (G4RegionStore::GetInstance())->GetRegion(regionName);
 // if(region && verboseLevel > 0) {
  G4cout << "### G4MicroElecCapture: Tracking cut E(MeV) = " 
         << kinEnergyThreshold/MeV << " is assigned to " << regionName 
         << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool G4MicroElecCapture::IsApplicable(const G4ParticleDefinition&)
{
  return true;
}

void G4MicroElecCapture::Initialise()
{
  if (isInitialised) { return; }

  G4ProductionCutsTable* theCoupleTable = G4ProductionCutsTable::GetProductionCutsTable();
  G4int numOfCouples = (G4int)theCoupleTable->GetTableSize();
  G4cout << numOfCouples << G4endl;

  for (G4int i = 0; i < numOfCouples; ++i)
  {
    const G4Material* material = theCoupleTable->GetMaterialCutsCouple(i)->GetMaterial();

    G4cout << "G4Capture, Material " << i + 1 << " / "
           << numOfCouples << " : " << material->GetName() << G4endl;
    if (material->GetName() == "Vacuum")
    {
      tableWF[material->GetName()] = 0;
      continue;
    }
    G4String mat = material->GetName();
    G4MicroElecMaterialStructure* str = new G4MicroElecMaterialStructure(mat);
    tableWF[mat] = str;
  }
  isInitialised = true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VParticleChange* G4MicroElecCapture::PostStepDoIt(const G4Track& aTrack, 
                                                    const G4Step&)
{
  if (!isInitialised) { Initialise(); }

  G4String mat = aTrack.GetMaterial()->GetName();
  G4int Ztarget = ((*(aTrack.GetMaterial()->GetElementVector()))[0])->GetZasInt();
  G4int Atarget = ((*(aTrack.GetMaterial()->GetElementVector()))[0])->GetAtomicMassAmu();
  G4double Nbelements = aTrack.GetMaterial()->GetNumberOfElements();
  G4double moleculeMass = aTrack.GetMaterial()->GetMassOfMolecule() / amu;
  auto FractionMass = aTrack.GetMaterial()->GetFractionVector();
  G4int Zinc = aTrack.GetParticleDefinition()->GetAtomicNumber();
  G4int Ainc = aTrack.GetParticleDefinition()->GetAtomicMass();
  G4String IncPartName = aTrack.GetParticleDefinition()->GetParticleName();
  G4double NIEdep = 0.0;

  for (G4int i = 0; i < Nbelements; ++i)
  {
    Ztarget = ((*(aTrack.GetMaterial()->GetElementVector()))[i])->GetZasInt();
    Atarget = ((*(aTrack.GetMaterial()->GetElementVector()))[i])->GetAtomicMassAmu();
    NIEdep = NIEdep + moleculeMass*FractionMass[i] / Atarget*G_Lindhard_Rob(aTrack.GetKineticEnergy(), Zinc, Ainc, Ztarget, Atarget);
  }

  WorkFunctionTable::iterator matWF;
  matWF = tableWF.find(mat);

  if (matWF == tableWF.end())
  {
    G4String str = "Material ";
    str += mat + " not found!";
    G4Exception("G4MicroElecCapture::PostStepGPIL", "em0002",
                FatalException, str);
    return nullptr;
  }
  else
  {
    G4MicroElecMaterialStructure* str = matWF->second;
    pParticleChange->Initialize(aTrack);
    pParticleChange->ProposeTrackStatus(fStopAndKill);

    G4double InitE = str->GetEnergyGap() + str->GetInitialEnergy();

    if (IncPartName == "e-")
    {
      // metals = Non ionizing deposited energy = 0.0
      if (((str->GetEnergyGap()) / eV)<(0.001))
      {
        pParticleChange->ProposeNonIonizingEnergyDeposit(0.0);
        pParticleChange->ProposeLocalEnergyDeposit(aTrack.GetKineticEnergy());
      }
      else // MicroElec materials Non ionizing deposited energy different from zero
      {
        G4int c = (G4int)((aTrack.GetKineticEnergy()) / (InitE));
        pParticleChange->ProposeNonIonizingEnergyDeposit(aTrack.GetKineticEnergy() - InitE*c);
        pParticleChange->ProposeLocalEnergyDeposit(aTrack.GetKineticEnergy());
      }
    }
    else
    {
      if ((IncPartName == "Genericion") || (IncPartName == "alpha")
       || (IncPartName == "He3") || (IncPartName == "deuteron")
       || (IncPartName == "triton") || (IncPartName == "proton"))
      {
        pParticleChange->ProposeNonIonizingEnergyDeposit(NIEdep);
        pParticleChange->ProposeLocalEnergyDeposit(aTrack.GetKineticEnergy());
      }
      else
      {
        pParticleChange->ProposeNonIonizingEnergyDeposit(0.0);
        pParticleChange->ProposeLocalEnergyDeposit(aTrack.GetKineticEnergy());
      }
    }
  } // matWF == tableWF.end())
		
  fParticleChange.SetProposedKineticEnergy(0.0);
  return pParticleChange;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4MicroElecCapture::GetMeanFreePath(const G4Track& aTrack, G4double,
                                             G4ForceCondition*)
{											
  G4String material = aTrack.GetMaterial()->GetName();
  // test particle type in order to applied the capture to both electrons, protons and heavy ions

  if ((aTrack.GetParticleDefinition()->GetParticleName()) == "e-")
  {
    if (material != "G4_ALUMINUM_OXIDE" && material != "G4_SILICON_DIOXIDE"
     && material != "G4_BORON_NITRIDE")
    {
      return DBL_MAX;
    }
    G4double    S = 0;
    G4double    y = 0;
    if (material == "G4_ALUMINUM_OXIDE")
    {
      S = 1 * (1 / nm);
      y = 0.25 * (1 / eV);
    }
    if (material == "G4_SILICON_DIOXIDE")
    {
      S = 0.3 * (1 / nm);
      y = 0.2 * (1 / eV);
    }
    if (material == "G4_BORON_NITRIDE")
    {
      S = 0 * (1 / nm);
      y = 1 * (1 / eV);
    }

    G4double P = S * G4Exp(-y * aTrack.GetKineticEnergy());
    if (P <= 0) { return DBL_MAX; }
           else { return 1 / P; }
  }
  else return DBL_MAX;
  
}    

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4MicroElecCapture::G_Lindhard_Rob(G4double Trecoil, G4int Zrecoil, G4int Arecoil, G4int Zcible, G4int Acible)
{
  G4double Lind =0.0;

  if (Arecoil <= 0 || Zrecoil == 0)
  {
    Lind = 0.0;
  }
  else
  {
    G4double El = 30.724 * Zcible * Zrecoil
                * std::pow((G4Pow::GetInstance()->Z23(Zcible) + G4Pow::GetInstance()->Z23(Zrecoil)), 0.5)
                * (Arecoil + Acible) / Acible;

    // multiplication by 1e6 to change El from eV to MeV
    G4double e = Trecoil / (El * CLHEP::eV);
    G4double Fl = (0.0793 * G4Pow::GetInstance()->Z23(Zrecoil) * std::pow(Zcible, 0.5) * std::pow((Arecoil + Acible), 1.5))
                / (std::pow((G4Pow::GetInstance()->Z23(Zcible) + G4Pow::GetInstance()->Z23(Zrecoil)), 3. / 4.) * std::pow(Arecoil, 3. / 2.) * std::pow(Acible, 1. / 2.));

    Lind = 1. / (1 + Fl * (3.4008 * std::pow(e, 1. / 6.) + 0.40244 * std::pow(e, 3. / 4.) + e));

    // to get the energie that go into displacement
    Lind = Lind * Trecoil;
  }
  return Lind;                                                                   
}
