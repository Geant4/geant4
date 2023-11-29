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
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:     G4EmLowEParameters
//
// Author:        Vladimir Ivanchenko
//
// Creation date: 07.05.2019
//
// -------------------------------------------------------------------
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4EmLowEParameters.hh"
#include "G4PhysicalConstants.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4VAtomDeexcitation.hh"
#include "G4EmLowEParametersMessenger.hh"
#include "G4RegionStore.hh"
#include "G4Region.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

G4EmLowEParameters::G4EmLowEParameters()
{
  theMessenger = new G4EmLowEParametersMessenger(this);
  Initialise();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

G4EmLowEParameters::~G4EmLowEParameters()
{
  delete theMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void G4EmLowEParameters::Initialise()
{
  fluo = false;
  auger = false;
  pixe = false;
  deexIgnoreCut = false;

  dnaFast = false;
  dnaStationary = false;
  dnaMsc = false;
  dnaElectronSolvation = fMeesungnoen2002eSolvation;

  fFluoDirectory = fluoDefault;
  namePIXE = "Empirical";
  nameElectronPIXE = "Livermore";
  livDataDir = "epics_2017";
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void G4EmLowEParameters::SetFluo(G4bool val)
{
  fluo = val;
}

G4bool G4EmLowEParameters::Fluo() const
{
  return fluo;
}

G4EmFluoDirectory G4EmLowEParameters::FluoDirectory() const
{
  return fFluoDirectory;
}

void G4EmLowEParameters::SetFluoDirectory(G4EmFluoDirectory val)
{
  fFluoDirectory = fluoDefault;
  if(fluoBearden == val) { fFluoDirectory = fluoBearden; }
  else if(fluoANSTO == val) { fFluoDirectory = fluoANSTO; }
  else if(fluoXDB_EADL == val) { fFluoDirectory = fluoXDB_EADL; }
}

void G4EmLowEParameters::SetBeardenFluoDir(G4bool val)
{
  fFluoDirectory = val ? fluoBearden : fluoDefault;
}

void G4EmLowEParameters::SetANSTOFluoDir(G4bool val)
{
  fFluoDirectory = val ? fluoANSTO : fluoDefault;
}

void G4EmLowEParameters::SetXDB_EADLFluoDir(G4bool val)
{
  fFluoDirectory = val ? fluoXDB_EADL : fluoDefault;
}

void G4EmLowEParameters::SetAuger(G4bool val)
{
  auger = val;
  if(val) { fluo = true; }
}

G4bool G4EmLowEParameters::Auger() const
{
  return auger;
}

void G4EmLowEParameters::SetPixe(G4bool val)
{
  pixe = val;
  if(val) { fluo = true; }
}

G4bool G4EmLowEParameters::Pixe() const
{
  return pixe;
}

void G4EmLowEParameters::SetDeexcitationIgnoreCut(G4bool val)
{
  deexIgnoreCut = val;
}

G4bool G4EmLowEParameters::DeexcitationIgnoreCut() const
{
  return deexIgnoreCut;
}

void G4EmLowEParameters::SetDNAFast(G4bool val)
{
  dnaFast = val;
}

G4bool G4EmLowEParameters::DNAFast() const
{
  return dnaFast;
}

void G4EmLowEParameters::SetDNAStationary(G4bool val)
{
  dnaStationary = val;
}

G4bool G4EmLowEParameters::DNAStationary() const
{
  return dnaStationary;
}

void G4EmLowEParameters::SetDNAElectronMsc(G4bool val)
{
  dnaMsc = val;
}

G4bool G4EmLowEParameters::DNAElectronMsc() const
{
  return dnaMsc;
}

void G4EmLowEParameters::SetDNAeSolvationSubType(G4DNAModelSubType val)
{
  dnaElectronSolvation = val;
}

G4DNAModelSubType G4EmLowEParameters::DNAeSolvationSubType() const
{
  return dnaElectronSolvation;
}

void G4EmLowEParameters::SetPIXECrossSectionModel(const G4String& sss)
{
  namePIXE = sss;
}

const G4String& G4EmLowEParameters::PIXECrossSectionModel()
{
  return namePIXE;
}

void G4EmLowEParameters::SetPIXEElectronCrossSectionModel(const G4String& sss)
{
  nameElectronPIXE = sss;
}

const G4String& G4EmLowEParameters::PIXEElectronCrossSectionModel()
{
  return nameElectronPIXE;
}

void G4EmLowEParameters::SetLivermoreDataDir(const G4String& sss)
{
  livDataDir = sss;
}

const G4String& G4EmLowEParameters::LivermoreDataDir()
{
  return livDataDir;
}

void G4EmLowEParameters::PrintWarning(G4ExceptionDescription& ed) const
{
  G4Exception("G4EmLowEParameters", "em0044", JustWarning, ed);
}

G4String G4EmLowEParameters::CheckRegion(const G4String& reg) const
{
  G4String r = reg;
  if(r == "" || r == "world" || r == "World") {
    r = "DefaultRegionForTheWorld";
  }
  return r;
}

void G4EmLowEParameters::AddMicroElec(const G4String& region)
{
  G4String r = CheckRegion(region);
  std::size_t nreg =  m_regnamesME.size();
  for(std::size_t i=0; i<nreg; ++i) {
    if(r == m_regnamesME[i]) { return; }
  }
  m_regnamesME.push_back(r);
}

const std::vector<G4String>& G4EmLowEParameters::RegionsMicroElec() const
{
  return m_regnamesME;
}

void G4EmLowEParameters::AddDNA(const G4String& region, const G4String& type)
{
  G4String r = CheckRegion(region);
  std::size_t nreg =  m_regnamesDNA.size();
  for(std::size_t i=0; i<nreg; ++i) {
    if(r == m_regnamesDNA[i]) { return; }
  }
  m_regnamesDNA.push_back(r);
  m_typesDNA.push_back(type);
}

const std::vector<G4String>& G4EmLowEParameters::RegionsDNA() const
{
  return m_regnamesDNA;
}

const std::vector<G4String>& G4EmLowEParameters::TypesDNA() const
{
  return m_typesDNA;
}

void 
G4EmLowEParameters::SetDeexActiveRegion(const G4String& region, G4bool fdeex,
                                        G4bool fauger, G4bool fpixe)
{
  if(fdeex) { fluo = true; }
  G4String r = CheckRegion(region);
  std::size_t nreg =  m_regnamesDeex.size();
  if(0 == nreg && r != "DefaultRegionForTheWorld") {
    m_regnamesDeex.push_back("DefaultRegionForTheWorld");
    m_fluo.push_back(false);
    m_auger.push_back(false);
    m_pixe.push_back(false);
    nreg = 1;
  }
  for(std::size_t i=0; i<nreg; ++i) {
    if(r == m_regnamesDeex[i]) { 
      m_fluo[i] = fdeex;
      m_auger[i]= fauger;
      m_pixe[i] = fpixe;
      return; 
    }
  }
  m_regnamesDeex.push_back(r);
  m_fluo.push_back(fdeex);
  m_auger.push_back(fauger);
  m_pixe.push_back(fpixe);
}

void G4EmLowEParameters::DefineRegParamForDeex(G4VAtomDeexcitation* ptr) const
{
  std::size_t n = m_regnamesDeex.size();
  for(std::size_t i=0; i<n; ++i) {
    ptr->SetDeexcitationActiveRegion(m_regnamesDeex[i],
				     m_fluo[i], m_auger[i], m_pixe[i]);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
