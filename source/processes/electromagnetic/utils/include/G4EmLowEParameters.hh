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
// GEANT4 Class header file
//
// File name:     G4EmLowEParameters
//
// Author:        Vladimir Ivanchenko
//                  
// Creation date: 06.05.2019
//
// Class Description:
//
// An internal utility class, responsable for keeping parameters
// for low-energy EM and DNA physics processes and models.
//
// It is initialized by the master thread but can be updated 
// at any moment via G4EmParameters interface. It is not assumed
// to be used for a direct initialisation
//
// -------------------------------------------------------------------
//

#ifndef G4EmLowEParameters_h
#define G4EmLowEParameters_h 1

#include "globals.hh"
#include "G4DNAModelSubType.hh"
#include "G4EmFluoDirectory.hh"
#include <vector>

class G4EmLowEParametersMessenger;
class G4VAtomDeexcitation;

class G4EmLowEParameters
{
public:

  explicit G4EmLowEParameters();

  ~G4EmLowEParameters();

  void Initialise();

  // boolean flags
  void SetFluo(G4bool val);
  G4bool Fluo() const;

  G4EmFluoDirectory FluoDirectory() const;

  void SetFluoDirectory(G4EmFluoDirectory val);
  void SetBeardenFluoDir(G4bool val);
  void SetANSTOFluoDir(G4bool val);
  void SetXDB_EADLFluoDir(G4bool val);

  void SetAuger(G4bool val);
  G4bool Auger() const;

  void SetPixe(G4bool val);
  G4bool Pixe() const;

  void SetDeexcitationIgnoreCut(G4bool val);
  G4bool DeexcitationIgnoreCut() const;

  void SetDNAFast(G4bool val);
  G4bool DNAFast() const;

  void SetDNAStationary(G4bool val);
  G4bool DNAStationary() const;

  void SetDNAElectronMsc(G4bool val);
  G4bool DNAElectronMsc() const;

  // integer parameters 
  void SetDNAeSolvationSubType(G4DNAModelSubType val);
  G4DNAModelSubType DNAeSolvationSubType() const;

  // string parameters 
  void SetPIXECrossSectionModel(const G4String&);
  const G4String& PIXECrossSectionModel();

  void SetPIXEElectronCrossSectionModel(const G4String&);
  const G4String& PIXEElectronCrossSectionModel();

  void SetLivermoreDataDir(const G4String&);
  const G4String& LivermoreDataDir();

  // parameters per region or per process 
  void AddMicroElec(const G4String& region);
  const std::vector<G4String>& RegionsMicroElec() const;

  void AddDNA(const G4String& region, const G4String& type);
  const std::vector<G4String>& RegionsDNA() const;
  const std::vector<G4String>& TypesDNA() const;

  void SetDeexActiveRegion(const G4String& region, G4bool fdeex,
			   G4bool fauger, G4bool fpixe);

  // initialisation methods
  void DefineRegParamForDeex(G4VAtomDeexcitation*) const;

  G4EmLowEParameters(G4EmLowEParameters &) = delete;
  G4EmLowEParameters & operator=
  (const G4EmLowEParameters &right) = delete;  

private:

  G4String CheckRegion(const G4String&) const;

  void PrintWarning(G4ExceptionDescription& ed) const;

  G4EmLowEParametersMessenger* theMessenger;

  G4bool fluo;
  G4bool auger;
  G4bool pixe;
  G4bool deexIgnoreCut;

  G4bool dnaFast;
  G4bool dnaStationary;
  G4bool dnaMsc;
  
  G4DNAModelSubType dnaElectronSolvation;

  G4EmFluoDirectory fFluoDirectory;

  G4String namePIXE;
  G4String nameElectronPIXE;
  G4String livDataDir;

  std::vector<G4String>  m_regnamesME;

  std::vector<G4String>  m_regnamesDNA;
  std::vector<G4String>  m_typesDNA;

  std::vector<G4String>  m_regnamesDeex;
  std::vector<G4bool>    m_fluo;
  std::vector<G4bool>    m_auger;
  std::vector<G4bool>    m_pixe;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

