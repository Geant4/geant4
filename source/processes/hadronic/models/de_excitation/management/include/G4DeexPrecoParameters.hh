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
// $Id: G4DeexPrecoParameters.hh 68028 2013-03-13 13:48:15Z gcosmo $
//
// 15.03.2016 V.Ivanchenko 
//
// List of parameters of the pre-compound model
// and the deexcitation module
//

#ifndef G4DeexPrecoParameters_h
#define G4DeexPrecoParameters_h 1

#include "globals.hh"
#include "G4Threading.hh"

enum G4DeexChannelType
{
  fEvaporation = 0,
  fGEM,
  fCombined
};

class G4StateManager;

class G4DeexPrecoParameters
{
public:

  explicit G4DeexPrecoParameters();

  void SetDefaults();

  // printing
  std::ostream& StreamInfo(std::ostream& os) const;
  void Dump() const;
  friend std::ostream& operator<< (std::ostream& os, 
				   const G4DeexPrecoParameters&);

  // inline access methods 

  inline G4double GetLevelDensity() const;

  inline G4double GetR0() const;

  inline G4double GetTransitionsR0() const;

  inline G4double GetFermiEnergy() const;

  inline G4double GetPrecoLowEnergy() const;

  inline G4double GetPhenoFactor() const;

  inline G4double GetMinExcitation() const;

  inline G4double GetMaxLifeTime() const;

  inline G4double GetMinExPerNucleounForMF() const;

  inline G4int GetInternalConversionID() const;

  inline G4int GetMinZForPreco() const;

  inline G4int GetMinAForPreco() const;

  inline G4int GetPrecoModelType() const;

  inline G4int GetDeexModelType() const;

  inline G4bool NeverGoBack() const;

  inline G4bool UseSoftCutoff() const;

  inline G4bool UseCEM() const;

  inline G4bool UseGNASH() const;

  inline G4bool UseHETC() const;

  inline G4bool UseAngularGen() const;

  inline G4bool UseFilesNEW() const;

  inline G4bool CorrelatedGamma() const;

  inline G4bool StoreAllLevels() const;

  inline G4DeexChannelType GetDeexChannelsType() const;

  // Set methods 

  void SetLevelDensity(G4double);

  void SetR0(G4double);

  void SetTransitionsR0(G4double);

  void SetFermiEnergy(G4double);

  void SetPrecoLowEnergy(G4double);

  void SetPhenoFactor(G4double);

  void SetMinExcitation(G4double);

  void SetMaxLifeTime(G4double);

  void SetMinExPerNucleounForMF(G4double);

  void SetMinEForMultiFrag(G4double);

  void SetMinZForPreco(G4int);

  void SetMinAForPreco(G4int);

  void SetPrecoModelType(G4int);

  void SetDeexModelType(G4int);

  void SetNeverGoBack(G4bool);

  void SetUseSoftCutoff(G4bool);

  void SetUseCEM(G4bool);

  void SetUseGNASH(G4bool);

  void SetUseHETC(G4bool);

  void SetUseAngularGen(G4bool);

  void SetUseFilesNEW(G4bool);

  void SetCorrelatedGamma(G4bool);

  void SetStoreAllLevels(G4bool);

  void SetDeexChannelsType(G4DeexChannelType);

private:

  G4bool IsLocked() const;

  G4DeexPrecoParameters(const G4DeexPrecoParameters & right) = delete;  
  const G4DeexPrecoParameters& operator=(const G4DeexPrecoParameters &right) = delete;
  G4bool operator==(const G4DeexPrecoParameters &right) const = delete;
  G4bool operator!=(const G4DeexPrecoParameters &right) const = delete;

  G4StateManager* fStateManager;

  // Level density parameter
  G4double fLevelDensity;

  // Nuclear radius r0 
  G4double fR0;

  // Nuclear radius r0 for transitions
  G4double fTransitionsR0;

  // Fermi energy level
  G4double fFermiEnergy;

  // Low-energy excitation per nucleon limit 
  G4double fPrecoLowEnergy;

  // Preco phenomenological factor
  G4double fPhenoFactor;

  // Excitation handler
  G4double fMinExcitation;
  G4double fMaxLifeTime;

  // Multi-fragmentation model
  G4double fMinExPerNucleounForMF;

  // Cross section type
  G4int fPrecoType;
  G4int fDeexType;

  // Internal conversion model ID
  G4int fInternalConversionID;

  // Preco model
  G4int fMinZForPreco;
  G4int fMinAForPreco;

  // Preco flags
  G4bool fNeverGoBack;
  G4bool fUseSoftCutoff;
  G4bool fUseCEM;
  G4bool fUseGNASH;
  G4bool fUseHETC;
  G4bool fUseAngularGen;

  // Deex flags
  G4bool fUseLongFiles;
  G4bool fCorrelatedGamma;
  G4bool fStoreAllLevels;

  // type of a set of e-exitation channels
  G4DeexChannelType fDeexChannelType;   

#ifdef G4MULTITHREADED
  static G4Mutex deexPrecoMutex;
#endif
};

inline G4double G4DeexPrecoParameters::GetLevelDensity() const
{ 
  return fLevelDensity; 
}
 
inline G4double G4DeexPrecoParameters::GetR0() const
{ 
  return fR0; 
}

inline G4double G4DeexPrecoParameters::GetTransitionsR0() const
{ 
  return fTransitionsR0; 
}

inline G4double G4DeexPrecoParameters::GetFermiEnergy() const
{ 
  return fFermiEnergy; 
}

inline G4double G4DeexPrecoParameters::GetPrecoLowEnergy() const
{ 
  return fPrecoLowEnergy; 
}

inline G4double G4DeexPrecoParameters::GetPhenoFactor() const
{ 
  return fPhenoFactor; 
}

inline G4double G4DeexPrecoParameters::GetMinExcitation() const
{
  return fMinExcitation;
}

inline G4double G4DeexPrecoParameters::GetMaxLifeTime() const
{
  return fMaxLifeTime;
}

inline G4double G4DeexPrecoParameters::GetMinExPerNucleounForMF() const
{
  return fMinExPerNucleounForMF;
}

inline G4int G4DeexPrecoParameters::GetInternalConversionID() const
{
  return fInternalConversionID;
}

inline G4int G4DeexPrecoParameters::GetMinZForPreco() const
{
  return fMinZForPreco;
}

inline G4int G4DeexPrecoParameters::GetMinAForPreco() const
{
  return fMinAForPreco;
}

inline G4int G4DeexPrecoParameters::GetPrecoModelType() const
{
  return fPrecoType;
}

inline G4int G4DeexPrecoParameters::GetDeexModelType() const
{
  return fDeexType;
}

inline G4bool G4DeexPrecoParameters::NeverGoBack() const
{
  return fNeverGoBack;
}

inline G4bool G4DeexPrecoParameters::UseSoftCutoff() const
{
  return fUseSoftCutoff;
}

inline G4bool G4DeexPrecoParameters::UseCEM() const
{
  return fUseCEM;
}

inline G4bool G4DeexPrecoParameters::UseGNASH() const
{
  return fUseGNASH;
}

inline G4bool G4DeexPrecoParameters::UseHETC() const
{
  return fUseHETC;
}

inline G4bool G4DeexPrecoParameters::UseAngularGen() const
{
  return fUseAngularGen;
}

inline G4bool G4DeexPrecoParameters::UseFilesNEW() const
{
  return fUseLongFiles;
}

inline G4bool G4DeexPrecoParameters::CorrelatedGamma() const
{
  return fCorrelatedGamma;
}

inline G4bool G4DeexPrecoParameters::StoreAllLevels() const
{
  return fStoreAllLevels;
}

inline G4DeexChannelType G4DeexPrecoParameters::GetDeexChannelsType() const
{
  return fDeexChannelType;
}

#endif
