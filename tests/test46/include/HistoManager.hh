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
// $Id: HistoManager.hh,v 1.5 2010-04-13 10:07:39 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//---------------------------------------------------------------------------
//
// ClassName:   HistoManager
//
// Description: Singleton class to hold parameters and build histograms.
//              User cannot access to the constructor.
//              The pointer of the only existing object can be got via
//              HistoManager::GetPointer() static method.
//              The first invokation of this static method makes
//              the singleton object.
//
// Author:      V.Ivanchenko 27/09/00
//
// Modified:
//
//----------------------------------------------------------------------------
//

#ifndef HistoManager_h
#define HistoManager_h 1

#include "globals.hh"
#include "G4Material.hh"
#include "G4Element.hh"

#include "G4SDManager.hh"
#include "HcalSD.hh"
#include "HcalAbsSD.hh"
#include "EcalSD.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class Histo;
class G4Track;
class G4Step;
class G4ParticleDefinition;
class AddHcalAbsorberHit;

class HistoManager
{
public:

  static HistoManager* GetPointer();

private:

  HistoManager();

public: 

  ~HistoManager();

  void bookHisto();

  void BeginOfRun();
  void EndOfRun();
  void BeginOfEvent();
  void EndOfEvent();
  void Fill(G4int id, G4double x, G4double w);

  void AddEcalHit(const G4ParticleDefinition*, G4int, G4double);
  void AddHcalHit(const G4ParticleDefinition*, G4int, G4double);
  void AddHcalAbsorberHit(const G4ParticleDefinition*, G4double);

  void AddStep(const G4ParticleDefinition*);

  void ScoreNewTrack(const G4Track*);

  void SetWorldLength(G4double val); 
  G4double GetWorldLength() const;

  void SetVerbose(G4int val);        
  G4int GetVerbose() const;

  //
  void SetNbins(G4int);
  void SetMaxEnergy(G4double);
  void SetMaxTotEnergy(G4double);
  void SetFactor1(G4double);
  void SetFactor2(G4double);

private:

  static HistoManager* fManager;

  const G4ParticleDefinition* primaryDef;
  const G4ParticleDefinition* neutron;

  G4double beamEnergy;
  G4double worldZ;
  G4double primaryKineticEnergy;
  G4double currentKinEnergy;

  G4double E[25];
  G4double ecal[6];
  G4double edep[6];
  G4double erms[6];
  G4double edeptr[6];
  G4double ermstr[6];
  G4int    stat[6];
  G4double Eecal;
  G4double eecal;
  G4double e9, e25, e19, e125, e925;

  G4double Ehcal;
  G4double Eehcal;
  G4double Eabshcal;

  G4double hcal;
  G4double ehcal;
  G4double abshcal;

  G4double edepSum;
  G4double edepSum2;

  G4double etotSum;
  G4double etotSum2;

  G4double factorEcal;
  G4double factorHcal;
  G4double factorHcal0;
 
  G4int verbose;

  G4double maxEnergy;
  G4double maxTotEnergy;
  G4double maxEnergyAbs;

  G4double m_gamma;
  G4double m_e;
  G4double m_h;
  G4double m_n;
  
  G4int n_evt;
  G4int n_elec;
  G4int n_posit;
  G4int n_gam;
  G4int n_gamph;
  G4int n_gam_tar;
  G4int n_lowe;
  G4int n_step;
  G4int nHisto;
  G4int nBins, nmax;

  Histo* histo;
};

#endif
