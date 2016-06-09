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
// $Id: HistoManager.hh,v 1.3 2006/06/29 17:23:42 gunter Exp $
// GEANT4 tag $Name: geant4-08-01 $
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
// 04.06.2006 Adoptation of hadr01 (V.Ivanchenko)
//
//----------------------------------------------------------------------------
//

#ifndef HistoManager_h
#define HistoManager_h 1

#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class Histo;
class G4Track;
class G4Step;
class G4ParticleDefinition;

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

  void ScoreNewTrack(const G4Track*);
  void AddTargetStep(const G4Step*);
  void AddLeakingParticle(const G4Track*);

  void SetTargetLength(G4double val) {length  = val;};
  void SetNumberOfSlices(G4int val)  {nSlices = val;};
  void SetNumberOfBinsE(G4int val)   {nBinsE  = val;};

  G4double Length()         const    {return length;};
  G4int    NumberOfSlices() const    {return nSlices;};

  void SetVerbose(G4int val);        
  G4int GetVerbose() const           {return verbose;};

  G4double CurrentKinEnergy()        {return currentKinEnergy;};
  const G4ParticleDefinition* CurrentDefinition() 
                                     {return currentDef;};

private:

  static HistoManager* fManager;

  const G4ParticleDefinition* primaryDef;
  const G4ParticleDefinition* currentDef;
  const G4ParticleDefinition* neutron;
  G4double currentKinEnergy;
 
  G4int verbose;
  G4int nBinsE;
  G4int nSlices;

  G4double beamEnergy;
  G4double length;
  G4double absZ0;
  G4double primaryKineticEnergy;

  G4int n_evt;
  G4int n_elec;
  G4int n_posit;
  G4int n_gam;
  G4int n_prot_leak;
  G4int n_pion_leak;
  G4int n_cpions;
  G4int n_pi0;
  G4int n_kaons;
  G4int n_muons;
  G4int n_ions;
  G4int n_deut;
  G4int n_alpha;
  G4int n_neutron;
  G4int n_proton;
  G4int n_aproton;
  G4int n_neu_forw;
  G4int n_neu_leak;
  G4int n_neu_back;
  G4int n_step;
  G4int nHisto;

  Histo* histo;
};

#endif
