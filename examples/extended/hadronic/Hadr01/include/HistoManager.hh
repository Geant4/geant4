//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
// $Id: HistoManager.hh,v 1.2 2006-06-06 19:48:38 vnivanch Exp $
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
