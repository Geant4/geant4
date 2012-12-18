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
// $Id: HistoManager.hh,v 1.1 2008-05-27 15:15:29 antoni Exp $
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
  void BeginOfEvent();
  void EndOfEvent();
  void Fill(G4int id, G4double x, G4double w);

  void ScoreNewTrack(const G4Track*);
  void AddLeakingParticle(const G4Track*);

  void SetTargetLength(G4double val)            {length  = val;};
  void SetNumberOfBinsE(G4int val)              {nbins   = val;};

  G4double Length()         const               {return length;};

  void SetVerbose(G4int val);        
  G4int GetVerbose() const                      {return verbose;};

private:

  static HistoManager* fManager;

  G4int IndexTheta(G4double theta);  

  const G4ParticleDefinition* primaryDef;
  const G4ParticleDefinition* neutron;

  G4double edepEvt;
  G4double edepSum;
  G4double edepSum2;
  G4double beamEnergy;
  G4double length;
  G4double primaryKineticEnergy;
  G4double currentKinEnergy;

  G4double angle[6];
  G4double dangle;
  G4double emin1, emin2, emax1, emax2;
  G4double de1, de2;
 
  G4int verbose;
  G4int nbins, nbins1, nbins2;

  G4int n_evt;
  G4int n_step;
  G4int nHisto;

  Histo* histo;
};

#endif
