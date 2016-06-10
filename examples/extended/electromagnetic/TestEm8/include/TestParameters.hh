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
/// \file electromagnetic/TestEm8/include/HistoManager.hh
/// \brief Definition of the HistoManager class
//
// $Id: HistoManager.hh 78549 2014-01-07 09:42:35Z gcosmo $
//
//---------------------------------------------------------------------------
//
// ClassName:   HistoManager
//
// Description: Singleton class to make analysis and build histograms.
//              User cannot access to the constructor.
//              The pointer of the only existing object can be got via
//              HistoManager::GetPointer() static method.
//              The first invokation of this static method makes
//              the singleton object.
//
// Author:      V.Ivanchenko 01.09.2010
//
//----------------------------------------------------------------------------
//

#ifndef TestParameters_h
#define TestParameters_h 1

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "globals.hh"
#include "G4DataVector.hh"
#include "G4StatDouble.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

//class RunAction;
//class G4Step;

class TestParameters
{

public:
  // With description

  static TestParameters* GetPointer();

private:

  TestParameters();

public: // Without description

  ~TestParameters();

  void SetMaxEnergy(G4double value);

  G4double GetMaxEnergy() const;

  void SetNumberBins(G4int value);

  G4int GetNumberBins() const;

  void SetNumberBinsCluster(G4int value);

  G4int GetNumberBinsCluster() const;

  void SetEnergyPerChannel(G4double value);

  G4double GetFactorALICE() const;

  void SetPositionZ(G4double val);

  G4double GetPositionZ() const;

private:

  // MEMBERS
  static TestParameters* fManager;

  //  G4int fNHisto;
  //  G4int fVerbose;

  G4double fMaxEnergy;
  G4double fFactorALICE;
  G4double fPositionZ;

  G4int fBinsE; 
  G4int fBinsCluster;
};

#endif
