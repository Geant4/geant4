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
//
// $Id: GammaRayTelAnalysisManager.hh,v 1.3 2001-07-11 09:56:56 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// ------------------------------------------------------------
//      GEANT 4 class header file
//      CERN Geneva Switzerland
//     
//
//      ------------ GammaRayTelAnalysisManager  ------
//           by R.Giannitrapani, F. Longo & G.Santin (30 nov 2000)
//
// ************************************************************

#ifdef G4ANALYSIS_USE
#ifndef GammaRayTelAnalysisManager_h
#define GammaRayTelAnalysisManager_h 1

#include "G4VAnalysisManager.hh"

#include "globals.hh"
#include "g4std/vector"
#include "G4ThreeVector.hh"

class GammaRayTelAnalysisMessenger;
class GammaRayTelDetectorConstruction;
class IHistogramFactory;
class IHistogram1D;
class IHistogram2D;
class IPlotter;
class IVectorFactory;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class GammaRayTelAnalysisManager: public G4VAnalysisManager
{
public:
  GammaRayTelAnalysisManager(GammaRayTelDetectorConstruction*);
  virtual ~GammaRayTelAnalysisManager();
  
public:
  G4bool RegisterAnalysisSystem(G4VAnalysisSystem*);
  IHistogramFactory* GetHistogramFactory(const G4String&);
  void Store(IHistogram* = 0, const G4String& = "");
  void Plot(IHistogram* = 0);
  void InsertPositionXZ(double x, double z);
  void InsertPositionYZ(double y, double z);
  void InsertEnergy(double en);
  void InsertHits(int nplane);
  void BeginOfRun();
  void EndOfRun(G4int n);
  void EndOfEvent(G4int flag);

  void SetHisto1DDraw(G4String str) {histo1DDraw = str;};
  void SetHisto1DSave(G4String str) {histo1DSave = str;};
  void SetHisto2DDraw(G4String str) {histo2DDraw = str;};
  void SetHisto2DSave(G4String str) {histo2DSave = str;};
  void SetHisto2DMode(G4String str) {histo2DMode = str;};

  G4String GetHisto2DMode() {return histo2DMode;};

private:
  G4VAnalysisSystem* analysisSystem;
  IPlotter* pl;
  IVectorFactory* fVectorFactory;
  IHistogramFactory* histoFactory;

  IHistogram1D* energy;
  IHistogram1D* hits;
  IHistogram2D* posXZ;
  IHistogram2D* posYZ;

  GammaRayTelDetectorConstruction*    GammaRayTelDetector;

  G4String histo1DDraw;
  G4String histo1DSave;
  G4String histo2DDraw;
  G4String histo2DSave;
  G4String histo2DMode;

  GammaRayTelAnalysisMessenger* analysisMessenger;
};


#endif
#endif









