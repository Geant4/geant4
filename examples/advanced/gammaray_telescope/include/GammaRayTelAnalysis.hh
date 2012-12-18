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
#ifdef G4ANALYSIS_USE
//
// $Id: GammaRayTelAnalysis.hh,v 1.18 2006-06-29 15:54:51 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// ------------------------------------------------------------
//      GEANT 4 class header file
//      CERN Geneva Switzerland
//     
//
//      ------------ GammaRayTelAnalysis  ------
//           by R.Giannitrapani, F. Longo & G.Santin (30 nov 2000)
//
// 07.12.2001 A.Pfeiffer
// - integrated Guy's addition of the ntuple
//
// 06.12.2001 A.Pfeiffer
// - updating to new design (singleton)
//
// 22.11.2001 G.Barrand
// - Adaptation to AIDA
//
// ************************************************************

#ifndef GammaRayTelAnalysis_h 
#define GammaRayTelAnalysis_h 1

#include "globals.hh"
#include <vector>
#include "G4ThreeVector.hh"

#include <AIDA/AIDA.h>

using namespace AIDA;

class GammaRayTelAnalysisMessenger;
class GammaRayTelDetectorConstruction;

class AIDA::IAnalysisFactory;
class AIDA::ITree;
class AIDA::IHistogramFactory;
class AIDA::ITupleFactory;
//class AIDA::IPlotter;
class AIDA::IHistogram1D;
class AIDA::IHistogram2D;

class GammaRayTelAnalysis {
public:
  virtual ~GammaRayTelAnalysis();
  
public:


  //  void BeginOfRun(G4int n);

  void BeginOfRun();
  void EndOfRun(G4int n);
  void EndOfEvent(G4int flag);

  void Init();
  void Finish();

  void SetHisto1DDraw(G4String str) {histo1DDraw = str;};
  void SetHisto1DSave(G4String str) {histo1DSave = str;};
  void SetHisto2DDraw(G4String str) {histo2DDraw = str;};
  void SetHisto2DSave(G4String str) {histo2DSave = str;};
  void SetHisto2DMode(G4String str) {histo2DMode = str;};

  G4String GetHisto2DMode() {return histo2DMode;};

  void InsertPositionXZ(double x, double z);
  void InsertPositionYZ(double y, double z);
  void InsertEnergy(double en);
  void InsertHits(int nplane);

  void setNtuple(float E, float p, float x, float y, float z);

  static GammaRayTelAnalysis* getInstance();

private:

  GammaRayTelAnalysis();

  void plot1D(IHistogram1D* histo);
  void plot2D(IHistogram2D* histo);
  void Plot();

private:
  static GammaRayTelAnalysis* instance;

  GammaRayTelDetectorConstruction*    GammaRayTelDetector;

  IAnalysisFactory* analysisFactory;
  ITree* tree;
  //IPlotter* plotter;
  ITuple* tuple;

  IHistogram1D* energy;
  IHistogram1D* hits;
  IHistogram2D* posXZ;
  IHistogram2D* posYZ;

  G4String histo1DDraw;
  G4String histo1DSave;
  G4String histo2DDraw;
  G4String histo2DSave;
  G4String histo2DMode;

  GammaRayTelAnalysisMessenger* analysisMessenger;
};


#endif
#endif
