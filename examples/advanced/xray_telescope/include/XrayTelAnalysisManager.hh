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
// **********************************************************************
// *                                                                    *
// *                    GEANT 4 xray_telescope advanced example         *
// *                                                                    *
// * MODULE:            XrayTelAnalysisManager.cc                       *     
// * -------                                                            *
// *                                                                    *
// * Version:           0.2                                             *
// * Date:              30/11/00                                        *
// * Author:            A. Pfeiffer, G. Barrand, MG Pia, R Nartallo     *
// * Organisation:      ESA/ESTEC, Noordwijk, THe Netherlands           *
// *                                                                    *
// **********************************************************************
//
// CHANGE HISTORY
// --------------
//
// 30.11.2000 M.G. Pia, R. Nartallo
// - Simplification of code
// - Inheritance directly from the base class G4VAnalysisManager instead
//   of the derived class G4AnalysisManager
//
// 16.10.2000 G. Barrand
// - First implementation of XrayAnalysisManager class
// - Provision of code for various AIDA and non-AIDA systems
//
// **********************************************************************

#ifndef XrayTelAnalysisManager_h
#define XrayTelAnalysisManager_h 1

#include "G4VAnalysisManager.hh"

#include "globals.hh"
#include "g4std/vector"
#include "G4ThreeVector.hh"                                                                         
class G4SteppingManager;
class G4VAnalysisSystem;

#ifdef G4ANALYSIS_USE
class IHistogramFactory;
class IHistogram1D;
class IHistogram2D;
#endif

#ifdef G4ANALYSIS_USE
class XrayTelAnalysisManager: public G4VAnalysisManager {
#endif
#ifndef G4ANALYSIS_USE
class XrayTelAnalysisManager {
#endif
public:
  XrayTelAnalysisManager(const G4String&);
  ~XrayTelAnalysisManager();

  G4bool RegisterAnalysisSystem(G4VAnalysisSystem*);

#ifdef G4ANALYSIS_USE
  IHistogramFactory* GetHistogramFactory(const G4String&);

  virtual void Store(IHistogram* = 0,const G4String& = "");
  virtual void Plot(IHistogram*);
#endif

  void BeginOfRun(); 
  void EndOfRun(); 
  void Step(const G4SteppingManager*);

private:
  G4VAnalysisSystem* analysisSystem;

#ifdef G4ANALYSIS_USE
  IHistogramFactory* hFactory;
  IHistogram1D* enteringEnergyHistogram;
  IHistogram2D* yzHistogram;
#endif
};

#endif





