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
// 07.12.2001 A.Pfeiffer
// - merged with Guy Barrand's AIDA 2.2 port
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

#ifndef XrayTelAnalysis_h
#define XrayTelAnalysis_h 1

#include "globals.hh"
#include "g4std/vector"
#include "G4ThreeVector.hh"                                                                         

class G4Track;

class IAnalysisFactory;
class ITree;
class IHistogram1D;
class IHistogram2D;
class IPlotter;

class XrayTelAnalysis {
public:

  virtual ~XrayTelAnalysis();

  void book(); 
  void finish(); 
  void analyseStepping(const G4Track& track, G4bool entering);

  static XrayTelAnalysis* getInstance(int = 0, char** = 0);

private:

  XrayTelAnalysis(int,char**);

  static XrayTelAnalysis* instance;

  IAnalysisFactory* analysisFactory;
  ITree* tree;
  IHistogram1D* enteringEnergyHistogram;
  IHistogram2D* yzHistogram;
  IPlotter* plotter;
};

#endif

