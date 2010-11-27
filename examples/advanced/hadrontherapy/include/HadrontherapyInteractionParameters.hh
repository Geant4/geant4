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
// This is the *BASIC* version of Hadrontherapy, a Geant4-based application
// See more at: http://g4advancedexamples.lngs.infn.it/Examples/hadrontherapy
//
// Visit the Hadrontherapy web site (http://www.lns.infn.it/link/Hadrontherapy) to request 
// the *COMPLETE* version of this program, together with its documentation;
// Hadrontherapy (both basic and full version) are supported by the Italian INFN
// Institute in the framework of the MC-INFN Group
//

#ifndef HadrontherapyInteractionParameters_H
#define HadrontherapyInteractionParameters_H 1

#include "G4EmCalculator.hh"
#include "G4NistMaterialBuilder.hh"
#include "G4NistElementBuilder.hh"

#ifdef G4ANALYSIS_USE_ROOT 
#include "TROOT.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TH1F.h" 
#include "TH2F.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TStyle.h"
#endif

class HadrontherapyDetectorConstruction;
class HadrontherapyParameterMessenger; 
class G4ParticleDefinition;
class G4Material;

class HadrontherapyInteractionParameters : public G4EmCalculator 
{
public:

  HadrontherapyInteractionParameters(G4bool);
  ~HadrontherapyInteractionParameters();

  // Get data for Mass SP    
  // G4NistMaterialBuilder class materials
  // User must provide: material kinetic energy lower limit, kinetic energy upper limit, number of points to retrieve,
  // [particle], [output filename].

  G4bool GetStoppingTable (const G4String& vararg);
  G4double GetStopping (G4double energy,
			const G4ParticleDefinition*, 
			const G4Material*, 
			G4double density = 0.);
#ifdef G4ANALYSIS_USE_ROOT 
  void PlotStopping(const G4String&);
#endif
  void ListOfNistMaterials (const G4String& vararg);
  void BeamOn();
  bool ParseArg (const G4String& vararg);	

private:
  G4Material* GetNistMaterial(G4String material);
  G4NistElementBuilder* nistEle;
  G4NistMaterialBuilder* nistMat;
  std::ofstream outfile;
  std::ostream data;
  G4Material* Pmaterial;
  HadrontherapyParameterMessenger* pMessenger; 
  bool beamFlag;

#ifdef G4ANALYSIS_USE_ROOT 
  TCanvas *theRootCanvas;
  TGraph *theRootGraph;
  TAxis *axisX, *axisY;
#endif
  G4double kinEmin, kinEmax, npoints;
  G4String particle, material, filename; 
  G4double dedxtot, density;
  std::vector<G4double> energy;
  std::vector<G4double> massDedx;

};
#endif

