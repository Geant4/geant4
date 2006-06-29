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
//
// $Id: G4ProcessTestAnalysis.hh,v 1.7 2006-06-29 19:48:28 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Author: A. Pfeiffer (Andreas.Pfeiffer@cern.ch) 
//         (copy of his UserAnalyser class)
//
// History:
// -----------
//  5 Nov 2001   MGP        Implemented according to A. Pfeiffer's design 
//                          and suggestions
//
// -------------------------------------------------------------------
// Class description:
// Analysis of the process test
// Further documentation available from http://www.ge.infn.it/geant4/lowE/index.html

// -------------------------------------------------------------------

#ifndef G4PROCESSTESTANALYSIS_HH
#define G4PROCESSTESTANALYSIS_HH

#include "globals.hh"
// #include "g4std/map"

// Histogramming (from AIDA)
#include "Interfaces/IHistogram1D.h"
#include "Interfaces/IHistogram2D.h"

// Histogramming (from Anaphe)
#include "Interfaces/IHistoManager.h"

// Ntuples (from Anaphe)
#include "NtupleTag/LizardNTupleFactory.h"

class G4ParticleChange;
class G4Track;

class G4ProcessTestAnalysis
{
public:

  ~G4ProcessTestAnalysis();

  void book(const G4String& storeName = "analysis");
  
  void finish();
  
  void analyseGeneral(const G4Track& track,
		      const G4ParticleChange* particleChange);

  void analyseSecondaries(const G4ParticleChange* particleChange);

  static G4ProcessTestAnalysis* getInstance();

private:

  G4ProcessTestAnalysis();

  static G4ProcessTestAnalysis* instance;
  IHistoManager* histoManager;
  Lizard::NTupleFactory* ntFactory;
  Lizard::NTuple* ntuple1;
  Lizard::NTuple* ntuple2;

  // ---- NOTE ----
  // Histograms are compliant to AIDA interfaces, ntuples are Lizard specific

  // std::map< G4String, IHistogram1D*, std::less<G4String> > histo1D;
  // std::map< G4String, IHistogram2D*, std::less<G4String> > histo2D;
  // std::map< G4String, IHistogram3D*, std::less<G4String> > histo3D;
  // std::map< G4String, Lizard::NTuple*, std::less<G4String> > ntuples;

  // Quantities for the general ntuple
  Lizard::Quantity<float> initialEnergy;
  Lizard::Quantity<float> energyChange;
  Lizard::Quantity<float> pxChange;
  Lizard::Quantity<float> pyChange;
  Lizard::Quantity<float> pzChange;
  Lizard::Quantity<float> eDeposit;
  Lizard::Quantity<float> thetaChange;
  Lizard::Quantity<float> phiChange;
  Lizard::Quantity<float> nElectrons;
  Lizard::Quantity<float> nPositrons;
  Lizard::Quantity<float> nPhotons;

  // Quantities for the secondary ntuple 
  //  Lizard::Quantity<float> initialEnergy;
  Lizard::Quantity<float> px;
  Lizard::Quantity<float> py;
  Lizard::Quantity<float> pz;
  Lizard::Quantity<float> p;
  Lizard::Quantity<float> e;
  Lizard::Quantity<float> eKin;
  Lizard::Quantity<float> theta;
  Lizard::Quantity<float> phi;
  Lizard::Quantity<float> partType;

};

#endif 
