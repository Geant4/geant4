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
// $Id: G4ProcessTestAnalysis.hh,v 1.1 2001-11-06 01:38:24 pia Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Author: A. Pfeiffer (Andreas.Pfeiffer@cern.ch) 
//         (copy of his UserAnalyser class)
//
// History:
// -----------
//  5 Nov 2001   MGP        Created file according to A. Pfeiffer's class
//
// -------------------------------------------------------------------
// Class description:
// Analysis of the process test
// Further documentation available from http://www.ge.infn.it/geant4/lowE/index.html

// -------------------------------------------------------------------

#ifndef G4PROCESSTESTANALYSIS_HH
#define G4PROCESSTESTANALYSIS_HH

#include "globals.hh"
#include "g4std/map"

// Histogramming (from AIDA and Anaphe):
#include "Interfaces/IHistoManager.h"
#include "Interfaces/IHistogram1D.h"
#include "Interfaces/IHistogram2D.h"

// For NtupleTag from Anaphe
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
  Lizard::NTupleFactory* factory;

  // ---- NOTE ----
  // Histograms are compliant to AIDA interfaces, ntuples are Lizard specific

  G4std::map< G4String, IHistogram*, G4std::less<G4String> > histo1D;
  //G4std::map< G4String, IHistogram*, G4std::less<G4String> > histo2D;
  //G4std::map< G4String, IHistogram*, G4std::less<G4String> > histo3D;
  G4std::map< G4String, Lizard::NTuple*, G4std::less<G4String> > ntuples;

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
