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
// $Id: XrayTelAnalysis.hh,v 1.1 2001-12-05 18:56:43 nartallo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Author: A. Pfeiffer (Andreas.Pfeiffer@cern.ch) 
//         (copy of his UserAnalyser class)
//
// History:
// -----------
//  7 Nov 2001   MGP  Implemented according to A. Pfeiffer's instructions
//
// -------------------------------------------------------------------
// Class description:
// Example of analysis in a simulation application (histograms, ntuples etc.)
// This class follows the singleton design pattern; 
// it is responsible for the analysis management and algorithms 
// Histograms are compliant with AIDA 1.0 (the header files are
// copies from the Anaphe Interfaces directory), except for the usage of
// IHistoManager which is Anaphe/Lizard-specific
// For ntuples, for which the AIDA interfaces and compliant implementations
// are still in progress at the time of the current Geant4 release, 
// an implementation with Anaphe/Lizard is shown
// Other implementations specific to an analysis system are possible too
// (see, for instance, JAS and OpenScientist documentation from the links
// in http://aida.freehep.org/)
// The implementation of the usage of ntuples shown in this example
// is expected to change in future Geant4 releases, when AIDA interfaces 
// and related implementations would be available
// Further documentation is available from: http://www.ge.infn.it/geant4/lowE/
//                                          http://aida.freehep.org/
//                                          http://cern.ch/anaphe/

// -------------------------------------------------------------------

#ifndef G4PROCESSTESTANALYSIS_HH
#define G4PROCESSTESTANALYSIS_HH

#include "globals.hh"

// Histogramming from AIDA (through their copy in Anaphe)
#include "Interfaces/IHistogram1D.h"
#include "Interfaces/IHistogram2D.h"

// Histogramming from Anaphe
#include "Interfaces/IHistoManager.h"

// Ntuples from Anaphe
#include "NtupleTag/LizardNTupleFactory.h"
#include "NtupleTag/LizardQuantity.h"

// Vectors from Anaphe
#include "Interfaces/IVector.h"
#include "Interfaces/IVectorFactory.h"

// Plotting from Anaphe
#include "Interfaces/IPlotter.h"

class G4Track;
class NTuple;

class XrayTelAnalysis
{
public:

  ~XrayTelAnalysis();

  void book();
  
  void finish();
  
  void analyseStepping(const G4Track& track, G4bool entering);

  static XrayTelAnalysis* getInstance();

private:

  XrayTelAnalysis();

  void plot1D(IHistogram1D* histo);
  void plot2D(IHistogram2D* histo);

  static XrayTelAnalysis* instance;

  IHistoManager* histoManager;
  IVectorFactory* vectorFactory;
  IPlotter* plotter;

  // ---- NOTE ----
  // Histograms are compliant to AIDA interfaces, ntuples are Lizard specific
 
  Lizard::NTuple* ntuple;
  Lizard::NTupleFactory* ntFactory;

  // Quantities for the ntuple
  Lizard::Quantity<float> eKin;
  Lizard::Quantity<float> y;
  Lizard::Quantity<float> z;
  Lizard::Quantity<float> dirX;
  Lizard::Quantity<float> dirY;
  Lizard::Quantity<float> dirZ;

};

#endif 
