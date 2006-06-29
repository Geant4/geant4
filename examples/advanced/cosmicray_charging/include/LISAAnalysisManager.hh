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

// ********************************************************************
// *                                                                  *
// * cosmicray_charging advanced example for Geant4                   *
// * (adapted simulation of test-mass charging in the LISA mission)   *
// *                                                                  *
// * Henrique Araujo (h.araujo@imperial.ac.uk) & Peter Wass           *
// * Imperial College London                                          *
// *                                                                  *
// ********************************************************************


#ifdef G4ANALYSIS_USE
#ifndef LISAAnalysisManager_h
#define LISAAnalysisManager_h 1

#include "globals.hh"

#include <AIDA/AIDA.h>

namespace AIDA {
  class IAnalysisFactory;
  class ITree;
  class IHistogramFactory;
  class ITupleFactory;
  class ITuple;
  class IHistogram1D;
  class IHistogram2D;
  class IPlotter;
  class IFitter;
  class IFitResult;
  class IFitData;
  class IRangeSet;
  class IFitParameterSettings;
  class IFunctionFactory;
  class IFunction;
  class IFitFactory;
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class LISAAnalysisManager {

  public:
    virtual ~LISAAnalysisManager();
    void Init();
    void Dispose();


    // Run analysis
    void bookRun(G4String);
    void FinishRun();
    void analyseRun
     (G4int evt, G4int tm, G4double energy, G4int charge, 
      G4int in, G4int out, long seed1, long seed2);

    // grab instance
    static LISAAnalysisManager* getInstance();
 

  private:
  
    // private constructor for singleton
    LISAAnalysisManager();
 
    static LISAAnalysisManager* instance;
  
 
    AIDA::IAnalysisFactory  *af;
    AIDA::ITreeFactory      *tf; 
    AIDA::ITree             *run_tree;
    AIDA::ITupleFactory     *run_tpf;
    AIDA::ITuple            *run_tuple;


};
#endif
#endif



