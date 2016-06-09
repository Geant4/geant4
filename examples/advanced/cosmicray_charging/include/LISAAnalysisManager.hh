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



