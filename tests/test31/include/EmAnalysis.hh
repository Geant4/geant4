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
#ifndef EmAnalysis_h
#define EmAnalysis_h 1

//---------------------------------------------------------------------------
//
// ClassName:   EmAnalysis
//
// Description: 
//
// Author:      V.Ivanchenko 20/08/04
//
//----------------------------------------------------------------------------
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "globals.hh"
#include <vector>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class Histo;
class EmAnalysisMessenger;
class G4EmCalculator;

class EmAnalysis
{

public:
  // With description
  EmAnalysis();
  
  ~EmAnalysis();

  void clear();

  G4int AddHistOnCrossSection(const G4String& part_name, 
			      const G4String& mat_name,
			      const G4String& proc_name, 
			      const G4String& proc_type = "dedx", 
			      const G4String& h_id="");
  
  void setHisto1D(G4int id, G4int nbins=90, G4double emin=0.001*MeV, 
                  G4double emax=1.0*TeV);

  void setHisto1D(const G4String& h_id, G4int nbins=90, G4double emin=0.001*MeV, 
                  G4double emax=1.0*TeV);

  void activate1D(G4int id, G4bool val=true);

  void activate1D(const G4String& h_id, G4bool val=true);

  void setVerbose(G4int val);  

  void PrintHist(G4int val);  

  void saveToFile();  
 
private: // Without description
    
  EmAnalysis(const EmAnalysis&);
  const EmAnalysis& operator=(const EmAnalysis& right);

  // MEMBERS

  Histo* histo;
  EmAnalysisMessenger* messenger;
  G4EmCalculator*      calculator;
  G4int nHisto;
  G4int verbose; 
  G4int nbins;
  G4double cut;
  G4double xmin;
  G4double xmax;

  std::vector<G4String> particle;
  std::vector<G4String> material;
  std::vector<G4String> process;
  std::vector<G4String> idstr;
  std::vector<G4String> atype;
  std::vector<G4int>    histid;

};

#endif
