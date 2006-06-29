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
