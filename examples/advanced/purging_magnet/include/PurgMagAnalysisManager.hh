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
// Code developed by:
//  S.Larsson
//
//    ************************************
//    *                                  *
//    *    PurgMagAnalysisManager.hh     *
//    *                                  *
//    ************************************
//
// $Id: PurgMagAnalysisManager.hh 84477 2014-10-16 08:44:04Z gcosmo $
//

#ifndef G4PROCESSTESTANALYSIS_HH
#define G4PROCESSTESTANALYSIS_HH

#include "globals.hh"
#include <vector>
#include "G4ThreeVector.hh"

//!Uncomment #include to switch to ROOT or XML output file
//#include "g4root.hh"
//#include "g4xml.hh"

//! Default here is a CSV file
#include "g4csv.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class PurgMagAnalysisManager
{
public:  
  ~PurgMagAnalysisManager();
  void book();
  void finish();
  static PurgMagAnalysisManager* getInstance();
  void fill_Tuple_Electrons(G4double,G4double,G4double,G4double,G4double,
			    G4double,G4double);
  void fill_Tuple_Gamma(G4double,G4double,G4double,G4double,G4double,
			G4double,G4double);
  void fill_Tuple_Positrons(G4double,G4double,G4double,G4double,G4double,
			    G4double,G4double);
private:  
  static PurgMagAnalysisManager* instance;
  
private:
  PurgMagAnalysisManager();
};

#endif




