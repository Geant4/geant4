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
// $Id: RunAction.hh,v 1.1 2007-08-16 10:32:04 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef RunAction_h
#define RunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"
#include "QGSP.hh"

class DetectorConstruction;
class RunActionMessenger;
class PrimaryGeneratorAction;
class G4Run;

namespace AIDA {
 class IAnalysisFactory;
 class ITree;
 class IHistogram1D;
} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class RunAction : public G4UserRunAction
{
public:
  RunAction(DetectorConstruction* ,PrimaryGeneratorAction*);
  virtual ~RunAction();

  void BeginOfRunAction(const G4Run*);
  void   EndOfRunAction(const G4Run*);
           
  G4double GetBinLength() {return binLength;};
  G4double GetOffsetX()   {return offsetX;} 
  void     SetBinSize(G4double size);
  void     FillHisto(G4int id, G4double x, G4double weight = 1.0);

  void     SetVerbose(G4int verbose) {verboseLevel = verbose;}
	void     SetHistoName(G4String name) {fname = name;}
  void     SetHistoType(G4String type) {ftype = type;}
  G4int    GetVerbose() {return verboseLevel;}
    
  void AddProjRange (G4double x) {projRange += x; projRange2 += x*x;};
                   
private:  
  void bookHisto();
  void saveHisto();
    
  DetectorConstruction*   detector;
  PrimaryGeneratorAction* kinematic;

  RunActionMessenger*	  runActionMessenger;

  G4int                   verboseLevel;

  G4double                binLength;
  G4double                offsetX;
  G4double                projRange, projRange2;
	G4String 								ftype, fname;

             
  AIDA::IAnalysisFactory* af;  
  AIDA::ITree*            tree;
  AIDA::IHistogram1D*     histo[5];        
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

