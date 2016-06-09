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
// $Id: RunAction.hh,v 1.3 2006-06-29 16:46:36 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef RunAction_h
#define RunAction_h 1

#include "G4UserRunAction.hh"
#include "ProcessesCount.hh"
#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class DetectorConstruction;
class PrimaryGeneratorAction;
class HistoManager;
class G4Run;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class RunAction : public G4UserRunAction
{
  public:
    RunAction(DetectorConstruction*, PrimaryGeneratorAction*, HistoManager*);
   ~RunAction();

  public:
    void BeginOfRunAction(const G4Run*);
    void   EndOfRunAction(const G4Run*);

    void CountProcesses(G4String);
    void SumPathLength (G4double truepl, G4double geompl) 
         {totalCount++; 
	  truePL += truepl; truePL2 += truepl*truepl;
	  geomPL += geompl; geomPL2 += geompl*geompl;
	 };
	 
    void SumLateralDisplacement (G4double displa)  
         {lDispl += displa; lDispl2 += displa*displa;}
	 
    void SumPsi (G4double psi)  
         {psiSpa += psi; psiSpa2 += psi*psi;}
	 
    void SumTetaPlane (G4double teta)  
         {tetPrj += teta; tetPrj2 += teta*teta;}
	 	 
    void SumPhiCorrel (G4double correl)  
         {phiCor += correl; phiCor2 += correl*correl;}
	 	 
   G4double ComputeMscHighland(G4double pathLength);
	 	 	           
  private:
    DetectorConstruction*   detector;
    PrimaryGeneratorAction* primary;
    ProcessesCount*         ProcCounter;
    HistoManager*           histoManager;
        
    G4int totalCount;
    G4double truePL, truePL2;
    G4double geomPL, geomPL2;
    G4double lDispl, lDispl2;
    G4double psiSpa, psiSpa2;
    G4double tetPrj, tetPrj2;
    G4double phiCor, phiCor2;     
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

