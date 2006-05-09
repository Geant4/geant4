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
// $Id: RunAction.hh,v 1.1 2006-05-09 14:03:03 maire Exp $
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
	 
    void SumPsiPlane (G4double psi)  
         {psiPrj += psi; psiPrj2 += psi*psi;}
	 
    void SumTetaPlane (G4double teta)  
         {tetPrj += teta; tetPrj2 += teta*teta;}
	 
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
    G4double psiPrj, psiPrj2;
    G4double tetPrj, tetPrj2;     
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

