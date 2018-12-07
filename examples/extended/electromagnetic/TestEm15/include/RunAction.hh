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
/// \file electromagnetic/TestEm15/include/RunAction.hh
/// \brief Definition of the RunAction class
//
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
    RunAction(DetectorConstruction*, PrimaryGeneratorAction*);
   ~RunAction();

  public:
    virtual void BeginOfRunAction(const G4Run*);
    virtual void   EndOfRunAction(const G4Run*);

    void CountProcesses(G4String);
    void SumPathLength (G4double truepl, G4double geompl) 
         {fTotalCount++; 
          fTruePL += truepl; fTruePL2 += truepl*truepl;
          fGeomPL += geompl; fGeomPL2 += geompl*geompl;
         };
         
    void SumLateralDisplacement (G4double displa)  
         {fLDispl += displa; fLDispl2 += displa*displa;}
         
    void SumPsi (G4double psi)  
         {fPsiSpa += psi; fPsiSpa2 += psi*psi;}
         
    void SumTetaPlane (G4double teta)  
         {fTetPrj += teta; fTetPrj2 += teta*teta;}
                  
    void SumPhiCorrel (G4double correl)  
         {fPhiCor += correl; fPhiCor2 += correl*correl;}
                  
   G4double ComputeMscHighland(G4double pathLength);
                                     
  private:
    DetectorConstruction*   fDetector;
    PrimaryGeneratorAction* fPrimary;
    ProcessesCount*         fProcCounter;
    HistoManager*           fHistoManager;
        
    G4int    fTotalCount;
    G4double fTruePL, fTruePL2;
    G4double fGeomPL, fGeomPL2;
    G4double fLDispl, fLDispl2;
    G4double fPsiSpa, fPsiSpa2;
    G4double fTetPrj, fTetPrj2;
    G4double fPhiCor, fPhiCor2;     
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

