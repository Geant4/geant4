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
// $Id: RunAction.hh,v 1.2 2003/10/10 16:21:27 maire Exp $
// GEANT4 tag $Name: geant4-06-00 $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef RunAction_h
#define RunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"

class DetectorConstruction;
class PhysicsList;
class PrimaryGeneratorAction;
class G4Run;

#ifdef G4ANALYSIS_USE
namespace AIDA {
 class ITree;
 class IHistogram1D;
} 
#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class RunAction : public G4UserRunAction
{
  public:
    RunAction(DetectorConstruction*, PhysicsList*,PrimaryGeneratorAction*);
   ~RunAction();

  public:
    void BeginOfRunAction(const G4Run*);
    void   EndOfRunAction(const G4Run*);
    
    G4double* GetTallyEdep() {return tallyEdep;};
    
#ifdef G4ANALYSIS_USE
    AIDA::IHistogram1D* GetHisto(G4int id) {return histo[id];}
    
    G4double GetBinLength() {return binLength;};
    G4double GetOffsetX()   {return offsetX;} 
#endif
           
  private:  
    void bookHisto();
    void cleanHisto();
    
  private:
    DetectorConstruction*   detector;
    PhysicsList*            physics;
    PrimaryGeneratorAction* kinematic;
    G4double*               tallyEdep;   

#ifdef G4ANALYSIS_USE
  private:        
    AIDA::ITree* tree;
    AIDA::IHistogram1D* histo[1];
    
    G4double binLength;
    G4double offsetX;
#endif                     
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

