#ifndef RunAction_h
#define RunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"
#include "G4ios.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class SteppingAction;
class DetectorConstruction;
class PrimaryGeneratorAction;

class G4Run;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class RunAction : public G4UserRunAction
{
  public:

    RunAction(DetectorConstruction*, PrimaryGeneratorAction*);
   ~RunAction();

    void BeginOfRunAction(const G4Run*);
    void EndOfRunAction(const G4Run*);

    void SumEnergy1(G4int k, G4double de){energyDeposit1[k] += de; };  	
    void SumEnergy2(G4int k, G4double de){energyDeposit2[k] += de; };  	
    void SumEnergy3(G4int k, G4double de){energyDeposit3[k] += de; };  	

    void AddEnergy1(G4double de) {energyDepositRun1 += de;};
    void AddEnergy2(G4double de) {energyDepositRun2 += de;};
    void AddEnergy3(G4double de) {energyDepositRun3 += de;};
    


  private:
    DetectorConstruction*   detector;
    PrimaryGeneratorAction* primary;
    G4String matName1     ;
    G4String matName2     ;
    G4String matName3     ;

    G4double              energyDeposit1[110];
    G4double              energyDeposit2[110];
    G4double              energyDeposit3[110];
    G4double              normalizedvalue1[110];
    G4double              normalizedvalue2[110];
    G4double              normalizedvalue3[110];
    G4double MFP1;        
    G4double MFP2;        
    G4double MFP3;        
    G4double density1;        
    G4double density2;        
    G4double density3;        
    G4double energyDepositRun1; 
    G4double energyDepositRun2;
    G4double energyDepositRun3;
    G4String asciiFileName;
     

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

