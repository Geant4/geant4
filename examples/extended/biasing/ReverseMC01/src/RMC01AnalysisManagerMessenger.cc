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
/// \file biasing/ReverseMC01/src/RMC01AnalysisManagerMessenger.cc
/// \brief Implementation of the RMC01AnalysisManagerMessenger class
//
//
//////////////////////////////////////////////////////////////
//      Class Name:        RMC01AnalysisManagerMessenger
//        Author:               L. Desorgher
//         Organisation:         SpaceIT GmbH
//        Contract:        ESA contract 21435/08/NL/AT
//         Customer:             ESA/ESTEC
//////////////////////////////////////////////////////////////

#include "RMC01AnalysisManagerMessenger.hh"

#include "RMC01AnalysisManager.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RMC01AnalysisManagerMessenger::RMC01AnalysisManagerMessenger(
                                    RMC01AnalysisManager* analysisManager)
: G4UImessenger(),
  fAnalysisManager(analysisManager),
  fAnalysisDir(0),
  fSetPrecisionForConvergenceTestCmd(0),
  fSetExpSpectrumToNormaliseAdjResCmd(0),
  fSetPowerLawSpectrumToNormaliseAdjResCmd(0)
{ 
  fAnalysisDir = new G4UIdirectory("/RMC01/analysis/");
  fAnalysisDir->SetGuidance("Analysis commands");
 
  G4UIparameter* fluence_par = new G4UIparameter("Fluence",'d',true);
  fluence_par->SetParameterRange("Fluence > 0");
  fluence_par->SetGuidance("Omnidirectional fluence for primary spectrum");
  
  G4UIparameter* fluence_unit_par = new G4UIparameter("Fluence_unit",'s',true);
  fluence_unit_par->SetParameterCandidates("1/cm2 1/m2 cm-2 m-2");
  
  G4UIparameter* alpha_par = new G4UIparameter("alpha",'d',true);
  
  G4UIparameter* e0_par = new G4UIparameter("E0",'d',true);
  e0_par->SetParameterRange("E0 > 0");
  
  G4UIparameter* e1_par = new G4UIparameter("E1",'d',true);
  e1_par->SetParameterRange("E1 > 0");
  
  G4UIparameter* e2_par = new G4UIparameter("E2",'d',true);
  e2_par->SetParameterRange("E2 > 0");
  
  G4UIparameter* e_unit_par = new G4UIparameter("E_unit",'s',true);
  e_unit_par->SetParameterCandidates("eV keV MeV GeV TeV");
  
  G4UIparameter* part_name_par = new G4UIparameter("particle_name",'s',true);
  part_name_par->SetParameterCandidates("e- gamma proton ");
  
  fSetPowerLawSpectrumToNormaliseAdjResCmd =
   new G4UIcommand("/RMC01/analysis/SetPowerLawPrimSpectrumForAdjointSim",this);
  fSetPowerLawSpectrumToNormaliseAdjResCmd
     ->SetGuidance("Set the primary spectrum to which adjoint simulation "
                  "results will be normalised as a power law (Ekin^-alpha).");
  fSetPowerLawSpectrumToNormaliseAdjResCmd->SetParameter(part_name_par);
  fSetPowerLawSpectrumToNormaliseAdjResCmd->SetParameter(fluence_par);
  fSetPowerLawSpectrumToNormaliseAdjResCmd->SetParameter(fluence_unit_par);
  fSetPowerLawSpectrumToNormaliseAdjResCmd->SetParameter(alpha_par);
  fSetPowerLawSpectrumToNormaliseAdjResCmd->SetParameter(e1_par);
  fSetPowerLawSpectrumToNormaliseAdjResCmd->SetParameter(e2_par);
  fSetPowerLawSpectrumToNormaliseAdjResCmd->SetParameter(e_unit_par);
  fSetPowerLawSpectrumToNormaliseAdjResCmd
                             ->AvailableForStates(G4State_PreInit,G4State_Idle);
  

  fSetExpSpectrumToNormaliseAdjResCmd = new G4UIcommand("/RMC01/analysis/"
                                  "SetExponentialSpectrumForAdjointSim",this);
  fSetExpSpectrumToNormaliseAdjResCmd
    ->SetGuidance("Set the primary spectrum to which adjoint simulation results"
                          "will be normalised as exponential (exp(-Ekin/E0)).");
  fSetExpSpectrumToNormaliseAdjResCmd
                             ->SetParameter(new G4UIparameter(*part_name_par));
  fSetExpSpectrumToNormaliseAdjResCmd
                               ->SetParameter(new G4UIparameter(*fluence_par));
  fSetExpSpectrumToNormaliseAdjResCmd
                          ->SetParameter(new G4UIparameter(*fluence_unit_par));
  fSetExpSpectrumToNormaliseAdjResCmd
                                 ->SetParameter(new G4UIparameter(*e0_par));
  fSetExpSpectrumToNormaliseAdjResCmd
                                 ->SetParameter(new G4UIparameter(*e1_par));
  fSetExpSpectrumToNormaliseAdjResCmd
                                 ->SetParameter(new G4UIparameter(*e2_par));
  fSetExpSpectrumToNormaliseAdjResCmd
                                ->SetParameter(new G4UIparameter(*e_unit_par));
  fSetExpSpectrumToNormaliseAdjResCmd
                             ->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  
  fSetPrecisionForConvergenceTestCmd = new G4UIcmdWithADouble("/RMC01/analysis/"
                                         "SetExpectedPrecisionOfResults",this);
  fSetPrecisionForConvergenceTestCmd
    ->SetGuidance("Set the precision in % that the computed energy deposited "
          "in the sensitive volume should reached. If this precision is reached"
          " before the end of the run, the run is aborted and the results are "
          "registered.");
  fSetPrecisionForConvergenceTestCmd->SetParameterName("Precision",true); 
  fSetPrecisionForConvergenceTestCmd
                             ->AvailableForStates(G4State_PreInit,G4State_Idle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RMC01AnalysisManagerMessenger::~RMC01AnalysisManagerMessenger()
{
  delete fAnalysisDir;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RMC01AnalysisManagerMessenger::SetNewValue(
                                      G4UIcommand* command,G4String newValue)
{
  if( command == fSetPowerLawSpectrumToNormaliseAdjResCmd){
    G4double  alpha,e1,e2,fluence;
    G4String  f_unit,e_unit,part_name;
    const char* nv = (const char*)newValue;
    std::istringstream is(nv);
    is >> part_name>>fluence>>f_unit>>alpha>>e1>>e2>>e_unit;
    
    G4double factor_f_unit=1/cm2;
    if (f_unit == "1/m2" ||   f_unit =="m-2") factor_f_unit=1/m2;
    fluence*=factor_f_unit;
    e1*= G4UnitDefinition::GetValueOf(e_unit);
    e2*= G4UnitDefinition::GetValueOf(e_unit);
    fAnalysisManager->SetPrimaryPowerLawSpectrumForAdjointSim(
                                            part_name, fluence, alpha, e1, e2);
  }
  else if( command == fSetExpSpectrumToNormaliseAdjResCmd){
          G4double  e0,e1,e2,fluence;
    G4String  f_unit,e_unit,part_name;
    const char* nv = (const char*)newValue;
    std::istringstream is(nv);
    is >> part_name>>fluence>>f_unit>>e0>>e1>>e2>>e_unit;
    
    G4double factor_f_unit=1/cm2;
    if (f_unit == "1/m2" || f_unit =="m-2") factor_f_unit=1/m2;
    
    fluence*=factor_f_unit;
    e0*= G4UnitDefinition::GetValueOf(e_unit);
    e1*= G4UnitDefinition::GetValueOf(e_unit);
    e2*= G4UnitDefinition::GetValueOf(e_unit);
        
    fAnalysisManager->SetPrimaryExpSpectrumForAdjointSim(part_name,
                                                          fluence, e0, e1, e2);
          
  }
  else if( command == fSetPrecisionForConvergenceTestCmd){ 
          fAnalysisManager->SetPrecision(
             fSetPrecisionForConvergenceTestCmd->GetNewDoubleValue(newValue));
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
