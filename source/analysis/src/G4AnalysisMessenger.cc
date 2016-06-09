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
// $Id$

// Author: Ivana Hrivnacova, 04/07/2012  (ivana@ipno.in2p3.fr)
//
// This messenger class is a generalization of the HistoMessenger class,
// originally developed for the extended/electromagnetic examples
// by Michel Maire (michel.maire@lapp.in2p3.fr)

#include "G4AnalysisMessenger.hh"
#include "G4VAnalysisManager.hh"
#include "G4HnInformation.hh"

#include "G4UIdirectory.hh"
#include "G4UIcommand.hh"
#include "G4UIparameter.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithABool.hh"

#include <iostream>

//_____________________________________________________________________________
G4AnalysisMessenger::G4AnalysisMessenger(G4VAnalysisManager* manager)
  : G4UImessenger(),
    fManager(manager),
    fAnalysisDir(0),  
    fSetFileNameCmd(0),
    fSetHistoDirNameCmd(0),
    fSetNtupleDirNameCmd(0),
    fSetActivationCmd(0),
    fVerboseCmd(0),
    fH1Dir(0),  
    fCreateH1Cmd(0),
    fSetH1Cmd(0),
    fSetH1AsciiCmd(0), 
    fSetH1TitleCmd(0), 
    fSetH1XAxisCmd(0), 
    fSetH1YAxisCmd(0), 
    fSetH1ActivationCmd(0),
    fSetH1ActivationAllCmd(0),
    fH2Dir(0),  
    fCreateH2Cmd(0),
    fSetH2Cmd(0),
    fSetH2AsciiCmd(0), 
    fSetH2TitleCmd(0), 
    fSetH2XAxisCmd(0), 
    fSetH2YAxisCmd(0), 
    fSetH2ActivationCmd(0),
    fSetH2ActivationAllCmd(0)
{  
  fAnalysisDir = new G4UIdirectory("/analysis/");
  fAnalysisDir->SetGuidance("analysis control");

  fSetFileNameCmd = new G4UIcmdWithAString("/analysis/setFileName",this);
  fSetFileNameCmd->SetGuidance("Set name for the histograms & ntuple file");
  fSetFileNameCmd->SetParameterName("Filename", false);
  fSetFileNameCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
  
  fSetHistoDirNameCmd = new G4UIcmdWithAString("/analysis/setHistoDirName",this);
  fSetHistoDirNameCmd->SetGuidance("Set name for the histograms directory");
  fSetHistoDirNameCmd->SetParameterName("HistoDirName", false);
  fSetHistoDirNameCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
  
  fSetNtupleDirNameCmd = new G4UIcmdWithAString("/analysis/setNtupleDirName",this);
  fSetNtupleDirNameCmd->SetGuidance("Set name for the ntuple directory");
  fSetNtupleDirNameCmd->SetParameterName("NtupleDirName", false);
  fSetNtupleDirNameCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
  
  fSetActivationCmd = new G4UIcmdWithABool("/analysis/setActivation",this);
  G4String guidance = "Set activation. \n";
  guidance += "When this option is enabled, only the histograms marked as activated\n";
  guidance += "are returned, filled or saved on file.\n";
  guidance += "No warning is issued when Get or Fill is called on inactive histogram.";
  fSetActivationCmd->SetGuidance(guidance);
  fSetActivationCmd->SetParameterName("Activation",false);

  fVerboseCmd = new G4UIcmdWithAnInteger("/analysis/verbose",this);
  fVerboseCmd->SetGuidance("Set verbose level");
  fVerboseCmd->SetParameterName("VerboseLevel",false);
  fVerboseCmd->SetRange("VerboseLevel>=0 && VerboseLevel<=4");

  fH1Dir = new G4UIdirectory("/analysis/h1/");
  fH1Dir->SetGuidance("1D histograms control");

  CreateH1Cmd();
  SetH1Cmd();
  
  fSetH1AsciiCmd = new G4UIcmdWithAnInteger("/analysis/h1/setAscii",this);
  fSetH1AsciiCmd->SetGuidance("Print 1D histogram of #Id on ascii file.");
  fSetH1AsciiCmd->SetParameterName("Id",false);
  fSetH1AsciiCmd->SetRange("Id>=0");
  fSetH1AsciiCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
  
  SetH1TitleCmd();
  SetH1XAxisCmd();
  SetH1YAxisCmd();
  SetH1ActivationCmd();

  fSetH1ActivationAllCmd = new G4UIcmdWithABool("/analysis/h1/setActivationToAll",this);
  fSetH1ActivationAllCmd->SetGuidance("Set activation to all 1D histograms.");
  fSetH1ActivationAllCmd->SetParameterName("Activation",false);

  fH2Dir = new G4UIdirectory("/analysis/h2/");
  fH2Dir->SetGuidance("2D histograms control");

  CreateH2Cmd();
  SetH2Cmd();
  
  fSetH2AsciiCmd = new G4UIcmdWithAnInteger("/analysis/h2/setAscii",this);
  fSetH2AsciiCmd->SetGuidance("Print 2D histogram of #Id on ascii file.");
  fSetH2AsciiCmd->SetParameterName("Id",false);
  fSetH2AsciiCmd->SetRange("Id>=0");
  fSetH2AsciiCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
  
  SetH2TitleCmd();
  SetH2XAxisCmd();
  SetH2YAxisCmd();
  SetH2ZAxisCmd();
  SetH2ActivationCmd();

  fSetH2ActivationAllCmd = new G4UIcmdWithABool("/analysis/h2/setActivationToAll",this);
  fSetH2ActivationAllCmd->SetGuidance("Set activation to all 2D histograms.");
  fSetH2ActivationAllCmd->SetParameterName("Activation",false);

}

//_____________________________________________________________________________
G4AnalysisMessenger::~G4AnalysisMessenger()
{
  delete fSetFileNameCmd;
  delete fSetHistoDirNameCmd;
  delete fSetNtupleDirNameCmd;
  delete fSetActivationCmd;
  delete fVerboseCmd;
  delete fCreateH1Cmd;
  delete fSetH1Cmd;
  delete fSetH1AsciiCmd;  
  delete fSetH1TitleCmd;  
  delete fSetH1XAxisCmd;  
  delete fSetH1YAxisCmd;  
  delete fSetH1ActivationCmd;
  delete fSetH1ActivationAllCmd;
  delete fH1Dir;
  delete fCreateH2Cmd;
  delete fSetH2Cmd;
  delete fSetH2AsciiCmd;  
  delete fSetH2TitleCmd;  
  delete fSetH2XAxisCmd;  
  delete fSetH2YAxisCmd;  
  delete fSetH2ZAxisCmd;  
  delete fSetH2ActivationCmd;
  delete fSetH2ActivationAllCmd;
  delete fH2Dir;
  delete fAnalysisDir;
}

//
// private functions
//

//_____________________________________________________________________________
void G4AnalysisMessenger::CreateH1Cmd()
{
  G4UIparameter* h1Name = new G4UIparameter("name", 's', false);
  h1Name->SetGuidance("Histogram name (label)");
  
  G4UIparameter* h1Title = new G4UIparameter("title", 's', false);
  h1Title->SetGuidance("Histogram title");

  G4UIparameter* h1Nbins0 = new G4UIparameter("nbins0", 'i', true);
  h1Nbins0->SetGuidance("Number of bins (default = 100)");
  h1Nbins0->SetGuidance("Can be reset with /analysis/h1/set command");
  h1Nbins0->SetDefaultValue(100);
  
  G4UIparameter* h1ValMin0 = new G4UIparameter("valMin0", 'd', true);
  h1ValMin0->SetGuidance("Minimum value, expressed in unit (default = 0.)");
  h1ValMin0->SetGuidance("Can be reset with /analysis/h1/set command");
  h1ValMin0->SetDefaultValue(0.);
  
  G4UIparameter* h1ValMax0 = new G4UIparameter("valMax0", 'd', true);
  h1ValMax0->SetGuidance("Maximum value, expressed in unit (default = 1.)");
  h1ValMax0->SetGuidance("Can be reset with /analysis/h1/set command");
  h1ValMax0->SetDefaultValue(1.);

  G4UIparameter* h1ValUnit0 = new G4UIparameter("valUnit0", 's', true);
  h1ValUnit0->SetGuidance("The unit of valMin0 and valMax0");
  h1ValUnit0->SetDefaultValue("none");
  
  G4UIparameter* h1ValFcn0 = new G4UIparameter("valFcn0", 's', true);
  G4String fcnGuidance = "The function applied to filled values (log, log10, exp).\n";
  fcnGuidance += "Note that the unit parameter cannot be omitted in this case,\n";
  fcnGuidance += "but none value should be used insted.";
  h1ValFcn0->SetGuidance(fcnGuidance);
  h1ValFcn0->SetParameterCandidates("log log10 exp none");
  h1ValFcn0->SetDefaultValue("none");
  
  fCreateH1Cmd = new G4UIcommand("/analysis/h1/create", this);
  fCreateH1Cmd->SetGuidance("Create 1D histogram");
  fCreateH1Cmd->SetParameter(h1Name);
  fCreateH1Cmd->SetParameter(h1Title);
  fCreateH1Cmd->SetParameter(h1Nbins0);
  fCreateH1Cmd->SetParameter(h1ValMin0);
  fCreateH1Cmd->SetParameter(h1ValMax0);
  fCreateH1Cmd->SetParameter(h1ValUnit0);
  fCreateH1Cmd->SetParameter(h1ValFcn0);
  fCreateH1Cmd->AvailableForStates(G4State_PreInit, G4State_Idle);
}  


//_____________________________________________________________________________
void G4AnalysisMessenger::SetH1Cmd()
{
  G4UIparameter* h1Id = new G4UIparameter("id", 'i', false);
  h1Id->SetGuidance("Histogram id");
  h1Id->SetParameterRange("id>=0");
  
  G4UIparameter* h1Nbins = new G4UIparameter("nbins", 'i', false);
  h1Nbins->SetGuidance("Number of bins");
  
  G4UIparameter* h1ValMin = new G4UIparameter("valMin", 'd', false);
  h1ValMin->SetGuidance("Minimum value, expressed in unit");
  
  G4UIparameter* h1ValMax = new G4UIparameter("valMax", 'd', false);
  h1ValMax->SetGuidance("Maximum value, expressed in unit");
  
  G4UIparameter* h1ValUnit = new G4UIparameter("valUnit", 's', true);
  h1ValUnit->SetGuidance("The unit of valMin and valMax");
  h1ValUnit->SetDefaultValue("none");
 
  G4UIparameter* h1ValFcn = new G4UIparameter("valFcn", 's', true);
  h1ValFcn->SetParameterCandidates("log log10 exp none");
  G4String fcnGuidance = "The function applied to filled values (log, log10, exp, none).\n";
  fcnGuidance += "Note that the unit parameter cannot be omitted in this case,\n";
  fcnGuidance += "but none value should be used insted.";
  h1ValFcn->SetGuidance(fcnGuidance);
  h1ValFcn->SetDefaultValue("none");
 
  fSetH1Cmd = new G4UIcommand("/analysis/h1/set", this);
  fSetH1Cmd->SetGuidance("Set parameters for the 1D histogram of #Id :");
  fSetH1Cmd->SetGuidance("  nbins; valMin; valMax; unit (of vmin and vmax)");
  fSetH1Cmd->SetParameter(h1Id);
  fSetH1Cmd->SetParameter(h1Nbins);
  fSetH1Cmd->SetParameter(h1ValMin);
  fSetH1Cmd->SetParameter(h1ValMax);
  fSetH1Cmd->SetParameter(h1ValUnit);
  fSetH1Cmd->SetParameter(h1ValFcn);
  fSetH1Cmd->AvailableForStates(G4State_PreInit, G4State_Idle);
}  

//_____________________________________________________________________________
void G4AnalysisMessenger::SetH1TitleCmd()
{
  G4UIparameter* h1Id = new G4UIparameter("idTitle", 'i', false);
  h1Id->SetGuidance("Histogram id");
  h1Id->SetParameterRange("idTitle>=0");

  G4UIparameter* h1Title = new G4UIparameter("h1Title", 's', true);
  h1Title->SetGuidance("Histogram title");
  h1Title->SetDefaultValue("none");

  fSetH1TitleCmd = new G4UIcommand("/analysis/h1/setTitle", this);
  fSetH1TitleCmd->SetGuidance("Set title for the 1D histogram of #Id");
  fSetH1TitleCmd->SetParameter(h1Id);
  fSetH1TitleCmd->SetParameter(h1Title);
  fSetH1TitleCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
}  

//_____________________________________________________________________________
void G4AnalysisMessenger::SetH1XAxisCmd()
{
  G4UIparameter* h1Id = new G4UIparameter("idXaxis", 'i', false);
  h1Id->SetGuidance("Histogram id");
  h1Id->SetParameterRange("idXaxis>=0");

  G4UIparameter* h1XAxis = new G4UIparameter("h1Xaxis", 's', true);
  h1XAxis->SetGuidance("Histogram x-axis title");
  h1XAxis->SetDefaultValue("none");

  fSetH1XAxisCmd = new G4UIcommand("/analysis/h1/setXaxis", this);
  fSetH1XAxisCmd->SetGuidance("Set x-axis title for the 1D histogram of #Id");
  fSetH1XAxisCmd->SetParameter(h1Id);
  fSetH1XAxisCmd->SetParameter(h1XAxis);
  fSetH1XAxisCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
}  

//_____________________________________________________________________________
void G4AnalysisMessenger::SetH1YAxisCmd()
{
  G4UIparameter* h1Id = new G4UIparameter("idYaxis", 'i', false);
  h1Id->SetGuidance("Histogram id");
  h1Id->SetParameterRange("idYaxis>=0");

  G4UIparameter* h1YAxis = new G4UIparameter("h1Yaxis", 's', true);
  h1YAxis->SetGuidance("Histogram y-axis title");
  h1YAxis->SetDefaultValue("none");

  fSetH1YAxisCmd = new G4UIcommand("/analysis/h1/setYaxis", this);
  fSetH1YAxisCmd->SetGuidance("Set y-axis title for the 1D histogram of #Id");
  fSetH1YAxisCmd->SetParameter(h1Id);
  fSetH1YAxisCmd->SetParameter(h1YAxis);
  fSetH1YAxisCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
}  

//_____________________________________________________________________________
void G4AnalysisMessenger::SetH1ActivationCmd()
{
  G4UIparameter* h1Id = new G4UIparameter("idActivation", 'i', false);
  h1Id->SetGuidance("Histogram id");
  h1Id->SetParameterRange("idActivation>=0");

  G4UIparameter* h1Activation = new G4UIparameter("h1Activation", 's', true);
  h1Activation->SetGuidance("Histogram activation");
  h1Activation->SetDefaultValue("none");

  fSetH1ActivationCmd = new G4UIcommand("/analysis/h1/setActivation", this);
  fSetH1ActivationCmd->SetGuidance("Set activation for the 1D histogram of #Id");
  fSetH1ActivationCmd->SetParameter(h1Id);
  fSetH1ActivationCmd->SetParameter(h1Activation);
  fSetH1ActivationCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
}  

//_____________________________________________________________________________
void G4AnalysisMessenger::CreateH2Cmd()
{
  G4UIparameter* h2Name = new G4UIparameter("name", 's', false);
  h2Name->SetGuidance("Histogram name (label)");
  
  G4UIparameter* h2Title = new G4UIparameter("title", 's', false);
  h2Title->SetGuidance("Histogram title");

  G4UIparameter* h2xNbins0 = new G4UIparameter("xnbins0", 'i', true);
  h2xNbins0->SetGuidance("Number of x-bins (default = 100)");
  h2xNbins0->SetGuidance("Can be reset with /analysis/h2/set command");
  h2xNbins0->SetDefaultValue(100);
  
  G4UIparameter* h2xValMin0 = new G4UIparameter("xvalMin0", 'd', true);
  h2xValMin0->SetGuidance("Minimum x-value, expressed in unit (default = 0.)");
  h2xValMin0->SetGuidance("Can be reset with /analysis/h2/set command");
  h2xValMin0->SetDefaultValue(0.);
  
  G4UIparameter* h2xValMax0 = new G4UIparameter("xvalMax0", 'd', true);
  h2xValMax0->SetGuidance("Maximum x-value, expressed in unit (default = 1.)");
  h2xValMax0->SetGuidance("Can be reset with /analysis/h2/set command");
  h2xValMax0->SetDefaultValue(1.);

  G4UIparameter* h2xValUnit0 = new G4UIparameter("xvalUnit0", 's', true);
  h2xValUnit0->SetGuidance("The unit of xvalMin0 and xvalMax0");
  h2xValUnit0->SetDefaultValue("none");
  
  G4UIparameter* h2xValFcn0 = new G4UIparameter("xvalFcn0", 's', true);
  G4String fcnxGuidance = "The function applied to filled x-values (log, log10, exp, none).\n";
  fcnxGuidance += "Note that the unit parameter cannot be omitted in this case,\n";
  fcnxGuidance += "but none value should be used insted.";
  h2xValFcn0->SetGuidance(fcnxGuidance);
  h2xValFcn0->SetParameterCandidates("log log10 exp none");
  h2xValFcn0->SetDefaultValue("none");
  
  G4UIparameter* h2yNbins0 = new G4UIparameter("ynbins0", 'i', true);
  h2yNbins0->SetGuidance("Number of y-bins (default = 100)");
  h2yNbins0->SetGuidance("Can be reset with /analysis/h2/set command");
  h2yNbins0->SetDefaultValue(100);
  
  G4UIparameter* h2yValMin0 = new G4UIparameter("yvalMin0", 'd', true);
  h2yValMin0->SetGuidance("Minimum y-value, expressed in unit (default = 0.)");
  h2yValMin0->SetGuidance("Can be reset with /analysis/h2/set command");
  h2yValMin0->SetDefaultValue(0.);
  
  G4UIparameter* h2yValMax0 = new G4UIparameter("yvalMax0", 'd', true);
  h2yValMax0->SetGuidance("Maximum y-value, expressed in unit (default = 1.)");
  h2yValMax0->SetGuidance("Can be reset with /analysis/h2/set command");
  h2yValMax0->SetDefaultValue(1.);

  G4UIparameter* h2yValUnit0 = new G4UIparameter("yvalUnit0", 's', true);
  h2yValUnit0->SetGuidance("The unit of xvalMin0 and yvalMax0");
  h2yValUnit0->SetDefaultValue("none");
  
  G4UIparameter* h2yValFcn0 = new G4UIparameter("yvalFcn0", 's', true);
  G4String fcnyGuidance = "The function applied to filled x-values (log, log10, exp, none).\n";
  fcnyGuidance += "Note that the unit parameter cannot be omitted in this case,\n";
  fcnyGuidance += "but none value should be used insted.";
  h2yValFcn0->SetGuidance(fcnyGuidance);
  h2yValFcn0->SetParameterCandidates("log log10 exp none");
  h2yValFcn0->SetDefaultValue("none");
  
  fCreateH2Cmd = new G4UIcommand("/analysis/h2/create", this);
  fCreateH2Cmd->SetGuidance("Create 2D histogram");
  fCreateH2Cmd->SetParameter(h2Name);
  fCreateH2Cmd->SetParameter(h2Title);
  fCreateH2Cmd->SetParameter(h2xNbins0);
  fCreateH2Cmd->SetParameter(h2xValMin0);
  fCreateH2Cmd->SetParameter(h2xValMax0);
  fCreateH2Cmd->SetParameter(h2xValUnit0);
  fCreateH2Cmd->SetParameter(h2xValFcn0);
  fCreateH2Cmd->SetParameter(h2yNbins0);
  fCreateH2Cmd->SetParameter(h2yValMin0);
  fCreateH2Cmd->SetParameter(h2yValMax0);
  fCreateH2Cmd->SetParameter(h2yValUnit0);
  fCreateH2Cmd->SetParameter(h2yValFcn0);
  fCreateH2Cmd->AvailableForStates(G4State_PreInit, G4State_Idle);
}  


//_____________________________________________________________________________
void G4AnalysisMessenger::SetH2Cmd()
{
  G4UIparameter* h2Id = new G4UIparameter("id", 'i', false);
  h2Id->SetGuidance("Histogram id");
  h2Id->SetParameterRange("id>=0");
  
  G4UIparameter* h2xNbins = new G4UIparameter("xnbins", 'i', false);
  h2xNbins->SetGuidance("Number of x-bins");
  
  G4UIparameter* h2xValMin = new G4UIparameter("xvalMin", 'd', false);
  h2xValMin->SetGuidance("Minimum x-value, expressed in unit");
  
  G4UIparameter* h2xValMax = new G4UIparameter("xvalMax", 'd', false);
  h2xValMax->SetGuidance("Maximum x-value, expressed in unit");
  
  G4UIparameter* h2xValFcn = new G4UIparameter("xvalFcn", 's', false);
  h2xValFcn->SetParameterCandidates("log log10 exp none");
  G4String fcnxGuidance = "The function applied to filled x-values (log, log10, exp, none).\n";
  fcnxGuidance += "Note that the unit parameter cannot be omitted in this case,\n";
  fcnxGuidance += "but none value should be used insted.";
  h2xValFcn->SetGuidance(fcnxGuidance);
  h2xValFcn->SetDefaultValue("none");
 
  G4UIparameter* h2yValUnit = new G4UIparameter("yvalUnit", 's', false);
  h2yValUnit->SetGuidance("The unit of yvalMin and yvalMax");
  h2yValUnit->SetDefaultValue("none");
 
  G4UIparameter* h2yNbins = new G4UIparameter("nybins", 'i', false);
  h2yNbins->SetGuidance("Number of y-bins");
  
  G4UIparameter* h2yValMin = new G4UIparameter("yvalMin", 'd', false);
  h2yValMin->SetGuidance("Minimum y-value, expressed in unit");
  
  G4UIparameter* h2yValMax = new G4UIparameter("yvalMax", 'd', false);
  h2yValMax->SetGuidance("Maximum y-value, expressed in unit");
  
  G4UIparameter* h2xValUnit = new G4UIparameter("xvalUnit", 's', true);
  h2xValUnit->SetGuidance("The unit of xvalMin and xvalMax");
  h2xValUnit->SetDefaultValue("none");
 
  G4UIparameter* h2yValFcn = new G4UIparameter("yvalFcn", 's', false);
  h2yValFcn->SetParameterCandidates("log log10 exp none");
  G4String fcnyGuidance = "The function applied to filled y-values (log, log10, exp, none).\n";
  fcnyGuidance += "Note that the unit parameter cannot be omitted in this case,\n";
  fcnyGuidance += "but none value should be used insted.";
  h2yValFcn->SetGuidance(fcnyGuidance);
  h2yValFcn->SetDefaultValue("none");
 
  fSetH2Cmd = new G4UIcommand("/analysis/h2/set", this);
  fSetH2Cmd->SetGuidance("Set parameters for the 2D histogram of #Id :");
  fSetH2Cmd->SetGuidance("  nbins; valMin; valMax; unit (of vmin and vmax)");
  fSetH2Cmd->SetParameter(h2Id);
  fSetH2Cmd->SetParameter(h2xNbins);
  fSetH2Cmd->SetParameter(h2xValMin);
  fSetH2Cmd->SetParameter(h2xValMax);
  fSetH2Cmd->SetParameter(h2xValUnit);
  fSetH2Cmd->SetParameter(h2xValFcn);
  fSetH2Cmd->SetParameter(h2yNbins);
  fSetH2Cmd->SetParameter(h2yValMin);
  fSetH2Cmd->SetParameter(h2yValMax);
  fSetH2Cmd->SetParameter(h2yValUnit);
  fSetH2Cmd->SetParameter(h2yValFcn);
  fSetH2Cmd->AvailableForStates(G4State_PreInit, G4State_Idle);
}  

//_____________________________________________________________________________
void G4AnalysisMessenger::SetH2TitleCmd()
{
  G4UIparameter* h2Id = new G4UIparameter("idTitle", 'i', false);
  h2Id->SetGuidance("Histogram id");
  h2Id->SetParameterRange("idTitle>=0");

  G4UIparameter* h2Title = new G4UIparameter("h2Title", 's', true);
  h2Title->SetGuidance("Histogram title");
  h2Title->SetDefaultValue("none");

  fSetH2TitleCmd = new G4UIcommand("/analysis/h2/setTitle", this);
  fSetH2TitleCmd->SetGuidance("Set title for the 2D histogram of #Id");
  fSetH2TitleCmd->SetParameter(h2Id);
  fSetH2TitleCmd->SetParameter(h2Title);
  fSetH2TitleCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
}  

//_____________________________________________________________________________
void G4AnalysisMessenger::SetH2XAxisCmd()
{
  G4UIparameter* h2Id = new G4UIparameter("idXaxis", 'i', false);
  h2Id->SetGuidance("Histogram id");
  h2Id->SetParameterRange("idXaxis>=0");

  G4UIparameter* h2XAxis = new G4UIparameter("h2Xaxis", 's', true);
  h2XAxis->SetGuidance("Histogram x-axis title");
  h2XAxis->SetDefaultValue("none");

  fSetH2XAxisCmd = new G4UIcommand("/analysis/h2/setXaxis", this);
  fSetH2XAxisCmd->SetGuidance("Set x-axis title for the 2D histogram of #Id");
  fSetH2XAxisCmd->SetParameter(h2Id);
  fSetH2XAxisCmd->SetParameter(h2XAxis);
  fSetH2XAxisCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
}  

//_____________________________________________________________________________
void G4AnalysisMessenger::SetH2YAxisCmd()
{
  G4UIparameter* h2Id = new G4UIparameter("idYaxis", 'i', false);
  h2Id->SetGuidance("Histogram id");
  h2Id->SetParameterRange("idYaxis>=0");

  G4UIparameter* h2YAxis = new G4UIparameter("h2Yaxis", 's', true);
  h2YAxis->SetGuidance("Histogram y-axis title");
  h2YAxis->SetDefaultValue("none");

  fSetH2YAxisCmd = new G4UIcommand("/analysis/h2/setYaxis", this);
  fSetH2YAxisCmd->SetGuidance("Set y-axis title for the 2D histogram of #Id");
  fSetH2YAxisCmd->SetParameter(h2Id);
  fSetH2YAxisCmd->SetParameter(h2YAxis);
  fSetH2YAxisCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
}  

//_____________________________________________________________________________
void G4AnalysisMessenger::SetH2ZAxisCmd()
{
  G4UIparameter* h2Id = new G4UIparameter("idYaxis", 'i', false);
  h2Id->SetGuidance("Histogram id");
  h2Id->SetParameterRange("idYaxis>=0");

  G4UIparameter* h2YAxis = new G4UIparameter("h2Yaxis", 's', true);
  h2YAxis->SetGuidance("Histogram y-axis title");
  h2YAxis->SetDefaultValue("none");

  fSetH2ZAxisCmd = new G4UIcommand("/analysis/h2/setYaxis", this);
  fSetH2ZAxisCmd->SetGuidance("Set y-axis title for the 2D histogram of #Id");
  fSetH2ZAxisCmd->SetParameter(h2Id);
  fSetH2ZAxisCmd->SetParameter(h2YAxis);
  fSetH2ZAxisCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
}  

//_____________________________________________________________________________
void G4AnalysisMessenger::SetH2ActivationCmd()
{
  G4UIparameter* h2Id = new G4UIparameter("idActivation", 'i', false);
  h2Id->SetGuidance("Histogram id");
  h2Id->SetParameterRange("idActivation>=0");

  G4UIparameter* h2Activation = new G4UIparameter("h2Activation", 's', true);
  h2Activation->SetGuidance("Histogram activation");
  h2Activation->SetDefaultValue("none");

  fSetH2ActivationCmd = new G4UIcommand("/analysis/h2/setActivation", this);
  fSetH2ActivationCmd->SetGuidance("Set activation for the 2D histogram of #Id");
  fSetH2ActivationCmd->SetParameter(h2Id);
  fSetH2ActivationCmd->SetParameter(h2Activation);
  fSetH2ActivationCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
}  

//_____________________________________________________________________________
void G4AnalysisMessenger::SetNewValue(G4UIcommand* command, G4String newValues)
{
  if ( command == fSetFileNameCmd ) {
    G4cout << "Set file name: " << newValues << G4endl;
    fManager->SetFileName(newValues);
  }  
  else if ( command == fSetHistoDirNameCmd ) {
    fManager->SetHistoDirectoryName(newValues);
  }  
  else if ( command == fSetNtupleDirNameCmd ) {
    fManager->SetNtupleDirectoryName(newValues);
  }  
  else if ( command == fSetActivationCmd ) {
    fManager->SetActivation(fSetActivationCmd->GetNewBoolValue(newValues));
  }  
  else if ( command == fVerboseCmd ) {
    fManager->SetVerboseLevel(fVerboseCmd->GetNewIntValue(newValues));
  }  
  else if ( command == fCreateH1Cmd ) { 
    G4String name, title;
    G4int nbins; 
    G4double vmin,vmax; 
    G4String sunit;
    G4String sfcn;
    std::istringstream is(newValues.data());
    is >> name >> title >> nbins >> vmin >> vmax >> sunit >> sfcn;
    fManager->CreateH1(name, title, nbins, vmin, vmax, sunit, sfcn);     
  }
  else if ( command == fSetH1Cmd ) {
    G4int id; 
    G4int nbins; 
    G4double vmin, vmax; 
    G4String sunit;
    G4String sfcn;
    std::istringstream is(newValues.data());
    is >> id >> nbins >> vmin >> vmax >> sunit >> sfcn;
    fManager->SetH1(id, nbins, vmin, vmax, sunit, sfcn);     
  }
  else if ( command == fSetH1AsciiCmd ) {
    G4int id = fSetH1AsciiCmd->GetNewIntValue(newValues);
    fManager->SetAscii(G4VAnalysisManager::kH1, id, true); 
  }      
  else if ( command == fSetH1TitleCmd ) {
    G4int id; 
    G4String title;
    std::istringstream is(newValues.data());
    is >> id;
    getline(is, title);
    fManager->SetH1Title(id, title);     
  }
  else if ( command == fSetH1XAxisCmd ) {
    G4int id; 
    G4String xaxis;
    std::istringstream is(newValues.data());
    is >> id;
    getline(is, xaxis);
    fManager->SetH1XAxisTitle(id, xaxis);     
  }
  else if ( command == fSetH1YAxisCmd ) {
    G4int id; 
    G4String yaxis;
    std::istringstream is(newValues.data());
    is >> id;
    getline(is, yaxis);
    fManager->SetH1YAxisTitle(id, yaxis);     
  }
  else if ( command == fSetH1ActivationCmd ) {
    G4int id; 
    G4String sactivation;
    std::istringstream is(newValues.data());
    is >> id >> sactivation;
    G4bool activation = G4UIcommand::ConvertToBool(sactivation);
    fManager->SetActivation(G4VAnalysisManager::kH1, id, activation);     
  }
  else if ( command == fSetH1ActivationAllCmd ) {
    G4bool activation = fSetH1ActivationAllCmd->GetNewBoolValue(newValues);
    fManager->SetActivation(G4VAnalysisManager::kH1, activation);
  }  
  else if ( command == fCreateH2Cmd ) { 
    G4String name, title;
    G4int xnbins, ynbins; 
    G4double xvmin, xvmax, yvmin, yvmax; 
    G4String xsunit,xsfcn, ysunit, ysfcn;
    std::istringstream is(newValues.data());
    is >> name >> title 
       >> xnbins >> xvmin >> xvmax >> xsunit >> xsfcn
       >> ynbins >> yvmin >> yvmax >> ysunit >> ysfcn;
    fManager->CreateH2(name, title, 
                       xnbins, xvmin, xvmax, ynbins, yvmin, yvmax, 
                       ysunit, ysfcn, ysunit, ysfcn);     
  }
  else if ( command == fSetH2Cmd ) {
    G4int id; 
    G4int xnbins, ynbins; 
    G4double xvmin, xvmax, yvmin, yvmax; 
    G4String xsunit,xsfcn, ysunit, ysfcn;
    std::istringstream is(newValues.data());
    is >> id 
       >> xnbins >> xvmin >> xvmax >> xsunit >> xsfcn  
       >> ynbins >> yvmin >> yvmax >> ysunit >> ysfcn;
    fManager->SetH2(id, 
                    xnbins, xvmin, xvmax, ynbins, yvmin, yvmax, 
                    ysunit, ysfcn, ysunit, ysfcn);     
  }
  else if ( command == fSetH2AsciiCmd ) {
    G4int id = fSetH2AsciiCmd->GetNewIntValue(newValues);
    fManager->SetAscii(G4VAnalysisManager::kH2, id, true); 
  }      
  else if ( command == fSetH2TitleCmd ) {
    G4int id; 
    G4String title;
    std::istringstream is(newValues.data());
    is >> id;
    getline(is, title);
    fManager->SetH2Title(id, title);     
  }
  else if ( command == fSetH2XAxisCmd ) {
    G4int id; 
    G4String xaxis;
    std::istringstream is(newValues.data());
    is >> id;
    getline(is, xaxis);
    fManager->SetH2XAxisTitle(id, xaxis);     
  }
  else if ( command == fSetH2YAxisCmd ) {
    G4int id; 
    G4String yaxis;
    std::istringstream is(newValues.data());
    is >> id;
    getline(is, yaxis);
    fManager->SetH2YAxisTitle(id, yaxis);     
  }
  else if ( command == fSetH2ZAxisCmd ) {
    G4int id; 
    G4String zaxis;
    std::istringstream is(newValues.data());
    is >> id;
    getline(is, zaxis);
    fManager->SetH2ZAxisTitle(id, zaxis);     
  }
  else if ( command == fSetH2ActivationCmd ) {
    G4int id; 
    G4String sactivation;
    std::istringstream is(newValues.data());
    is >> id >> sactivation;
    G4bool activation = G4UIcommand::ConvertToBool(sactivation);
    fManager->SetActivation(G4VAnalysisManager::kH2, id, activation);     
  }
  else if ( command == fSetH2ActivationAllCmd ) {
    G4bool activation = fSetH2ActivationAllCmd->GetNewBoolValue(newValues);
    fManager->SetActivation(G4VAnalysisManager::kH2, activation);
  }  
}  
