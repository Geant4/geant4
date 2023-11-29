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
// Guy Barrand 12 October 2021
//

#include "G4PlotterManager.hh"
#include "G4ios.hh"

#include <tools/forit>
#include <tools/tokenize>

G4PlotterManager& G4PlotterManager::GetInstance () {
  static G4PlotterManager s_instance;
  return s_instance;
}

G4PlotterManager::G4PlotterManager():fMessenger(0) {
  fMessenger = new Messenger(*this);
}

G4PlotterManager::~G4PlotterManager() {
  delete fMessenger;
}

G4Plotter& G4PlotterManager::GetPlotter(const G4String& a_name) {
  tools_vforit(NamedPlotter,fPlotters,it) {
    if((*it).first==a_name) {
      return (*it).second;
    }
  }
  fPlotters.push_back(NamedPlotter(a_name,G4Plotter()));
  return fPlotters.back().second;
}

void G4PlotterManager::List() const {
  tools_vforcit(NamedPlotter,fPlotters,it) {
    G4cout << (*it).first << G4endl;
  }
}

//////////////////////////////////////////////////////////////////
/// styles: //////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////

void G4PlotterManager::ListStyles() const {
  tools_vforcit(NamedStyle,fStyles,it) {
    G4cout << (*it).first << G4endl;
  }
}  

G4PlotterManager::Style* G4PlotterManager::FindStyle(const G4String& a_name) {
  tools_vforit(NamedStyle,fStyles,it){
    if((*it).first==a_name) return &((*it).second);
  }
  return 0;
}

void G4PlotterManager::SelectStyle(const G4String& a_name) {
  if(!FindStyle(a_name)) {
    fStyles.push_back(NamedStyle(a_name,Style()));
  }
  fCurrentStyle = a_name;
}  

void G4PlotterManager::RemoveStyle(const G4String& a_name) {
  tools_vforit(NamedStyle,fStyles,it) {
    if((*it).first==a_name) {
      fStyles.erase(it);
      if(fCurrentStyle==a_name) fCurrentStyle.clear();
      return;
    }
  }
}  

void G4PlotterManager::PrintStyle(const G4String& a_name) const {
  tools_vforcit(NamedStyle,fStyles,it) {
    if((*it).first==a_name) {
      G4cout << (*it).first << ":" << G4endl;
      tools_vforcit(StyleItem,(*it).second,its) {
        G4cout << " " << (*its).first << " " << (*its).second << G4endl;
      }
    }
  }
}  

void G4PlotterManager::AddStyleParameter(const G4String& a_parameter,const G4String& a_value) {
  Style* _style = FindStyle(fCurrentStyle);
  if(!_style) {
    G4cout << "G4PlotterManager::AddStyleParameter: style " << fCurrentStyle << " not found." << G4endl;
    return;
  }
  tools_vforit(StyleItem,(*_style),it) {
    if((*it).first==a_parameter) {
      (*it).second = a_value;
      return;
    }
  }
  _style->push_back(StyleItem(a_parameter,a_value));
}  

void G4PlotterManager::Messenger::SetNewValue(G4UIcommand* a_cmd,G4String a_value) {
  std::vector<std::string> args;
  tools::double_quotes_tokenize(a_value,args);
  if(args.size()!=a_cmd->GetParameterEntries()) return;
  if(a_cmd==select_style) {
    fPlotterManager.SelectStyle(args[0]);
  } else if(a_cmd==add_style_parameter) {
    fPlotterManager.AddStyleParameter(args[0],args[1]);
  } else if(a_cmd==remove_style) {
    fPlotterManager.RemoveStyle(args[0]);
  } else if(a_cmd==list_styles) {
    G4cout << "default (embedded)." << G4endl;
    G4cout << "ROOT_default (embedded)." << G4endl;
    G4cout << "hippodraw (embedded)." << G4endl;
    fPlotterManager.ListStyles();
  } else if(a_cmd==print_style) {
    fPlotterManager.PrintStyle(args[0]);
  }
}
