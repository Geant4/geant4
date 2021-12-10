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

#ifndef G4PLOTTERMANAGER_HH
#define G4PLOTTERMANAGER_HH

#include "G4Plotter.hh"
#include "G4UImessenger.hh"

#include <vector>
#include <utility>

class G4PlotterManager {
public:
  static G4PlotterManager& GetInstance ();
  G4Plotter& GetPlotter(const G4String& a_name);
  void List() const;

  using StyleItem = std::pair<G4String,G4String>;
  using Style = std::vector<StyleItem>;
  using NamedStyle = std::pair<G4String,Style>;
  using Styles = std::vector<NamedStyle>;
  const Styles& GetStyles() const {return fStyles;}
  Styles& GetStyles() {return fStyles;}

private:
  G4PlotterManager();
  virtual ~G4PlotterManager();
  G4PlotterManager (const G4PlotterManager&);
  G4PlotterManager& operator = (const G4PlotterManager&);

  void ListStyles() const;
  Style* FindStyle(const G4String& name);
  void SelectStyle(const G4String& style);
  void RemoveStyle(const G4String& name);
  void PrintStyle(const G4String&) const;
  void AddStyleParameter(const G4String& param,const G4String& value);

  typedef std::pair<G4String,G4Plotter> NamedPlotter;
  std::vector<NamedPlotter> fPlotters;

  G4String fCurrentStyle;
  Styles fStyles;

  class Messenger: public G4UImessenger {
  public:  
    Messenger(G4PlotterManager& aPlotterManager):fPlotterManager(aPlotterManager) {
      G4UIparameter* parameter;
      //////////////////////////////////////////////////////////
      //////////////////////////////////////////////////////////
      //////////////////////////////////////////////////////////
      remove_style = new G4UIcommand("/vis/plotter/style/remove",this);
      remove_style->SetGuidance("Remove a named style.");

      parameter = new G4UIparameter("name",'s',false);
      remove_style->SetParameter(parameter);

      //////////////////////////////////////////////////////////
      //////////////////////////////////////////////////////////
      select_style = new G4UIcommand("/vis/plotter/style/select",this);
      select_style->SetGuidance("Select a named style for further style/add commands.");
      select_style->SetGuidance("If not existing, the named style is created.");

      parameter = new G4UIparameter("name",'s',false);
      select_style->SetParameter(parameter);

      //////////////////////////////////////////////////////////
      //////////////////////////////////////////////////////////
      add_style_parameter = new G4UIcommand("/vis/plotter/style/add",this);
      add_style_parameter->SetGuidance("Add a (parameter,value) to the current named style.");

      parameter = new G4UIparameter("parameter",'s',false);
      add_style_parameter->SetParameter (parameter);
  
      parameter = new G4UIparameter("value",'s',false);
      add_style_parameter->SetParameter (parameter);

      //////////////////////////////////////////////////////////
      //////////////////////////////////////////////////////////
      list_styles = new G4UIcommand("/vis/plotter/style/list", this);
      list_styles->SetGuidance("List known not embedded styles.");
      
      //////////////////////////////////////////////////////////
      //////////////////////////////////////////////////////////
      print_style = new G4UIcommand("/vis/plotter/style/print", this);
      print_style->SetGuidance("Print a style.");

      parameter = new G4UIparameter("style",'s',false);
      print_style->SetParameter (parameter);
    }
    virtual ~Messenger() {
      delete remove_style;
      delete select_style;
      delete add_style_parameter;
      delete list_styles;
      delete print_style;
    }
    virtual void SetNewValue(G4UIcommand*,G4String);
  private:  
    G4PlotterManager& fPlotterManager;
    G4UIcommand* remove_style; 
    G4UIcommand* select_style;
    G4UIcommand* add_style_parameter;
    G4UIcommand* list_styles;
    G4UIcommand* print_style;
  };
  
  Messenger* fMessenger;
};

#endif
