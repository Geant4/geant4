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
// John Allison  6th October 2020

#ifndef G4TOOLSSGSCENEHANDLER_HH
#define G4TOOLSSGSCENEHANDLER_HH

#include "G4VSceneHandler.hh"

#include "G4VVisCommand.hh"

#include <vector>

#include <tools/sg/separator>

namespace tools {namespace sg {class base_freetype;}}
namespace tools {namespace sg {class plots;}}

class G4ToolsSGNode;

class G4ToolsSGSceneHandler: public G4VSceneHandler {
  typedef G4VSceneHandler parent;
public:  
  virtual void AddPrimitive(const G4Polyline&);
  virtual void AddPrimitive(const G4Text&);
  virtual void AddPrimitive(const G4Circle&);
  virtual void AddPrimitive(const G4Square&);
  virtual void AddPrimitive(const G4Polymarker&);
  virtual void AddPrimitive(const G4Polyhedron&);
  virtual void AddPrimitive(const G4Plotter&);

  using G4VSceneHandler::AddCompound;
  virtual void AddCompound(const G4Mesh&);

  virtual void ClearStore ();
  virtual void ClearTransientStore ();

  G4ToolsSGSceneHandler(G4VGraphicsSystem& system,const G4String& name);
  virtual ~G4ToolsSGSceneHandler();

  tools::sg::separator& GetTransient2DObjects() {return fpTransient2DObjects;}
  tools::sg::separator& GetPersistent2DObjects() {return fpPersistent2DObjects;}
  tools::sg::separator& GetTransient3DObjects() {return fpTransient3DObjects;}
  tools::sg::separator& GetPersistent3DObjects() {return fpPersistent3DObjects;}
  
  void TouchPlotters(tools::sg::node&);

protected:  
  G4ToolsSGSceneHandler(const G4ToolsSGSceneHandler&);
  G4ToolsSGSceneHandler& operator=(const G4ToolsSGSceneHandler&){return *this;}

protected:
  void CreateSG();
  void EstablishBaseNodes();
  tools::sg::separator* GetOrCreateNode();  // For next solid or primitive

  void SetPlotterHistograms(tools::sg::plots&);

  static G4int fSceneIdCount;

  tools::sg::separator fpTransient2DObjects;
  tools::sg::separator fpPersistent2DObjects;
  
  tools::sg::separator fpTransient3DObjects;
  tools::sg::separator fpPersistent3DObjects;

  std::vector<G4ToolsSGNode*> fpPhysicalVolumeObjects;  // Multiple worlds

  tools::sg::base_freetype* fFreetypeNode;
  
  using Region_h1 = std::pair<unsigned int,int>;
  using Region_h2 = std::pair<unsigned int,int>;
  std::vector<Region_h1> fRegionH1s;
  std::vector<Region_h2> fRegionH2s;

  class Messenger: public G4VVisCommand {
  public:  
    static void Create() {static Messenger s_messenger;}
    virtual void SetNewValue(G4UIcommand*,G4String);
  private:  
    Messenger() {
      //////////////////////////////////////////////////////////
      //////////////////////////////////////////////////////////
      print_plotter_params = new G4UIcommand("/vis/tsg/plotter/printParameters", this);
      print_plotter_params->SetGuidance("Print available tools::sg::plotter parameters.");
    }
    virtual ~Messenger() {
      delete print_plotter_params;
    }
    G4UIcommand* print_plotter_params;
  };

};

#endif
