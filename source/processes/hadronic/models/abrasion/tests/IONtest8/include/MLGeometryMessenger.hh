#ifndef MLGeometryMessenger_h
#define MLGeometryMessenger_h 1
////////////////////////////////////////////////////////////////////////////////
//
#include "globals.hh"
#include "G4UImessenger.hh"

#include "G4UIdirectory.hh"
#include "G4UIcommand.hh"
#include "G4UIparameter.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithAString.hh"

class MLGeometryConstruction;
////////////////////////////////////////////////////////////////////////////////
//
class MLGeometryMessenger: public G4UImessenger
{
public:
  MLGeometryMessenger (MLGeometryConstruction*);
  ~MLGeometryMessenger ();

  void SetNewValue (G4UIcommand*, G4String);
    
private:
  MLGeometryConstruction    *MLGeometry;

  G4UIdirectory             *MLdetDir;

  G4UIdirectory             *LayerDir;
  G4UIcmdWithAString        *ShapeCmd;
  G4UIcmdWithAnInteger      *DeleteCmd;
  G4UIcmdWithAnInteger      *ListCmd;
  G4UIcommand               *AddCmd;

  G4UIcmdWithoutParameter   *UseDefaultCmd;
  G4UIcmdWithoutParameter   *UpdateCmd;

  G4UIcmdWithAnInteger      *AddELayerCmd;
  G4UIcmdWithAnInteger      *DeleteELayerCmd;
  G4UIcmdWithoutParameter   *ListELayerCmd;

  G4UIcmdWithAnInteger      *AddFLayerCmd;
  G4UIcmdWithAnInteger      *DeleteFLayerCmd;
  G4UIcmdWithoutParameter   *ListFLayerCmd;
};
////////////////////////////////////////////////////////////////////////////////
#endif

