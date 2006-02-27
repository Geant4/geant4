// $Id: NISTmaterials.cc,v 1.1 2006-02-27 09:48:56 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   NISTmaterials.cc
//
//                                         2005 Q
// ====================================================================
#include "G4NistManager.hh"

////////////////
void Construct()
////////////////
{
  G4NistManager* NistMgr = G4NistManager::Instance();

  // define typically used materials
  NistMgr-> FindOrBuildMaterial("G4_Al");
  NistMgr-> FindOrBuildMaterial("G4_Si");
  NistMgr-> FindOrBuildMaterial("G4_Ar");
  NistMgr-> FindOrBuildMaterial("G4_Cu");
  NistMgr-> FindOrBuildMaterial("G4_Fe");
  NistMgr-> FindOrBuildMaterial("G4_Ge");
  NistMgr-> FindOrBuildMaterial("G4_Ag");
  NistMgr-> FindOrBuildMaterial("G4_W");
  NistMgr-> FindOrBuildMaterial("G4_Au");
  NistMgr-> FindOrBuildMaterial("G4_Pb");

  NistMgr-> FindOrBuildMaterial("G4_AIR");
  NistMgr-> FindOrBuildMaterial("G4_Galactic");
  NistMgr-> FindOrBuildMaterial("G4_WATER");
  NistMgr-> FindOrBuildMaterial("G4_CESIUM_IODIDE");
  NistMgr-> FindOrBuildMaterial("G4_SODIUM_IODIDE");
  NistMgr-> FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");
  NistMgr-> FindOrBuildMaterial("G4_MYLAR");
  //NistMgr-> FindOrBuildMaterial("G4_GRAPHITE"); 
  // bug in G4 (element isnot defined. should be fixed.)

}
