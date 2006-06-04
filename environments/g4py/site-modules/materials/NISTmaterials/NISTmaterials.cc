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
// $Id: NISTmaterials.cc,v 1.3 2006-06-04 21:36:35 kmura Exp $
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
