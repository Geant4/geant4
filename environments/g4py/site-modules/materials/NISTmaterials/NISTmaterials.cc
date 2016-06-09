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
// $Id: NISTmaterials.cc,v 1.4 2006-06-29 15:29:32 gunter Exp $
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
