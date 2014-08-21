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

// Author: Ivana Hrivnacova, 24/07/2014  (ivana@ipno.in2p3.fr)

#include "G4H3DummyManager.hh"

//_____________________________________________________________________________
G4H3DummyManager::G4H3DummyManager(const G4AnalysisManagerState& state)
 : G4VH3Manager(state)
{
}

//_____________________________________________________________________________
G4H3DummyManager::~G4H3DummyManager()
{  
}

// 
// protected methods
//

//_____________________________________________________________________________
G4int G4H3DummyManager::CreateH3(const G4String& /*name*/, 
                               const G4String& /*title*/, 
                               G4int /*nxbins*/, 
                               G4double /*xmin*/, G4double /*xmax*/,
                               G4int /*nybins*/, 
                               G4double /*ymin*/, G4double /*ymax*/,
                               G4int /*nzbins*/, 
                               G4double /*zmin*/, G4double /*zmax*/,
                               const G4String& /*xunitName*/, 
                               const G4String& /*yunitName*/, 
                               const G4String& /*zunitName*/, 
                               const G4String& /*xfcnName*/,
                               const G4String& /*yfcnName*/,
                               const G4String& /*zfcnName*/,
                               const G4String& /*xbinScheme*/,
                               const G4String& /*ybinScheme*/,
                               const G4String& /*zbinScheme*/)
{
  ExceptionForHistograms("CreateH3");
  return 0;
}                                         

//_____________________________________________________________________________
G4int G4H3DummyManager::CreateH3(const G4String& /*name*/, 
                               const G4String& /*title*/,
                               const std::vector<G4double>& /*xedges*/,
                               const std::vector<G4double>& /*yedges*/,
                               const std::vector<G4double>& /*zedges*/,
                               const G4String& /*xunitName*/, 
                               const G4String& /*yunitName*/, 
                               const G4String& /*zunitName*/, 
                               const G4String& /*xfcnName*/,
                               const G4String& /*yfcnName*/,
                               const G4String& /*zfcnName*/)
                               
{
  ExceptionForHistograms("CreateH3");
  return 0;
}                                         

//_____________________________________________________________________________
G4bool G4H3DummyManager::SetH3(G4int /*id*/,
                               G4int /*nxbins*/, 
                               G4double /*xmin*/, G4double /*xmax*/,
                               G4int /*nybins*/, 
                               G4double /*ymin*/, G4double /*ymax*/,
                               G4int /*nzbins*/, 
                               G4double /*zmin*/, G4double /*zmax*/,
                               const G4String& /*xunitName*/, 
                               const G4String& /*yunitName*/, 
                               const G4String& /*zunitName*/, 
                               const G4String& /*xfcnName*/,
                               const G4String& /*yfcnName*/,
                               const G4String& /*zfcnName*/,
                               const G4String& /*xbinScheme*/,
                               const G4String& /*ybinScheme*/,
                               const G4String& /*zbinScheme*/)
{                                
  ExceptionForHistograms("SetH3");
  return false;
}
                                  
//_____________________________________________________________________________
G4bool G4H3DummyManager::SetH3(G4int /*id*/,
                               const std::vector<G4double>& /*xedges*/,
                               const std::vector<G4double>& /*yedges*/,
                               const std::vector<G4double>& /*zedges*/,
                               const G4String& /*xunitName*/, 
                               const G4String& /*yunitName*/, 
                               const G4String& /*zunitName*/, 
                               const G4String& /*xfcnName*/,
                               const G4String& /*yfcnName*/,
                               const G4String& /*zfcnName*/)
{                                
  ExceptionForHistograms("SetH3");
  return false;
}
                                  
//_____________________________________________________________________________
G4bool G4H3DummyManager::ScaleH3(G4int /*id*/, G4double /*factor*/)
{
  ExceptionForHistograms("ScaleH3");
  return false;
}  
                           
//_____________________________________________________________________________
G4bool G4H3DummyManager::FillH3(G4int /*id*/, 
                                G4double /*xvalue*/, G4double /*yvalue*/,
                                G4double /*zvalue*/,
                                G4double /*weight*/)
{
  ExceptionForHistograms("FillH3");
  return false;
}

//_____________________________________________________________________________
G4int  G4H3DummyManager::GetH3Id(const G4String& /*name*/, G4bool /*warn*/) const
{
  ExceptionForHistograms("GetH3Id");
  return 0;
}  
                                      
//_____________________________________________________________________________
G4int G4H3DummyManager::GetH3Nxbins(G4int /*id*/) const
{
  ExceptionForHistograms("GetH3Nxbins");
  return 0;
}  

//_____________________________________________________________________________
G4double G4H3DummyManager::GetH3Xmin(G4int /*id*/) const
{
  ExceptionForHistograms("GetH3Xmin");
  return 0;
}  

//_____________________________________________________________________________
G4double G4H3DummyManager::GetH3Xmax(G4int /*id*/) const
{
  ExceptionForHistograms("GetH3Xmax");
  return 0;
}  

//_____________________________________________________________________________
G4double G4H3DummyManager::GetH3XWidth(G4int /*id*/) const
{
  ExceptionForHistograms("GetH3XWidth");
  return 0;
}  

//_____________________________________________________________________________
G4int G4H3DummyManager::GetH3Nybins(G4int /*id*/) const
{
  ExceptionForHistograms("GetH3Nybins");
  return 0;
}  

//_____________________________________________________________________________
G4double G4H3DummyManager::GetH3Ymin(G4int /*id*/) const
{
  ExceptionForHistograms("");
  return 0;
}  

//_____________________________________________________________________________
G4double G4H3DummyManager::GetH3Ymax(G4int /*id*/) const
{
  ExceptionForHistograms("GetH3Ymax");
  return 0;
}  

//_____________________________________________________________________________
G4double G4H3DummyManager::GetH3YWidth(G4int /*id*/) const
{
  ExceptionForHistograms("GetH3YWidth");
  return 0;
}  

//_____________________________________________________________________________
G4int G4H3DummyManager::GetH3Nzbins(G4int /*id*/) const
{
  ExceptionForHistograms("GetH3Nzbins");
  return 0;
}  

//_____________________________________________________________________________
G4double G4H3DummyManager::GetH3Zmin(G4int /*id*/) const
{
  ExceptionForHistograms("GetH3Zmin");
  return 0;
}  

//_____________________________________________________________________________
G4double G4H3DummyManager::GetH3Zmax(G4int /*id*/) const
{
  ExceptionForHistograms("GetH3Zmax");
  return 0;
}  

//_____________________________________________________________________________
G4double G4H3DummyManager::GetH3ZWidth(G4int /*id*/) const
{
  ExceptionForHistograms("GetH3ZWidth");
  return 0;
}  

//_____________________________________________________________________________
G4bool G4H3DummyManager::SetH3Title(G4int /*id*/, const G4String& /*title*/)
{
  ExceptionForHistograms("SetH3Title");
  return false;
}  

//_____________________________________________________________________________
G4bool G4H3DummyManager::SetH3XAxisTitle(G4int /*id*/, const G4String& /*title*/)
{
  ExceptionForHistograms("SetH3XAxisTitle");
  return false;
}  

//_____________________________________________________________________________
G4bool G4H3DummyManager::SetH3YAxisTitle(G4int /*id*/, const G4String& /*title*/)
{
  ExceptionForHistograms("SetH3YAxisTitle");
  return false;
}  

//_____________________________________________________________________________
G4bool G4H3DummyManager::SetH3ZAxisTitle(G4int /*id*/, const G4String& /*title*/)
{
  ExceptionForHistograms("SetH3ZAxisTitle");
  return false;
}  

//_____________________________________________________________________________
G4String G4H3DummyManager::GetH3Title(G4int /*id*/) const
{
  ExceptionForHistograms("GetH3Title");
  return "";
}  

//_____________________________________________________________________________
G4String G4H3DummyManager::GetH3XAxisTitle(G4int /*id*/) const 
{
  ExceptionForHistograms("GetH3XAxisTitle");
  return "";
} 

//_____________________________________________________________________________
G4String G4H3DummyManager::GetH3YAxisTitle(G4int /*id*/) const 
{
  ExceptionForHistograms("GetH3YAxisTitle");
  return "";
}  

//_____________________________________________________________________________
G4String G4H3DummyManager::GetH3ZAxisTitle(G4int /*id*/) const 
{
  ExceptionForHistograms("GetH3ZAxisTitle");
  return "";
}  

//_____________________________________________________________________________
G4bool G4H3DummyManager::WriteOnAscii(std::ofstream& /*output*/)
{
  ExceptionForHistograms("WriteOnAscii");
  return false;
} 

