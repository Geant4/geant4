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

#include "G4P2DummyManager.hh"

//_____________________________________________________________________________
G4P2DummyManager::G4P2DummyManager(const G4AnalysisManagerState& state)
 : G4VP2Manager(state)
{
}

//_____________________________________________________________________________
G4P2DummyManager::~G4P2DummyManager()
{  
}

// 
// protected methods
//

//_____________________________________________________________________________
G4int G4P2DummyManager::CreateP2(const G4String& /*name*/, 
                               const G4String& /*title*/, 
                               G4int /*nxbins*/, 
                               G4double /*xmin*/, G4double /*xmax*/,
                               G4int /*nybins*/, 
                               G4double /*ymin*/, G4double /*ymax*/,
                               G4double /*zmin*/, G4double /*zmax*/,
                               const G4String& /*xunitName*/, 
                               const G4String& /*yunitName*/, 
                               const G4String& /*zunitName*/, 
                               const G4String& /*xfcnName*/,
                               const G4String& /*yfcnName*/,
                               const G4String& /*zfcnName*/,
                               const G4String& /*xbinScheme*/,
                               const G4String& /*ybinScheme*/)
{
  ExceptionForHistograms("CreateP2");
  return 0;
}                                         

//_____________________________________________________________________________
G4int G4P2DummyManager::CreateP2(const G4String& /*name*/, 
                               const G4String& /*title*/,
                               const std::vector<G4double>& /*xedges*/,
                               const std::vector<G4double>& /*yedges*/,
                               G4double /*zmin*/, G4double /*zmax*/,
                               const G4String& /*xunitName*/, 
                               const G4String& /*yunitName*/, 
                               const G4String& /*zunitName*/, 
                               const G4String& /*xfcnName*/,
                               const G4String& /*yfcnName*/,
                               const G4String& /*zfcnName*/)
                               
{
  ExceptionForHistograms("CreateP2");
  return 0;
}                                         

//_____________________________________________________________________________
G4bool G4P2DummyManager::SetP2(G4int /*id*/,
                               G4int /*nxbins*/, 
                               G4double /*xmin*/, G4double /*xmax*/,
                               G4int /*nybins*/, 
                               G4double /*ymin*/, G4double /*ymax*/,
                               G4double /*zmin*/, G4double /*zmax*/,
                               const G4String& /*xunitName*/, 
                               const G4String& /*yunitName*/, 
                               const G4String& /*zunitName*/, 
                               const G4String& /*xfcnName*/,
                               const G4String& /*yfcnName*/,
                               const G4String& /*zfcnName*/,
                               const G4String& /*xbinScheme*/,
                               const G4String& /*ybinScheme*/)
{                                
  ExceptionForHistograms("SetP2");
  return false;
}
                                  
//_____________________________________________________________________________
G4bool G4P2DummyManager::SetP2(G4int /*id*/,
                               const std::vector<G4double>& /*xedges*/,
                               const std::vector<G4double>& /*yedges*/,
                               G4double /*zmin*/, G4double /*zmax*/,
                               const G4String& /*xunitName*/, 
                               const G4String& /*yunitName*/, 
                               const G4String& /*zunitName*/, 
                               const G4String& /*xfcnName*/,
                               const G4String& /*yfcnName*/,
                               const G4String& /*zfcnName*/)
{                                
  ExceptionForHistograms("SetP2");
  return false;
}
                                  
//_____________________________________________________________________________
G4bool G4P2DummyManager::ScaleP2(G4int /*id*/, G4double /*factor*/)
{
  ExceptionForHistograms("ScaleP2");
  return false;
}  
                           
//_____________________________________________________________________________
G4bool G4P2DummyManager::FillP2(G4int /*id*/, 
                                G4double /*xvalue*/, G4double /*yvalue*/,
                                G4double /*zvalue*/,
                                G4double /*weight*/)
{
  ExceptionForHistograms("FillP2");
  return false;
}

//_____________________________________________________________________________
G4int  G4P2DummyManager::GetP2Id(const G4String& /*name*/, G4bool /*warn*/) const
{
  ExceptionForHistograms("GetP2Id");
  return 0;
}  
                                      
//_____________________________________________________________________________
G4int G4P2DummyManager::GetP2Nxbins(G4int /*id*/) const
{
  ExceptionForHistograms("GetP2Nxbins");
  return 0;
}  

//_____________________________________________________________________________
G4double G4P2DummyManager::GetP2Xmin(G4int /*id*/) const
{
  ExceptionForHistograms("GetP2Xmin");
  return 0;
}  

//_____________________________________________________________________________
G4double G4P2DummyManager::GetP2Xmax(G4int /*id*/) const
{
  ExceptionForHistograms("GetP2Xmax");
  return 0;
}  

//_____________________________________________________________________________
G4double G4P2DummyManager::GetP2XWidth(G4int /*id*/) const
{
  ExceptionForHistograms("GetP2XWidth");
  return 0;
}  

//_____________________________________________________________________________
G4int G4P2DummyManager::GetP2Nybins(G4int /*id*/) const
{
  ExceptionForHistograms("GetP2Nybins");
  return 0;
}  

//_____________________________________________________________________________
G4double G4P2DummyManager::GetP2Ymin(G4int /*id*/) const
{
  ExceptionForHistograms("");
  return 0;
}  

//_____________________________________________________________________________
G4double G4P2DummyManager::GetP2Ymax(G4int /*id*/) const
{
  ExceptionForHistograms("GetP2Ymax");
  return 0;
}  

//_____________________________________________________________________________
G4double G4P2DummyManager::GetP2YWidth(G4int /*id*/) const
{
  ExceptionForHistograms("GetP2YWidth");
  return 0;
}  

//_____________________________________________________________________________
G4double G4P2DummyManager::GetP2Zmin(G4int /*id*/) const
{
  ExceptionForHistograms("GetP2Zmin");
  return 0;
}  

//_____________________________________________________________________________
G4double G4P2DummyManager::GetP2Zmax(G4int /*id*/) const
{
  ExceptionForHistograms("GetP2Zmax");
  return 0;
}  

//_____________________________________________________________________________
G4bool G4P2DummyManager::SetP2Title(G4int /*id*/, const G4String& /*title*/)
{
  ExceptionForHistograms("SetP2Title");
  return false;
}  

//_____________________________________________________________________________
G4bool G4P2DummyManager::SetP2XAxisTitle(G4int /*id*/, const G4String& /*title*/)
{
  ExceptionForHistograms("SetP2XAxisTitle");
  return false;
}  

//_____________________________________________________________________________
G4bool G4P2DummyManager::SetP2YAxisTitle(G4int /*id*/, const G4String& /*title*/)
{
  ExceptionForHistograms("SetP2YAxisTitle");
  return false;
}  

//_____________________________________________________________________________
G4bool G4P2DummyManager::SetP2ZAxisTitle(G4int /*id*/, const G4String& /*title*/)
{
  ExceptionForHistograms("SetP2ZAxisTitle");
  return false;
}  

//_____________________________________________________________________________
G4String G4P2DummyManager::GetP2Title(G4int /*id*/) const
{
  ExceptionForHistograms("GetP2Title");
  return "";
}  

//_____________________________________________________________________________
G4String G4P2DummyManager::GetP2XAxisTitle(G4int /*id*/) const 
{
  ExceptionForHistograms("GetP2XAxisTitle");
  return "";
} 

//_____________________________________________________________________________
G4String G4P2DummyManager::GetP2YAxisTitle(G4int /*id*/) const 
{
  ExceptionForHistograms("GetP2YAxisTitle");
  return "";
}  

//_____________________________________________________________________________
G4String G4P2DummyManager::GetP2ZAxisTitle(G4int /*id*/) const 
{
  ExceptionForHistograms("GetP2ZAxisTitle");
  return "";
}  

//_____________________________________________________________________________
G4bool G4P2DummyManager::WriteOnAscii(std::ofstream& /*output*/)
{
  ExceptionForHistograms("WriteOnAscii");
  return false;
} 

