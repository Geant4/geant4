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
// $Id: G4H2DummyManager.cc 70604 2013-06-03 11:27:06Z ihrivnac $

// Author: Ivana Hrivnacova, 18/06/2013  (ivana@ipno.in2p3.fr)

#include "G4H2DummyManager.hh"

//_____________________________________________________________________________
G4H2DummyManager::G4H2DummyManager(const G4AnalysisManagerState& state)
 : G4VH2Manager(state)
{
}

//_____________________________________________________________________________
G4H2DummyManager::~G4H2DummyManager()
{  
}

// 
// protected methods
//

//_____________________________________________________________________________
G4int G4H2DummyManager::CreateH2(const G4String& /*name*/, 
                               const G4String& /*title*/, 
                               G4int /*nxbins*/, 
                               G4double /*xmin*/, G4double /*xmax*/,
                               G4int /*nybins*/, 
                               G4double /*ymin*/, G4double /*ymax*/,
                               const G4String& /*xunitName*/, 
                               const G4String& /*yunitName*/, 
                               const G4String& /*xfcnName*/,
                               const G4String& /*yfcnName*/,
                               const G4String& /*xbinScheme*/,
                               const G4String& /*ybinScheme*/)
{
  ExceptionForHistograms("CreateH2");
  return 0;
}                                         

//_____________________________________________________________________________
G4int G4H2DummyManager::CreateH2(const G4String& /*name*/, 
                               const G4String& /*title*/, 
                               const std::vector<G4double>& /*xedges*/,
                               const std::vector<G4double>& /*yedges*/,
                               const G4String& /*xunitName*/, 
                               const G4String& /*yunitName*/, 
                               const G4String& /*xfcnName*/,
                               const G4String& /*yfcnName*/)
{
  ExceptionForHistograms("CreateH2");
  return 0;
}                                         

//_____________________________________________________________________________
G4bool G4H2DummyManager::SetH2(G4int /*id*/,
                                G4int /*nxbins*/, 
                                G4double /*xmin*/, G4double /*xmax*/, 
                                G4int /*nybins*/, 
                                G4double /*ymin*/, G4double /*ymax*/,
                                const G4String& /*xunitName*/, 
                                const G4String& /*yunitName*/, 
                                const G4String& /*xfcnName*/,
                                const G4String& /*yfcnName*/,
                                const G4String& /*xbinScheme*/,
                                const G4String& /*ybinScheme*/)
{                                
  ExceptionForHistograms("SetH2");
  return false;
}

//_____________________________________________________________________________
G4bool G4H2DummyManager::SetH2(G4int /*id*/,
                                const std::vector<G4double>& /*xedges*/,
                                const std::vector<G4double>& /*yedges*/,
                                const G4String& /*xunitName*/, 
                                const G4String& /*yunitName*/, 
                                const G4String& /*xfcnName*/,
                                const G4String& /*yfcnName*/)
{                                
  ExceptionForHistograms("SetH2");
  return false;
}

//_____________________________________________________________________________
G4bool G4H2DummyManager::ScaleH2(G4int /*id*/, G4double /*factor*/)
{
  ExceptionForHistograms("ScaleH2");
  return false;
}

//_____________________________________________________________________________
G4bool G4H2DummyManager::FillH2(G4int /*id*/, 
                                G4double /*xvalue*/, G4double /*yvalue*/,
                                G4double /*weight*/)
{
  G4ExceptionDescription description;
  description << "      " 
              << "Histograms are not supported." ;
  G4Exception("G4H2DummyManager::FillH2()",
            "Analysis_W007", JustWarning, description);
  return false;
}

//_____________________________________________________________________________
G4int G4H2DummyManager::GetNofH2s() const
{
  ExceptionForHistograms("GetNofH2s");
  return 0;
}

//_____________________________________________________________________________
G4int  G4H2DummyManager::GetH2Id(const G4String& /*name*/, G4bool /*warn*/) const
{
  ExceptionForHistograms("GetH2Id");
  return 0;
}


//_____________________________________________________________________________
G4int G4H2DummyManager::GetH2Nxbins(G4int /*id*/) const
{
  ExceptionForHistograms("GetH2NXbins");
  return 0;
}
  
//_____________________________________________________________________________
G4double G4H2DummyManager::GetH2Xmin(G4int /*id*/) const
{
  ExceptionForHistograms("GetH2Xmin");
  return 0;
}
  
//_____________________________________________________________________________
G4double G4H2DummyManager::GetH2Xmax(G4int /*id*/) const
{
  ExceptionForHistograms("GetH2Xmin");
  return 0;
}
  
//_____________________________________________________________________________
G4double G4H2DummyManager::GetH2XWidth(G4int /*id*/) const
{
  ExceptionForHistograms("GetH2XWidth");
  return 0;
}
  
//_____________________________________________________________________________
G4int G4H2DummyManager::GetH2Nybins(G4int /*id*/) const
{
  ExceptionForHistograms("GetH2NYbins");
  return 0;
}
  
//_____________________________________________________________________________
G4double G4H2DummyManager::GetH2Ymin(G4int /*id*/) const
{
  ExceptionForHistograms("GetH2Ymin");
  return 0;
}
  
//_____________________________________________________________________________
G4double G4H2DummyManager::GetH2Ymax(G4int /*id*/) const
{
  ExceptionForHistograms("GetH2Ymax");
  return 0;
}
  
//_____________________________________________________________________________
G4double G4H2DummyManager::GetH2YWidth(G4int /*id*/) const
{
  ExceptionForHistograms("GetH2YWidth");
  return 0;
}

//_____________________________________________________________________________
G4bool G4H2DummyManager::SetH2Title(G4int /*id*/, 
                                    const G4String& /*title*/)
{
  ExceptionForHistograms("SetH2Title");
  return false;
}

//_____________________________________________________________________________
G4bool G4H2DummyManager::SetH2XAxisTitle(G4int /*id*/, 
                                         const G4String& /*title*/)
{
  ExceptionForHistograms("SetH2XAxisTitle");
  return false;
}

//_____________________________________________________________________________
G4bool G4H2DummyManager::SetH2YAxisTitle(G4int /*id*/, 
                                         const G4String& /*title*/)
{
  ExceptionForHistograms("SetH2YAxisTitle");
  return false;
}

//_____________________________________________________________________________
G4bool G4H2DummyManager::SetH2ZAxisTitle(G4int /*id*/, 
                                         const G4String& /*title*/)
{
  ExceptionForHistograms("SetH2ZAxisTitle");
  return false;
}

//_____________________________________________________________________________
G4String G4H2DummyManager::GetH2Title(G4int /*id*/) const
{
  ExceptionForHistograms("GetH2Title");
  return "";
}
  
//_____________________________________________________________________________
G4String G4H2DummyManager::GetH2XAxisTitle(G4int /*id*/) const
{
  ExceptionForHistograms("GetH2XAxisTitle");
  return "";
}

//_____________________________________________________________________________
G4String G4H2DummyManager::GetH2YAxisTitle(G4int /*id*/) const
{
  ExceptionForHistograms("GetH2YAxisTitle");
  return "";
}

//_____________________________________________________________________________
G4String G4H2DummyManager::GetH2ZAxisTitle(G4int /*id*/) const
{
  ExceptionForHistograms("GetH2ZAxisTitle");
  return "";
}
  
//_____________________________________________________________________________
G4bool G4H2DummyManager::WriteOnAscii(std::ofstream& /*output*/)
{
  ExceptionForHistograms("WriteOnAscii");
  return false;
}
