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
// $Id: G4H1DummyManager.cc 70604 2013-06-03 11:27:06Z ihrivnac $

// Author: Ivana Hrivnacova, 18/06/2013  (ivana@ipno.in2p3.fr)

#include "G4H1DummyManager.hh"

//_____________________________________________________________________________
G4H1DummyManager::G4H1DummyManager(const G4AnalysisManagerState& state)
 : G4VH1Manager(state)
{
}

//_____________________________________________________________________________
G4H1DummyManager::~G4H1DummyManager()
{  
}

// 
// protected methods
//

//_____________________________________________________________________________
G4int G4H1DummyManager::CreateH1(const G4String& /*name*/, 
                               const G4String& /*title*/, 
                               G4int /*nbins*/, 
                               G4double /*xmin*/, G4double /*xmax*/,
                               const G4String& /*unitName*/, 
                               const G4String& /*fcnName*/,
                               const G4String& /*binScheme*/)
{
  ExceptionForHistograms("CreateH1");
  return 0;
}                                         

//_____________________________________________________________________________
G4int G4H1DummyManager::CreateH1(const G4String& /*name*/, 
                               const G4String& /*title*/, 
                               const std::vector<G4double>& /*edges*/,
                               const G4String& /*unitName*/, 
                               const G4String& /*fcnName*/)
{
  ExceptionForHistograms("CreateH1");
  return 0;
}                                         

//_____________________________________________________________________________
G4bool G4H1DummyManager::SetH1(G4int /*id*/,
                                G4int /*nbins*/, 
                                G4double /*xmin*/, G4double /*xmax*/,
                                const G4String& /*unitName*/, 
                                const G4String& /*fcnName*/,
                                const G4String& /*binScheme*/)
{                                
  ExceptionForHistograms("SetH1");
  return false;
}
  
//_____________________________________________________________________________
G4bool G4H1DummyManager::SetH1(G4int /*id*/,
                                const std::vector<G4double>& /*edges*/,
                                const G4String& /*unitName*/, 
                                const G4String& /*fcnName*/)
{                                
  ExceptionForHistograms("SetH1");
  return false;
}
  
//_____________________________________________________________________________
G4bool G4H1DummyManager::ScaleH1(G4int /*id*/, G4double /*factor*/)
{
  ExceptionForHistograms("ScaleH1");
  return false;
}
                                  
//_____________________________________________________________________________
G4bool G4H1DummyManager::FillH1(G4int /*id*/, 
                                    G4double /*value*/, G4double /*weight*/)
{
  G4ExceptionDescription description;
  description << "      " 
              << "Histograms are not supported." ;
  G4Exception("G4H1DummyManager::FillH1()",
            "Analysis_W007", JustWarning, description);
  return false;
}

//_____________________________________________________________________________
G4int G4H1DummyManager::GetNofH1s() const
{
  ExceptionForHistograms("GetNofH1s");
  return 0;
}

//_____________________________________________________________________________
G4int  G4H1DummyManager::GetH1Id(const G4String& /*name*/, G4bool /*warn*/) const
{
  ExceptionForHistograms("GetH1Id");
  return 0;
}


//_____________________________________________________________________________
G4int G4H1DummyManager::GetH1Nbins(G4int /*id*/) const
{
  ExceptionForHistograms("GetH1Nbins");
  return 0;
}
  
//_____________________________________________________________________________
G4double G4H1DummyManager::GetH1Xmin(G4int /*id*/) const
{
  ExceptionForHistograms("GetH1Xmin");
  return 0;
}
  
//_____________________________________________________________________________
G4double G4H1DummyManager::GetH1Xmax(G4int /*id*/) const
{
  ExceptionForHistograms("GetH1Xmax");
  return 0;
}
  
//_____________________________________________________________________________
G4double G4H1DummyManager::GetH1Width(G4int /*id*/) const
{
  ExceptionForHistograms("GetH1Xwidth");
  return 0;
}

//_____________________________________________________________________________
G4bool G4H1DummyManager::SetH1Title(G4int /*id*/, 
                                    const G4String& /*title*/)
{
  ExceptionForHistograms("SetH1Title");
  return false;
}

//_____________________________________________________________________________
G4bool G4H1DummyManager::SetH1XAxisTitle(G4int /*id*/, 
                                         const G4String& /*title*/)
{
  ExceptionForHistograms("SetH1XAxisTitle");
  return false;
}

//_____________________________________________________________________________
G4bool G4H1DummyManager::SetH1YAxisTitle(G4int /*id*/, 
                                         const G4String& /*title*/)
{
  ExceptionForHistograms("SetH1YAxisTitle");
  return false;
}

//_____________________________________________________________________________
G4String G4H1DummyManager::GetH1XAxisTitle(G4int /*id*/) const
{
  ExceptionForHistograms("GetH1XAxisTitle");
  return "";
}

//_____________________________________________________________________________
G4String G4H1DummyManager::GetH1Title(G4int /*id*/) const
{
  ExceptionForHistograms("GetH1Title");
  return "";
}
  
//_____________________________________________________________________________
G4String G4H1DummyManager::GetH1YAxisTitle(G4int /*id*/) const
{
  ExceptionForHistograms("GetH1YAxisTitle");
  return "";
}

//_____________________________________________________________________________
G4bool G4H1DummyManager::WriteOnAscii(std::ofstream& /*output*/)
{
  ExceptionForHistograms("GetH1YAxisTitle");
  return false;
}

