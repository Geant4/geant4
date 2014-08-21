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

#include "G4P1DummyManager.hh"

//
// Constructors, destructor
//

//_____________________________________________________________________________
G4P1DummyManager::G4P1DummyManager(const G4AnalysisManagerState& state)
 : G4VP1Manager(state)
{
}

//_____________________________________________________________________________
G4P1DummyManager::~G4P1DummyManager()
{  
}

// 
// protected methods
//

//_____________________________________________________________________________
G4int G4P1DummyManager::CreateP1(const G4String& /*name*/,  
                               const G4String& /*title*/,
                               G4int /*nbins*/, 
                               G4double /*xmin*/, G4double /*xmax*/,
                               G4double /*ymin*/, G4double /*ymax*/,
                               const G4String& /*xunitName*/, 
                               const G4String& /*yunitName*/, 
                               const G4String& /*xfcnName*/,
                               const G4String& /*yfcnName*/,
                               const G4String& /*xbinScheme*/)
{
  ExceptionForProfiles("CreateP1");
  return 0;
}                                         

//_____________________________________________________________________________
G4int G4P1DummyManager::CreateP1(const G4String& /*name*/,  
                                 const G4String& /*title*/,
                                 const std::vector<G4double>& /*edges*/,
                                 G4double /*ymin*/, G4double /*ymax*/,
                                 const G4String& /*xunitName*/, 
                                 const G4String& /*yunitName*/, 
                                 const G4String& /*xfcnName*/,
                                 const G4String& /*yfcnName*/)
{
  ExceptionForProfiles("CreateP1");
  return 0;
}                                         

//_____________________________________________________________________________
G4bool G4P1DummyManager::SetP1(G4int /*id*/,
                               G4int /*nbins*/, 
                               G4double /*xmin*/, G4double /*xmax*/,
                               G4double /*ymin*/, G4double /*ymax*/,
                               const G4String& /*xunitName*/, 
                               const G4String& /*yunitName*/, 
                               const G4String& /*xfcnName*/,
                               const G4String& /*yfcnName*/,
                               const G4String& /*xbinScheme*/)
{                                
  ExceptionForProfiles("SetP1");
  return false;
}

//_____________________________________________________________________________
G4bool G4P1DummyManager::SetP1(G4int /*id*/,
                               const std::vector<G4double>& /*edges*/,
                               G4double /*ymin*/, G4double /*ymax*/,
                               const G4String& /*xunitName*/, 
                               const G4String& /*yunitName*/, 
                               const G4String& /*xfcnName*/,
                               const G4String& /*yfcnName*/)
{
  ExceptionForProfiles("SetP1");
  return false;
}
                           
  
//_____________________________________________________________________________
G4bool G4P1DummyManager::ScaleP1(G4int /*id*/, G4double /*factor*/)
{
  ExceptionForProfiles("ScaleP1");
  return false;
}  

//_____________________________________________________________________________
G4bool G4P1DummyManager::FillP1(G4int /*id*/, 
                                G4double /*xvalue*/, G4double /*yvalue*/,
                                G4double /*weight*/)
{
  ExceptionForHistograms("FillP1");
  return false;
}

//_____________________________________________________________________________
G4int  G4P1DummyManager::GetP1Id(const G4String& /*name*/,  
                                 G4bool /*warn*/) const
{
  ExceptionForProfiles("GetP1Id");
  return 0;
}  

//_____________________________________________________________________________
G4int G4P1DummyManager::GetP1Nbins(G4int /*id*/) const
{
  ExceptionForProfiles("GetP1Nbins");
  return 0;
}  

//_____________________________________________________________________________
G4double G4P1DummyManager::GetP1Xmin(G4int /*id*/) const
{
// Returns xmin value with applied unit and histogram function

  ExceptionForProfiles("GetP1Xmin");
  return 0;
}  

//_____________________________________________________________________________
G4double G4P1DummyManager::GetP1Xmax(G4int /*id*/) const
{
  ExceptionForProfiles("GetP1Xmax");
  return 0;
}  

//_____________________________________________________________________________
G4double G4P1DummyManager::GetP1XWidth(G4int /*id*/) const
{
  ExceptionForProfiles("GetP1XWidth");
  return 0;
}  

//_____________________________________________________________________________
G4double G4P1DummyManager::GetP1Ymin(G4int /*id*/) const
{
// Returns xmin value with applied unit and histogram function

  ExceptionForProfiles("GetP1Ymin");
  return 0;
}  

//_____________________________________________________________________________
G4double G4P1DummyManager::GetP1Ymax(G4int /*id*/) const
{
  ExceptionForProfiles("GetP1Ymax");
  return 0;
}  

//_____________________________________________________________________________
G4bool G4P1DummyManager::SetP1Title(G4int /*id*/, const G4String& /*title*/)
{
  ExceptionForProfiles("SetP1Title");
  return false;
}  

//_____________________________________________________________________________
G4bool G4P1DummyManager::SetP1XAxisTitle(G4int /*id*/, 
                                         const G4String& /*title*/)
{
  ExceptionForProfiles("SetP1XAxisTitle");
  return false;
}  

//_____________________________________________________________________________
G4bool G4P1DummyManager::SetP1YAxisTitle(G4int /*id*/,
                                         const G4String& /*title*/)
{
  ExceptionForProfiles("SetP1YAxisTitle");
  return false;
}  

//_____________________________________________________________________________
G4String G4P1DummyManager::GetP1Title(G4int /*id*/) const
{
  ExceptionForProfiles("GetP1Title");
  return "";
}  


//_____________________________________________________________________________
G4String G4P1DummyManager::GetP1XAxisTitle(G4int /*id*/) const 
{
  ExceptionForProfiles("GetP1XAxisTitle");
  return "";
}  

//_____________________________________________________________________________
G4String G4P1DummyManager::GetP1YAxisTitle(G4int /*id*/) const 
{
  ExceptionForProfiles("GetP1YAxisTitle");
  return "";
}  
