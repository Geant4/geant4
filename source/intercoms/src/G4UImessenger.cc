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
//
// $Id: G4UImessenger.cc 84281 2014-10-13 07:21:50Z gcosmo $
//

#include "G4UImessenger.hh"
#include "G4UImanager.hh"
#include "G4UIcommand.hh"
#include "G4UIdirectory.hh"
#include "G4UIcommandTree.hh"
#include "G4ios.hh"
#include <sstream>

G4UImessenger::G4UImessenger()
  : baseDir(NULL), baseDirName(""), commandsShouldBeInMaster(false)
{ 
}

G4UImessenger::G4UImessenger(const G4String& path, const G4String& dsc,
                             G4bool commandsToBeBroadcasted)
  : baseDir(NULL), baseDirName(""), commandsShouldBeInMaster(false)
{
  CreateDirectory(path, dsc, commandsToBeBroadcasted);
}

G4UImessenger::~G4UImessenger()
{
  if(baseDir) delete baseDir;
}

G4String G4UImessenger::GetCurrentValue(G4UIcommand*) 
{ 
  G4String nullString;
  return nullString;
}

void G4UImessenger::SetNewValue(G4UIcommand*,G4String) 
{ ; }

G4bool G4UImessenger::operator == (const G4UImessenger& messenger) const {
  return this == &messenger;
}

G4String G4UImessenger::ItoS(G4int i)
{
  std::ostringstream os;
  os << i;
  return G4String(os.str());
}

G4String G4UImessenger::DtoS(G4double a)
{
  std::ostringstream os;
  os << a;
  return G4String(os.str());
}

G4String G4UImessenger::BtoS(G4bool b)
{
  G4String vl = "0";
  if(b) vl = "true";
  return vl;
}

G4int G4UImessenger::StoI(G4String str)
{
  G4int vl;
  const char* t = str;
  std::istringstream is(t);
  is >> vl;
  return vl;
}

G4double G4UImessenger::StoD(G4String str)
{
  G4double vl;
  const char* t = str;
  std::istringstream is(t);
  is >> vl;
  return vl;
}

G4bool G4UImessenger::StoB(G4String str)
{
  G4String v = str;
  v.toUpper();
  G4bool vl = false;
  if( v=="Y" || v=="YES" || v=="1" || v=="T" || v=="TRUE" )
  { vl = true; }
  return vl;
}


void G4UImessenger::AddUIcommand(G4UIcommand * newCommand)
{
  G4cerr << "Warning : Old style definition of G4UIcommand <" 
         << newCommand->GetCommandPath() << ">." << G4endl;
}

void G4UImessenger::CreateDirectory(const G4String& path, const G4String& dsc,
                             G4bool commandsToBeBroadcasted)
{
  G4UImanager* ui = G4UImanager::GetUIpointer();

  G4String fullpath = path;
  if(fullpath(fullpath.length()-1) != '/') fullpath.append("/");

  G4UIcommandTree* tree= ui-> GetTree()-> FindCommandTree(fullpath.c_str());
  if (tree) {
    baseDirName = tree-> GetPathName();
  } else {
    baseDir = new G4UIdirectory(fullpath.c_str(),commandsToBeBroadcasted);
    baseDirName = fullpath;
    baseDir-> SetGuidance(dsc.c_str());
  }
}
