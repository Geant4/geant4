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
//
// $Id: G4UImessenger.cc,v 1.6 2005/10/26 06:10:22 kmura Exp $
// GEANT4 tag $Name: geant4-08-00 $
//

#include "G4UImessenger.hh"
#include "G4UImanager.hh"
#include "G4UIcommand.hh"
#include "G4ios.hh"
#include <sstream>


G4UImessenger::G4UImessenger() { }

G4UImessenger::~G4UImessenger() { }

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

G4int G4UImessenger::StoI(G4String s)
{
  G4int vl;
  const char* t = s;
  std::istringstream is(t);
  is >> vl;
  return vl;
}

G4double G4UImessenger::StoD(G4String s)
{
  G4double vl;
  const char* t = s;
  std::istringstream is(t);
  is >> vl;
  return vl;
}

G4bool G4UImessenger::StoB(G4String s)
{
  G4String v = s;
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

