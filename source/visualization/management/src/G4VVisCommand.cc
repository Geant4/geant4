// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VVisCommand.cc,v 1.8 2001-02-06 23:36:55 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

// Base class for visualization commands - John Allison  9th August 1998
// It is really a messenger - we have one command per messenger.

#include "G4VVisCommand.hh"

#include "G4UnitsTable.hh"
#include "g4std/strstream"

G4VVisCommand::~G4VVisCommand () {}

G4VisManager* G4VVisCommand::fpVisManager = 0;

G4std::vector <G4UIcommand*> G4VVisCommand::sceneNameCommands;

G4std::vector <G4UIcommand*> G4VVisCommand::sceneHandlerNameCommands;

G4std::vector <G4UIcommand*> G4VVisCommand::viewerNameCommands;

G4double G4VVisCommand::ValueOf(G4String unitName) {
   return G4UnitDefinition::GetValueOf(unitName);
}

G4String G4VVisCommand::ConvertToString(G4bool blValue)
{
  G4String vl = "false";
  if(blValue) vl = "true";
  return vl;
}

G4String G4VVisCommand::ConvertToString
(G4double x, G4double y, const char * unitName)
{
  G4double uv = ValueOf(unitName);
  
  char st[50];
  G4std::ostrstream os(st,50);
  os << x/uv << " " << y/uv << " " << unitName << G4std::ends;
  G4String vl = st;
  return vl;
}

G4String G4VVisCommand::ConvertToString(const G4ThreeVector& vec)
{
  char st[100];
  G4std::ostrstream os(st,100);
  os << vec.x() << " " << vec.y() << " " << vec.z() << G4std::ends;
  G4String vl = st;
  return vl;
}

G4bool G4VVisCommand::GetNewBoolValue(const G4String& paramString)
{
  G4String v = paramString;
  v.toUpper();
  G4bool vl = false;
  if( v=="Y" || v=="YES" || v=="1" || v=="T" || v=="TRUE" )
  { vl = true; }
  return vl;
}

G4int G4VVisCommand::GetNewIntValue(const G4String& paramString)
{
  G4int vl;
  const char* t = paramString;
  G4std::istrstream is((char*)t);
  is >> vl;
  return vl;
}

G4double G4VVisCommand::GetNewDoubleValue(const G4String& paramString)
{
  G4double vl;
  char unts[30];
  
  const char* t = paramString;
  G4std::istrstream is((char*)t);
  is >> vl >> unts;
  G4String unt = unts;
  
  return (vl*ValueOf(unt));
}

G4ThreeVector G4VVisCommand::GetNew3VectorValue(const G4String& paramString)
{
  G4double vx;
  G4double vy;
  G4double vz;
  const char* t = paramString;
  G4std::istrstream is((char*)t);
  is >> vx >> vy >> vz;
  return G4ThreeVector(vx,vy,vz);
}

void G4VVisCommand::GetNewDoublePairValue(const G4String& paramString,
					  G4double& xval,
					  G4double& yval)
{
  G4double x, y;
  char unts[30];
  
  const char* t = paramString;
  G4std::istrstream is((char*)t);
  is >> x >> y >> unts;
  G4String unt = unts;

  xval = x*ValueOf(unt);
  yval = y*ValueOf(unt);

  return;
}
