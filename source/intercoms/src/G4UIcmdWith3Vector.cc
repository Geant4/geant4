// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4UIcmdWith3Vector.cc,v 1.1 1999-01-07 16:09:25 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//

#include "G4UIcmdWith3Vector.hh"
#ifdef WIN32
#  include <Strstrea.h>
#else
#  include <strstream.h>
#endif

G4UIcmdWith3Vector::G4UIcmdWith3Vector
(const char * theCommandPath,G4UImessenger * theMessenger)
:G4UIcommand(theCommandPath,theMessenger)
{
  G4UIparameter * dblParamX = new G4UIparameter('d');
  SetParameter(dblParamX);
  G4UIparameter * dblParamY = new G4UIparameter('d');
  SetParameter(dblParamY);
  G4UIparameter * dblParamZ = new G4UIparameter('d');
  SetParameter(dblParamZ);
}

G4ThreeVector G4UIcmdWith3Vector::GetNew3VectorValue(G4String paramString)
{
  G4double vx;
  G4double vy;
  G4double vz;
  const char* t = paramString;
  istrstream is((char*)t);
  is >> vx >> vy >> vz;
  return G4ThreeVector(vx,vy,vz);
}

G4String G4UIcmdWith3Vector::ConvertToString(G4ThreeVector vec)
{
  char st[100];
  ostrstream os(st,100);
  os << vec.x() << " " << vec.y() << " " << vec.z() << '\0';
  G4String vl = st;
  return vl;
}

void G4UIcmdWith3Vector::SetParameterName
(const char * theNameX,const char * theNameY,const char * theNameZ,
G4bool omittable,G4bool currentAsDefault)
{
  G4UIparameter * theParamX = GetParameter(0);
  theParamX->SetParameterName(theNameX);
  theParamX->SetOmittable(omittable);
  theParamX->SetCurrentAsDefault(currentAsDefault);
  G4UIparameter * theParamY = GetParameter(1);
  theParamY->SetParameterName(theNameY);
  theParamY->SetOmittable(omittable);
  theParamY->SetCurrentAsDefault(currentAsDefault);
  G4UIparameter * theParamZ = GetParameter(2);
  theParamZ->SetParameterName(theNameZ);
  theParamZ->SetOmittable(omittable);
  theParamZ->SetCurrentAsDefault(currentAsDefault);
}

void G4UIcmdWith3Vector::SetDefaultValue(G4ThreeVector vec)
{
  G4UIparameter * theParamX = GetParameter(0);
  theParamX->SetDefaultValue(vec.x());
  G4UIparameter * theParamY = GetParameter(1);
  theParamY->SetDefaultValue(vec.y());
  G4UIparameter * theParamZ = GetParameter(2);
  theParamZ->SetDefaultValue(vec.z());
}

