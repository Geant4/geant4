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
// $Id: G4Pstring.cc,v 1.3 2002-08-29 15:30:51 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// GEANT 4 class source file
//
// G4Pstring.cc
//
// ----------------------------------------------------------------------

#include "G4Pstring.hh"

G4String str(const int &i)
{
  char s[200];
  sprintf(s,"%d",i);
  return s;
}

G4String str(const double &d)
{
  char s[200];
  sprintf(s,"%f",d);
  return s;
}

G4String str(const G4ThreeVector &v)
{
  G4String out = "(";
  out += str(v.x());
  out += ", "; 
  out += str(v.y()); 
  out += ", ";
  out += str(v.z());
  out += ")";
  return out;
}


