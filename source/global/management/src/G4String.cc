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
// $Id: G4String.cc,v 1.1 2008-12-08 10:07:22 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
//---------------------------------------------------------------
//  GEANT 4 class implementation file
//
//  G4String
//---------------------------------------------------------------

#include "G4String.hh"

G4String::G4String ()
{
}

G4String::~G4String ()
{
}

G4int G4String::strcasecompare(const char* s1, const char* s2) const
{
  char* buf1 = new char[strlen(s1)+1];
  char* buf2 = new char[strlen(s2)+1];

  for (str_size i=0; i<=strlen(s1); i++)
    { buf1[i] = tolower(char(s1[i])); }
  for (str_size j=0; j<=strlen(s2); j++)
    { buf2[j] = tolower(char(s2[j])); }

  G4int res = strcmp(buf1, buf2);
  delete [] buf1;
  delete [] buf2;
  return res;
}

G4String G4String::strip (G4int strip_Type, char c)
{
  G4String retVal = *this;
  if(length()==0) { return retVal; }
  str_size i=0;
  switch ( strip_Type )
  {
    case leading: 
    {
      for(i=0;i<length();i++)
	{ if (std_string::operator[](i) != c) { break; } }
      retVal = substr(i,length()-i);
    }
    break;
    case trailing:
    {
      G4int j=0;
      for(j=length()-1;j>=0;j--)
	{ if (std_string::operator[](j) != c) { break; } }
      retVal = substr(0,j+1);
    }
    break;
    case both:
    { 
      for(i=0;i<length();i++)
	{ if (std_string::operator[](i) != c) { break; } }
      G4String tmp(substr(i,length()-i));
      G4int k=0;
      for(k=tmp.length()-1;k>=0;k--)
	{ if (tmp.std_string::operator[](k) != c) { break; } }
      retVal = tmp.substr(0,k+1);
    }
    break;
    default:
    break;
  }
  return retVal;
}

void G4String::toLower ()
{
  for (str_size i=0; i<size();i++)
  {
    std_string::operator[](i) = tolower(char(std_string::operator[](i)));
    //at(i) = tolower(at(i)); 
  } 
}

void G4String::toUpper ()
{
  for (str_size i=0; i<size();i++)
  {
    std_string::operator[](i) = toupper(char(std_string::operator[](i)));
    //at(i) = toupper(at(i)); 
  }
}

unsigned int G4String::hash( caseCompare ) const
{
  const char*s=c_str();
  unsigned long h = 0;
  for ( ; *s; ++s)
  {
    h = 5*h + *s;
  }
  return str_size(h);
}

unsigned int G4String::stlhash() const
{
  const char*s=c_str();
  unsigned long h = 0;
  for ( ; *s; ++s)
  {
    h = 5*h + *s;
  }
  return str_size(h);
}
