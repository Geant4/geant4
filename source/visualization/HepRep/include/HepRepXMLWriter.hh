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
//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: HepRepXMLWriter.hh,v 1.10 2003-01-24 21:18:48 perl Exp $
//
// Description:
//	Create a HepRep XML File (HepRep version 1).
//
// Environment:
//	Software developed for the general High Energy Physics community.
//
// Author :
//       J. Perl                    Original Author
//
// Copyright Information:
//      Copyright (C) 2001          Stanford Linear Accelerator Center
//------------------------------------------------------------------------
#ifndef HepRepXMLWriter_hh
#define HepRepXMLWriter_hh

//#define G4HEPREPFILEDEBUG  // Comment this out to suppress debug code.

#include "globals.hh"
#include "g4std/fstream"

class HepRepXMLWriter
{
public:
  HepRepXMLWriter();

  void addType(const char* name, int newTypeDepth);
  void addInstance();
  void addPrimitive();
  void addPoint(double x, double y, double z);

  void addAttDef(const char* name,
		 const char* desc,
		 const char* type,
		 const char* extra);

  void addAttValue(const char* name,
		   const char* value);

  void addAttValue(const char* name,
		   double value);

  void addAttValue(const char* name,
		   int value);

  void addAttValue(const char* name,
		   bool value);

  void addAttValue(const char* name,
		   double value1,
		   double value2,
		   double value3);

  void open(const char* filespec);
  void close();

  void endTypes();

  bool isOpen;
  int typeDepth;
  bool inType[50];
  bool inInstance[50];
  char* prevTypeName[50];
  
private:
  G4std::ofstream fout;

  void init();

  bool inPrimitive;
  bool inPoint;

  void endType();
  void endInstance();
  void endPrimitive();
  void endPoint();

  void indent();
};
#endif
