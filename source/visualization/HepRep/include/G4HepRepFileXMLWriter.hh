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
//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: G4HepRepFileXMLWriter.hh 66373 2012-12-18 09:41:34Z gcosmo $
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
#ifndef G4HepRepFileXMLWriter_hh
#define G4HepRepFileXMLWriter_hh

//#define G4HEPREPFILEDEBUG  // Comment this out to suppress debug code.

#include "globals.hh"
#include <fstream>

class G4HepRepFileXMLWriter
{
public:
  G4HepRepFileXMLWriter();

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
  std::ofstream fout;

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
