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
// 	$Id: HepRepXMLWriter.cc,v 1.7 2002-02-02 04:00:30 perl Exp $
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

#include "HepRepXMLWriter.hh"

#include "G4ios.hh"

HepRepXMLWriter::HepRepXMLWriter()
{
  isOpen = false;
  init();
}

void HepRepXMLWriter::init()
{
  typeDepth = -1;

  int i = -1;
  while (i++<49) {
    prevTypeName[i] = new char[1];
    strcpy(prevTypeName[i],"");

    inType[i] = false;
    inInstance[i] = false;
  }

  inPrimitive = false;
  inPoint = false;
}

void HepRepXMLWriter::addType(const char* name,int newTypeDepth)
{
  if (fout.good())
  {
    // Flatten structure if it exceeds maximum allowed typeDepth of 49.
    if (newTypeDepth > 49)
      newTypeDepth = 49;

    if (newTypeDepth < 0)
      newTypeDepth = 0;

    // If moving closer to the root, close previously open types.
    while (newTypeDepth<typeDepth)
      endType();

    // Close any remaining primitives of the current instance.
    endPrimitive();

    // If this is a new type name for the current depth, declare the
    // new Type.  Otherwise, it is just another Instance of the current Type.
    if (strcmp(name,prevTypeName[newTypeDepth])!=0)
    {
      if (inType[newTypeDepth])
	endType();

      prevTypeName[newTypeDepth] = new char[strlen(name)+1];
      strcpy(prevTypeName[newTypeDepth],name);

      inType[newTypeDepth] = true;
      indent();
      fout << "<heprep:type version=\"null\" name=\"" << name << "\">"
	 << G4endl;

      typeDepth = newTypeDepth;
    }
  } else {
    G4cout << "HepRepXMLWriter:addType No file is currently open." << G4endl;
  }
}

void HepRepXMLWriter::addInstance()
{
  if (fout.good())
  {
    if (inType[typeDepth])
    {
      endInstance();
      inInstance[typeDepth] = true;
      indent();
      fout << "<heprep:instance>" << G4endl;
    } else {
      G4cout <<
	"HepRepXMLWriter:addInstance No HepRep Type is currently open"
	     << G4endl;
    }
  } else {
    G4cout << "HepRepXMLWriter:addInstance No file is currently open"
	   << G4endl;
  }
}

void HepRepXMLWriter::addPrimitive()
{
  if (fout.good())
  {
    if (inInstance[typeDepth])
    {
      endPrimitive();
      inPrimitive = true;
      indent();
      fout << "<heprep:primitive>" << G4endl;
    } else {
      G4cout <<
	"HepRepXMLWriter:addPrimitive No HepRep Instance is currently open"
	     << G4endl;
    }
  } else {
    G4cout << "HepRepXMLWriter:addPrimitive No file is currently open"
	   << G4endl;
  }
}

void HepRepXMLWriter::addPoint(double x, double y, double z)
{
  if (fout.good())
  {
    if (inPrimitive)
    {
      endPoint();
      inPoint = true;
      indent();
      fout << "<heprep:point x=\"" << x << "\" y=\"" << y << "\" z=\"" << z << "\">" << G4endl;
    } else {
      G4cout <<
	"HepRepXMLWriter:addPoint No HepRep Primitive is currently open"
	     << G4endl;
    }
  } else {
    G4cout << "HepRepXMLWriter:addPoint No file is currently open" << G4endl;
  }
}

void HepRepXMLWriter::addAttDef(const char* name,
	       const char* desc,
	       const char* type,
	       const char* extra)
{
  if (fout.good())
  {
    indent();
    fout << "<heprep:attdef extra=\"" << extra << "\" name=\"" << name << "\" type=\"" << type << "\"" << G4endl;
    indent();
    fout << "  desc=\"" << desc << "\"/>" << G4endl;
  } else {
    G4cout << "HepRepXMLWriter:addAttDef No file is currently open" << G4endl;
  }
}

// Four methods to fill attValues
void HepRepXMLWriter::addAttValue (const char* name,
		  const char* value)
{
  if (fout.good())
  {
    indent();
    fout << "  <heprep:attvalue showLabel=\"NONE\" name=\"" << name << "\"" << G4endl;
    indent();
    fout << "    value=\"" << value << "\"/>" << G4endl;
  } else {
    G4cout << "HepRepXMLWriter:addAttValue No file is currently open" << G4endl;
  }
}

void HepRepXMLWriter::addAttValue (const char* name,
		  double value)
{
  if (fout.good())
  {
    indent();
    fout << "  <heprep:attvalue showLabel=\"NONE\" name=\"" << name << "\"" << G4endl;
    indent();
    fout << "    value=\"" << value << "\"/>" << G4endl;
  } else {
    G4cout << "HepRepXMLWriter:addAttValue No file is currently open" << G4endl;
  }
}

void HepRepXMLWriter::addAttValue (const char* name,
		  int value)
{
  if (fout.good())
  {
    indent();
    fout << "  <heprep:attvalue showLabel=\"NONE\" name=\"" << name << "\"" << G4endl;
    indent();
    fout << "    value=\"" << value << "\"/>" << G4endl;
  } else {
    G4cout << "HepRepXMLWriter:addAttValue No file is currently open" << G4endl;
  }
}

void HepRepXMLWriter::addAttValue (const char* name,
		  bool value)
{
  if (fout.good())
  {
    indent();
    fout << "  <heprep:attvalue showLabel=\"NONE\" name=\"" << name << "\"" << G4endl;
    indent();
    if (value)
      fout << "    value=\"True\"/>" << G4endl;
    else
      fout << "    value=\"False\"/>" << G4endl;
  } else {
    G4cout << "HepRepXMLWriter:addAttValue No file is currently open" << G4endl;
  }
}

void HepRepXMLWriter::addAttValue (const char* name,
				   double value1,
				   double value2,
				   double value3)
{
  if (fout.good())
  {
    int redness = int(value1*255.);
    int greenness = int(value2*255.);
    int blueness = int(value3*255.);
    indent();
    fout << "  <heprep:attvalue showLabel=\"NONE\" name=\"" << name << "\"" << G4endl;
    indent();
    fout << "    value=\"" << redness << "," << greenness << "," << blueness << "\"/>" << G4endl;
  } else {
    G4cout << "HepRepXMLWriter:addAttValue No file is currently open" << G4endl;
  }
}

void HepRepXMLWriter::open(const char* fileSpec)
{
  if (isOpen)
    close();
  
  fout.open(fileSpec);
    
  if (fout.good()) {
    fout << "<?xml version=\"1.0\" ?>" << G4endl;
    fout << "<heprep:heprep xmlns:heprep=\"http://www.freehep.org/HepRep\"" << G4endl;
    fout << "  xmlns:xsi=\"http://www.w3.org/1999/XMLSchema-instance\" xsi:schemaLocation=\"HepRep.xsd\">" << G4endl;
    
    isOpen = true;
    init();
  } else {
    G4cout << "HepRepXMLWriter:open Unable to write to file " << fileSpec << G4endl;
  }
}

void HepRepXMLWriter::close()
{
  // Close any remaining open Types
  endTypes();

  if (fout.good()) {      
    fout << "</heprep:heprep>" << G4endl;
    fout.close( );
    isOpen = false;
  } else {
    G4cout << "HepRepXMLWriter:close No file is currently open" << G4endl;
  }
}

void HepRepXMLWriter::endTypes()
{
  // Close any remaining open Types
    while(typeDepth>-1)
      endType();
}

void HepRepXMLWriter::endType()
{
  endInstance();
  indent();
  fout << "</heprep:type>" << G4endl;
  inType[typeDepth] = false;
  delete [] prevTypeName[typeDepth];
  prevTypeName[typeDepth] = new char[1];
  strcpy(prevTypeName[typeDepth],"");
  typeDepth--;
}

void HepRepXMLWriter::endInstance()
{
  if (inInstance[typeDepth])
  {
    endPrimitive();
    indent();
    fout << "</heprep:instance>" << G4endl;
    inInstance[typeDepth] = false;
  }
}

void HepRepXMLWriter::endPrimitive()
{
  if (inPrimitive)
  {
    endPoint();
    indent();
    fout << "</heprep:primitive>" << G4endl;
    inPrimitive = false;
  }
}

void HepRepXMLWriter::endPoint()
{
  if (inPoint)
  {
    indent();
    fout << "</heprep:point>" << G4endl;
    inPoint = false;
  }
}

void HepRepXMLWriter::indent()
{
  if (fout.good())
  {
    int i = 0;
    while (inType[i] && i<12) {
      fout << "  ";
      if (inInstance[i])
        fout << "  ";
      i++;
    }

    if (inPrimitive)
      fout << "  ";
    if (inPoint)
      fout << "  ";
  }
}
