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
// 	$Id: HepRepXMLWriter.cc,v 1.4 2002-01-14 22:31:28 perl Exp $
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
}

void HepRepXMLWriter::addType(const char* name)
{
  if (fout.good())
  {
    endType();
    fout << "  <heprep:type version=\"null\" name=\"" << name << "\">"
	 << G4endl;
    inType = true;
  } else {
    G4cout << "HepRepXMLWriter:addType No file is currently open." << G4endl;
  }
}

void HepRepXMLWriter::addInstance()
{
  if (fout.good())
  {
    if (inType)
    {
      endInstance();
      fout << "    <heprep:instance>" << G4endl;
      inInstance = true;
    } else {
      G4cout << "HepRepXMLWriter:addInstance No HepRep Type is currently open"
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
    if (inInstance)
    {
      endPrimitive();
      fout << "      <heprep:primitive>" << G4endl;
      inPrimitive = true;
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
      fout << "        <heprep:point x=\"" << x << "\" y=\"" << y << "\" z=\"" << z << "\">" << G4endl;
      inPoint = true;
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
    fout << "  <heprep:attdef extra=\"" << extra << "\" name=\"" << name << "\" type=\"" << type << "\"" << G4endl;
    indent();
    fout << "    desc=\"" << desc << "\"/>" << G4endl;
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
    fout << "    value=\"" << value << "\"/>" << G4endl;
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
    indent();
    int redness = int(value1*255.);
    int greenness = int(value2*255.);
    int blueness = int(value3*255.);
    fout << "  <heprep:attvalue showLabel=\"NONE\" name=\"" << name << "\"" << G4endl;
    indent();
    fout << "    value=\"" << redness << "," << greenness << "," << blueness << "\"/>" << G4endl;
  } else {
    G4cout << "HepRepXMLWriter:addAttValue No file is currently open" << G4endl;
  }
}

void HepRepXMLWriter::open(const char* filespec)
{
  if (isOpen)
    fout.close();

  fout.open(filespec);

  if (fout.good())
  {
    fout << "<?xml version=\"1.0\" ?>" << G4endl;
    fout << "<heprep:heprep xmlns:heprep=\"http://www.freehep.org/HepRep\"" << G4endl;
    fout << "  xmlns:xsi=\"http://www.w3.org/1999/XMLSchema-instance\" xsi:schemaLocation=\"HepRep.xsd\">" << G4endl;

    isOpen = true;
    inType = false;
    inInstance = false;
    inPrimitive = false;
    inPoint = false;
  } else {
    G4cout << "HepRepXMLWriter:open Unable to write to file " << filespec << G4endl;
  }
}

void HepRepXMLWriter::close()
{
  if (fout.good())
  {
    endType();
    fout << "</heprep:heprep>" << G4endl;
    fout.close( );
  } else {
    G4cout << "HepRepXMLWriter:close No file is currently open" << G4endl;
  }
}

void HepRepXMLWriter::endType()
{
  if (inType)
  {
    endInstance();
    fout << "  </heprep:type>" << G4endl;
    inType = false;
  }
}

void HepRepXMLWriter::endInstance()
{
  if (inInstance)
  {
    endPrimitive();
    fout << "    </heprep:instance>" << G4endl;
    inInstance = false;
  }
}

void HepRepXMLWriter::endPrimitive()
{
  if (inPrimitive)
  {
    endPoint();
    fout << "      </heprep:primitive>" << G4endl;
    inPrimitive = false;
  }
}

void HepRepXMLWriter::endPoint()
{
  if (inPoint)
  {
    fout << "        </heprep:point>" << G4endl;
    inPoint = false;
  }
}

void HepRepXMLWriter::indent()
{
  if (fout.good())
  {
    if (inType)
      fout << "  ";
    if (inInstance)
      fout << "  ";
    if (inPrimitive)
      fout << "  ";
    if (inPoint)
      fout << "  ";
  }
}
