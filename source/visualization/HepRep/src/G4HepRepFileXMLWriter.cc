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
//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: G4HepRepFileXMLWriter.cc 66373 2012-12-18 09:41:34Z gcosmo $
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

#include "G4HepRepFileXMLWriter.hh"

#include "G4HepRepMessenger.hh"
#include "G4ios.hh"

G4HepRepFileXMLWriter::G4HepRepFileXMLWriter()
{
  isOpen = false;
  init();
}

void G4HepRepFileXMLWriter::init()
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

void G4HepRepFileXMLWriter::addType(const char* name,int newTypeDepth)
{
  if (fout.good())
  {
    // Flatten structure if it exceeds maximum allowed typeDepth of 49.
    if (newTypeDepth > 49)
      newTypeDepth = 49;

    if (newTypeDepth < 0)
      newTypeDepth = 0;

    // Insert any layers that are missing from the hierarchy (protects against
    // callers that skip from, say, layer 1 to layer 3 with no layer 2).
    while (typeDepth < (newTypeDepth-1)) {
      addType("Layer Inserted by G4HepRepFileXMLWriter", typeDepth + 1);
      addInstance();
    }

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
#ifdef G4HEPREPFILEDEBUG
    G4cout << "G4HepRepFileXMLWriter:addType No file is currently open." << G4endl;
#endif
  }
}

void G4HepRepFileXMLWriter::addInstance()
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
#ifdef G4HEPREPFILEDEBUG
      G4cout << "G4HepRepFileXMLWriter:addInstance No HepRep Type is currently open" << G4endl;
#endif
    }
  } else {
#ifdef G4HEPREPFILEDEBUG
    G4cout << "G4HepRepFileXMLWriter:addInstance No file is currently open" << G4endl;
#endif
  }
}

void G4HepRepFileXMLWriter::addPrimitive()
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
#ifdef G4HEPREPFILEDEBUG
      G4cout << "G4HepRepFileXMLWriter:addPrimitive No HepRep Instance is currently open" << G4endl; 
#endif
    }
  } else {
#ifdef G4HEPREPFILEDEBUG
    G4cout << "G4HepRepFileXMLWriter:addPrimitive No file is currently open" << G4endl;
#endif
  }
}

void G4HepRepFileXMLWriter::addPoint(double x, double y, double z)
{
  if (fout.good())
  {
    if (inPrimitive)
    {
      endPoint();
      inPoint = true;
      indent();
		
		// Include scale and center values
		G4HepRepMessenger* messenger = G4HepRepMessenger::GetInstance();
		G4double scale = messenger->getScale();
		G4ThreeVector center = messenger->getCenter();
		G4double xNew = scale * ( x - center.x());
		G4double yNew = scale * ( y - center.y());
		G4double zNew = scale * ( z - center.z());
		
      fout << "<heprep:point x=\"" << xNew << "\" y=\"" << yNew << "\" z=\"" << zNew << "\">" << G4endl;
    } else {
#ifdef G4HEPREPFILEDEBUG
      G4cout <<	"G4HepRepFileXMLWriter:addPoint No HepRep Primitive is currently open" << G4endl;
#endif
    }
  } else {
#ifdef G4HEPREPFILEDEBUG
    G4cout << "G4HepRepFileXMLWriter:addPoint No file is currently open" << G4endl;
#endif
  }
}

void G4HepRepFileXMLWriter::addAttDef(const char* name,
	       const char* desc,
	       const char* type,
	       const char* extra)
{
  if (fout.good())
  {
    indent();
    fout << "  <heprep:attdef extra=\"" << extra << "\" name=\"" << name << "\" type=\"" << type << "\"" << G4endl;
    indent();
    fout << "  desc=\"" << desc << "\"/>" << G4endl;
  } else {
#ifdef G4HEPREPFILEDEBUG
    G4cout << "G4HepRepFileXMLWriter:addAttDef No file is currently open" << G4endl;
#endif
  }
}

// Four methods to fill attValues
void G4HepRepFileXMLWriter::addAttValue (const char* name,
		  const char* value)
{
  if (fout.good())
  {
    indent();
    fout << "  <heprep:attvalue showLabel=\"NONE\" name=\"" << name << "\"" << G4endl;
    indent();
    fout << "    value=\"" << value << "\"/>" << G4endl;
  } else {
#ifdef G4HEPREPFILEDEBUG
    G4cout << "G4HepRepFileXMLWriter:addAttValue No file is currently open" << G4endl;
#endif
  }
}

void G4HepRepFileXMLWriter::addAttValue (const char* name,
		  double value)
{
  if (fout.good())
  {
    indent();
    fout << "  <heprep:attvalue showLabel=\"NONE\" name=\"" << name << "\"" << G4endl;
    indent();
    fout << "    value=\"" << value << "\"/>" << G4endl;
  } else {
#ifdef G4HEPREPFILEDEBUG
    G4cout << "G4HepRepFileXMLWriter:addAttValue No file is currently open" << G4endl;
#endif
  }
}

void G4HepRepFileXMLWriter::addAttValue (const char* name,
		  int value)
{
  if (fout.good())
  {
    indent();
    fout << "  <heprep:attvalue showLabel=\"NONE\" name=\"" << name << "\"" << G4endl;
    indent();
    fout << "    value=\"" << value << "\"/>" << G4endl;
  } else {
#ifdef G4HEPREPFILEDEBUG
    G4cout << "G4HepRepFileXMLWriter:addAttValue No file is currently open" << G4endl;
#endif
  }
}

void G4HepRepFileXMLWriter::addAttValue (const char* name,
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
#ifdef G4HEPREPFILEDEBUG
    G4cout << "G4HepRepFileXMLWriter:addAttValue No file is currently open" << G4endl;
#endif
  }
}

void G4HepRepFileXMLWriter::addAttValue (const char* name,
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
#ifdef G4HEPREPFILEDEBUG
    G4cout << "G4HepRepFileXMLWriter:addAttValue No file is currently open" << G4endl;
#endif
  }
}

void G4HepRepFileXMLWriter::open(const char* fileSpec)
{
  if (isOpen)
    close();
  
  fout.open(fileSpec);
    
  if (fout.good()) {
    fout << "<?xml version=\"1.0\" ?>" << G4endl;
    fout << "<heprep:heprep xmlns:heprep=\"http://www.slac.stanford.edu/~perl/heprep/\"" << G4endl;
    fout << "  xmlns:xsi=\"http://www.w3.org/1999/XMLSchema-instance\" xsi:schemaLocation=\"HepRep.xsd\">" << G4endl;
    
    isOpen = true;
    init();
  } else {
    G4cout << "G4HepRepFileXMLWriter:open Unable to write to file " << fileSpec << G4endl;
  }
}

void G4HepRepFileXMLWriter::close()
{
  // Close any remaining open Types
  endTypes();

  if (fout.good()) {      
    fout << "</heprep:heprep>" << G4endl;
    fout.close( );
    isOpen = false;
  } else {
    G4cout << "G4HepRepFileXMLWriter:close No file is currently open" << G4endl;
  }
}

void G4HepRepFileXMLWriter::endTypes()
{
  // Close any remaining open Types
    while(typeDepth>-1)
      endType();
}

void G4HepRepFileXMLWriter::endType()
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

void G4HepRepFileXMLWriter::endInstance()
{
  if (inInstance[typeDepth])
  {
    endPrimitive();
    indent();
    fout << "</heprep:instance>" << G4endl;
    inInstance[typeDepth] = false;
  }
}

void G4HepRepFileXMLWriter::endPrimitive()
{
  if (inPrimitive)
  {
    endPoint();
    indent();
    fout << "</heprep:primitive>" << G4endl;
    inPrimitive = false;
  }
}

void G4HepRepFileXMLWriter::endPoint()
{
  if (inPoint)
  {
    indent();
    fout << "</heprep:point>" << G4endl;
    inPoint = false;
  }
}

void G4HepRepFileXMLWriter::indent()
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
