//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: HepRepXMLWriter.cc,v 1.1 2001-11-06 11:48:09 johna Exp $
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
#include <fstream.h>
#include <unistd.h>

#include "HepRepXMLWriter.hh"

HepRepXMLWriter::HepRepXMLWriter()
{
}

void HepRepXMLWriter::addType(const char* name)
{
  if (fout.good())
  {
    endType();
    fout << "  <heprep:type version=\"null\" name=\"" << name << "\">" << endl;
    inType = true;
  } else {
    cout << "HepRepXMLWriter:addType No file is currently open." << endl;
  }
}

void HepRepXMLWriter::addInstance()
{
  if (fout.good())
  {
    if (inType)
    {
      endInstance();
      fout << "    <heprep:instance>" << endl;
      inInstance = true;
    } else {
      cout << "HepRepXMLWriter:addInstance No HepRep Type is currently open" << endl;
    }
  } else {
    cout << "HepRepXMLWriter:addInstance No file is currently open" << endl;
  }
}

void HepRepXMLWriter::addPrimitive()
{
  if (fout.good())
  {
    if (inInstance)
    {
      endPrimitive();
      fout << "      <heprep:primitive>" << endl;
      inPrimitive = true;
    } else {
      cout << "HepRepXMLWriter:addPrimitive No HepRep Instance is currently open" << endl;
    }
  } else {
    cout << "HepRepXMLWriter:addPrimitive No file is currently open" << endl;
  }
}

void HepRepXMLWriter::addPoint(double x, double y, double z)
{
  if (fout.good())
  {
    if (inPrimitive)
    {
      endPoint();
      fout << "        <heprep:point x=\"" << x << "\" y=\"" << y << "\" z=\"" << z << "\">" << endl;
      inPoint = true;
    } else {
      cout << "HepRepXMLWriter:addPoint No HepRep Primitive is currently open" << endl;
    }
  } else {
    cout << "HepRepXMLWriter:addPoint No file is currently open" << endl;
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
    fout << "  <heprep:attdef extra=\"" << extra << "\" name=\"" << name << "\" type=\"" << type << "\"" << endl;
    indent();
    fout << "    desc=\"" << desc << "\"/>" << endl;
  } else {
    cout << "HepRepXMLWriter:addAttDef No file is currently open" << endl;
  }
}

// Four methods to fill attValues
void HepRepXMLWriter::addAttValue (const char* name,
		  const char* value)
{
  if (fout.good())
  {
    indent();
    fout << "  <heprep:attvalue showLabel=\"NONE\" name=\"" << name << "\"" << endl;
    indent();
    fout << "    value=\"" << value << "\"/>" << endl;
  } else {
    cout << "HepRepXMLWriter:addAttValue No file is currently open" << endl;
  }
}

void HepRepXMLWriter::addAttValue (const char* name,
		  double value)
{
  if (fout.good())
  {
    indent();
    fout << "  <heprep:attvalue showLabel=\"NONE\" name=\"" << name << "\"" << endl;
    indent();
    fout << "    value=\"" << value << "\"/>" << endl;
  } else {
    cout << "HepRepXMLWriter:addAttValue No file is currently open" << endl;
  }
}

void HepRepXMLWriter::addAttValue (const char* name,
		  int value)
{
  if (fout.good())
  {
    indent();
    fout << "  <heprep:attvalue showLabel=\"NONE\" name=\"" << name << "\"" << endl;
    indent();
    fout << "    value=\"" << value << "\"/>" << endl;
  } else {
    cout << "HepRepXMLWriter:addAttValue No file is currently open" << endl;
  }
}

void HepRepXMLWriter::addAttValue (const char* name,
		  bool value)
{
  if (fout.good())
  {
    indent();
    fout << "  <heprep:attvalue showLabel=\"NONE\" name=\"" << name << "\"" << endl;
    indent();
    fout << "    value=\"" << value << "\"/>" << endl;
  } else {
    cout << "HepRepXMLWriter:addAttValue No file is currently open" << endl;
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
    fout << "  <heprep:attvalue showLabel=\"NONE\" name=\"" << name << "\"" << endl;
    indent();
    fout << "    value=\"" << redness << "," << greenness << "," << blueness << "\"/>" << endl;
  } else {
    cout << "HepRepXMLWriter:addAttValue No file is currently open" << endl;
  }
}

void HepRepXMLWriter::open(const char* filespec)
{
  if (fout.good())
    fout.close();

  fout.open(filespec);

  if (fout.good())
  {
    fout << "<?xml version=\"1.0\" ?>" << endl;
    fout << "<heprep:heprep xmlns:heprep=\"http://www.freehep.org/HepRep\"" << endl;
    fout << "  xmlns:xsi=\"http://www.w3.org/1999/XMLSchema-instance\" xsi:schemaLocation=\"HepRep.xsd\">" << endl;

    inType = false;
    inInstance = false;
    inPrimitive = false;
    inPoint = false;
  } else {
    cout << "HepRepXMLWriter:open Unable to write to file " << filespec << endl;
  }
}

void HepRepXMLWriter::close()
{
  if (fout.good())
  {
    endType();
    fout << "</heprep:heprep>" << endl;
    fout.close( );
  } else {
    cout << "HepRepXMLWriter:close No file is currently open" << endl;
  }
}

void HepRepXMLWriter::endType()
{
  if (inType)
  {
    endInstance();
    fout << "  </heprep:type>" << endl;
    inType = false;
  }
}

void HepRepXMLWriter::endInstance()
{
  if (inInstance)
  {
    endPrimitive();
    fout << "    </heprep:instance>" << endl;
    inInstance = false;
  }
}

void HepRepXMLWriter::endPrimitive()
{
  if (inPrimitive)
  {
    endPoint();
    fout << "      </heprep:primitive>" << endl;
    inPrimitive = false;
  }
}

void HepRepXMLWriter::endPoint()
{
  if (inPoint)
  {
    fout << "        </heprep:point>" << endl;
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
