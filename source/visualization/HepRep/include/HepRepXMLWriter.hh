//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: HepRepXMLWriter.hh,v 1.1 2001-11-06 11:48:06 johna Exp $
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

#include <fstream.h>

class HepRepXMLWriter
{
public:
  HepRepXMLWriter();

  void addType(const char* name);
  void addInstance();
  void addPrimitive();
  void addPoint(double x, double y, double z);

  void addAttDef(const char* name,
		 const char* desc,
		 const char* type,
		 const char* extra);

  // Four methods to fill attValues
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
  
private:
  ofstream fout;

  bool inType;
  bool inInstance;
  bool inPrimitive;
  bool inPoint;

  void endType();
  void endInstance();
  void endPrimitive();
  void endPoint();

  void indent();
};
#endif
