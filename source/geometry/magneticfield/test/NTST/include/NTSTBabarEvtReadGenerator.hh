//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: NTSTBabarEvtReadGenerator.hh,v 1.1 2003-11-07 21:30:28 japost Exp $
//
// Description:
//	Class NTSTBabarEvtReadGenerator
//
// Environment:
//	Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//	Bill Lockman
//
// Copyright Information:
//	Copyright (C) 2000         SCIPP, U.C. Santa Cruz
//
//------------------------------------------------------------------------

#ifndef NTSTBabarEvtReadGenerator_hh
#define NTSTBabarEvtReadGenerator_hh 1

// #include <fstream>
#include <fstream.h>
//#include <g4rw/tpordvec.h>
#include "globals.hh"
#include "G4VPrimaryGenerator.hh"

class G4Event;

class NTSTBabarEvtReadGenerator:public G4VPrimaryGenerator
{
public:
  NTSTBabarEvtReadGenerator(const char* evfile);
  ~NTSTBabarEvtReadGenerator();
  
  void GeneratePrimaryVertex(G4Event* evt);
  
private:
  G4String fileName;
  ifstream inputFile;
};

#endif



