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
#ifndef _NTSTFileRead_
#define _NTSTFileRead_ 1

#include <sstream>
#include "globals.hh"

#include <fstream>

class NTSTFileRead{

public:
  NTSTFileRead(const char* FileName, G4bool echo=false);
  ~NTSTFileRead();
  char* ReadLine();
  std::istringstream &StreamLine();
  
private:
  char _Line[255];
  int _LineLength;
  std::ifstream* _Istr;
  G4bool _echo;
  
  std::istringstream *stuff;
};

#endif
