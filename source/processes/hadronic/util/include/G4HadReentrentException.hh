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
#ifndef G4HadReentrentException_h
#define G4HadReentrentException_h

#include <exception>
#include <iostream>
#include "G4HadronicException.hh"

class G4HadReentrentException : public G4HadronicException
{
  public:
  G4HadReentrentException(G4String in, G4int at, G4String mess)
  : G4HadronicException(in, at, mess) {}
  virtual ~G4HadReentrentException() throw () {}
  
  void Report(std::ostream & aS)
  {
    aS<< "G4HadReentrentException:" <<std::endl;
    G4HadronicException::Report(aS);
  }
};

#endif
