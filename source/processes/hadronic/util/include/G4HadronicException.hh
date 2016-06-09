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
#ifndef G4HadronicException_h
#define G4HadronicException_h

#include <exception>
#include <iostream>
#include "globals.hh"

class G4HadronicException : public std::exception
{
  public:
  G4HadronicException(G4String in, G4int at, G4String mess)
  {
    theMessage = mess;
    theName = in;
    theLine = at;
    
    if(getenv("DumpCoreOnHadronicException") )
    {
      Report(G4cout);
      G4Exception("G4HadronicException", "007", FatalException,
                  "Fatal problem in above location");
    }
    
  }
  virtual ~G4HadronicException() throw () {}
  
  void Report(std::ostream & aS)
  {
    aS<< "In " <<theName<<", line "<<theLine<<": "<<std::endl;
    aS<< "===> "<<theMessage<<std::endl;
  }
  
  private:
  G4String theMessage;
  G4String theName;
  G4int theLine;
};

#endif
