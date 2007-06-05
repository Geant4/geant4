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
//
// $Id: B01RunAction.hh,v 1.1 2007-06-05 18:20:09 ahoward Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// 

#ifndef B01RunAction_h
#define B01RunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"
#include <vector>

class G4Run;

//=======================================================================
// B01RunAction
//   
//
//
//=======================================================================
//
class B01RunAction : public G4UserRunAction
{
public:
  // constructor and destructor
  B01RunAction();
  virtual ~B01RunAction();

public:
  // virtual method from G4UserRunAction.
  virtual G4Run* GenerateRun();
  virtual void BeginOfRunAction(const G4Run*);
  virtual void EndOfRunAction(const G4Run*);

public:
  void PrintHeader(std::ostream *out);
  std::string FillString(const std::string &name, char c, G4int n, G4bool back=true);

private:
  // Data member 
  // - vector of MultiFunctionalDetecor names.
  std::vector<G4String> theSDName;  
  G4int FieldName;
  G4int FieldValue;


};

//

#endif
