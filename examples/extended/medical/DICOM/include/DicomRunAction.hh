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
//
// $Id: DicomRunAction.hh,v 1.1 2008-11-27 21:55:27 arce Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// 

#ifndef DicomRunAction_h
#define DicomRunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"
#include <vector>

class G4Run;

//=======================================================================
// DicomRunAction
//   
//
//
//=======================================================================
//
class DicomRunAction : public G4UserRunAction
{
public:
  // constructor and destructor
  DicomRunAction();
  virtual ~DicomRunAction();

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
