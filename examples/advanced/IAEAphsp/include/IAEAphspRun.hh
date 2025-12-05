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

#ifndef IAEAphspRun_h
#define IAEAphspRun_h 1

#include "G4Run.hh"

class G4Event;

class G4IAEAphspWriter;
class G4IAEAphspWriterStack;


class IAEAphspRun : public G4Run
{

public:

  // constructor and destructor.
  IAEAphspRun();
  IAEAphspRun(G4IAEAphspWriterStack* iaeaStack);
  virtual ~IAEAphspRun() override;

  // virtual method from G4Run. 
  // The method is overriden in this class for scoring.
  virtual void RecordEvent(const G4Event*) override;

  // virtual method from G4Run.
  // To merge local G4Run object into the global G4Run object.
  virtual void Merge(const G4Run*) override;

  // method to dump info into IAEAphsp output files
  void DumpToIAEAphspFiles(const G4IAEAphspWriterStack*);

  // Get/Set methods
  G4IAEAphspWriter* GetIAEAphspWriter() const   { return fIAEAphspWriter; }
  G4IAEAphspWriterStack* GetIAEAphspWriterStack() const
  { return fIAEAphspWriterStack; }


private:

  // DATA MEMBERS
  
  G4IAEAphspWriter* fIAEAphspWriter = nullptr;
  G4IAEAphspWriterStack* fIAEAphspWriterStack = nullptr;

};

#endif
