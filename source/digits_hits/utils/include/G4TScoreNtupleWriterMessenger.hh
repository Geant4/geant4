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
// Author: Ivana Hrivnacova, 30/10/2018  (ivana@ipno.in2p3.fr)

#ifndef G4TScoreNtupleWriterMessenger_h
#define G4TScoreNtupleWriterMessenger_h 1

#include "G4UImessenger.hh"

template <typename T>
class G4TScoreNtupleWriter;

class G4UIdirectory;
class G4UIcmdWithoutParameter;
class G4UIcmdWithAString;
class G4UIcmdWithAnInteger;

// class description:
//
// This is a concrete class of G4UImessenger which handles the commands for
// G4TScoreNtupleWriter. This class has the following commands:
//   /score/writerFileName filename
//   /score/writerVerbose value
//

template <typename T>
class G4TScoreNtupleWriterMessenger: public G4UImessenger
{
  public:
    G4TScoreNtupleWriterMessenger(G4TScoreNtupleWriter<T>* scoreNtupleWriter);
    ~G4TScoreNtupleWriterMessenger();
    void SetNewValue(G4UIcommand * command,G4String newValues);
  
  private:
    G4TScoreNtupleWriter<T>*  fScoreNtupleWriter;
    G4UIdirectory*        fDirectory;
    G4UIcmdWithAString*   fWriterFileNameCmd;   
    G4UIcmdWithAnInteger* fWriterVerboseCmd;
};

#include "G4TScoreNtupleWriterMessenger.icc"

#endif

