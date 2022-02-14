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
<<<<<<< HEAD:source/processes/hadronic/models/parton_string/management/include/G4VertexCode.hh
//
// $Id: G4VertexCode.hh 67999 2013-03-13 11:14:32Z gcosmo $
//
#ifndef G4VertexCode_h
#define G4VertexCode_h 1


class G4VertexCode
{

  public:
    G4VertexCode(G4String & aCode);

    void SetCode(G4String & aCode);
    G4String GetCode();
  private:

    G4String theCode;
};
=======
// 
#pragma once
#include "G4UserRunAction.hh"
#include "G4String.hh"

class G4Run;
>>>>>>> 5baee230e93612916bcea11ebf822756cfa7282c:examples/extended/medical/dna/dnadamage1/include/RunAction.hh

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class RunAction : public G4UserRunAction
{
public:
    RunAction();
    ~RunAction() override;

<<<<<<< HEAD:source/processes/hadronic/models/parton_string/management/include/G4VertexCode.hh
#endif
=======
    void BeginOfRunAction(const G4Run*) override;
    void EndOfRunAction(const G4Run*) override;

private:
    void CreateNtuple(const G4String& fileName);
    void WriteNtuple();
};
>>>>>>> 5baee230e93612916bcea11ebf822756cfa7282c:examples/extended/medical/dna/dnadamage1/include/RunAction.hh
