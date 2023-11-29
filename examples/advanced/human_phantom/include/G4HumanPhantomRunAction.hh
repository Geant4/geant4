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
// Previous authors: G. Guerrieri, S. Guatelli and M. G. Pia, INFN Genova, Italy
// Authors (since 2007): S. Guatelli,University of Wollongong, Australia
// Contributions by F. Ambroglini INFN Perugia, Italy
//

#ifndef G4HumanPhantomRunAction_h
#define G4HumanPhantomRunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"
#include "G4HumanPhantomAnalysisManager.hh"
#include <map>

class G4HumanPhantomRunAction : public G4UserRunAction
{
  public:
    explicit G4HumanPhantomRunAction(G4HumanPhantomAnalysisManager* analysis);
   ~G4HumanPhantomRunAction();

  public:
    void BeginOfRunAction(const G4Run*);
    void EndOfRunAction(const G4Run*);
    void Fill(G4String bodypartName, G4double energyDeposit);
private:
    void totalRunEnergyDeposit();
    std::map<std::string,G4double> fEnergyTotal;
    G4String fBodypartName;

    G4HumanPhantomAnalysisManager* fAnalysisMan;
};
#endif





