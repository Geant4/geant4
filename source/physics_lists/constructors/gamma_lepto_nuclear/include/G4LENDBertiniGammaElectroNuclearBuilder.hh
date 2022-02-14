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
#ifndef G4LENDBertiniGammaElectroNuclearBuilder_h
#define G4LENDBertiniGammaElectroNuclearBuilder_h 1

#include "G4BertiniElectroNuclearBuilder.hh"

//A. Dotti (June2013): No need to change this class for MT
// Since each thread owns its own instance (created by G4EmExtraPhysics)

class G4LENDBertiniGammaElectroNuclearBuilder
:public G4BertiniElectroNuclearBuilder
{
  using base = G4BertiniElectroNuclearBuilder;
  public: 
    G4LENDBertiniGammaElectroNuclearBuilder(G4bool eNucl);
    virtual ~G4LENDBertiniGammaElectroNuclearBuilder();

  public: 
    virtual void Build();

};
#endif

