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
// $Id: G4EventGenerator.hh 100828 2016-11-02 15:25:59Z gcosmo $
//
#ifndef G4EventGenerator_h
#define G4EventGenerator_h 1

#include "G4HadronicInteraction.hh"
#include "G4VertexCode.hh"
#include "G4InteractionCode.hh"

class G4EventGenerator : public G4HadronicInteraction
{
  public:
    G4EventGenerator();
    ~G4EventGenerator();

  private:
    G4EventGenerator(const G4EventGenerator &right);
    const G4EventGenerator & operator=(const G4EventGenerator &right);
    int operator==(const G4EventGenerator &right) const;
    int operator!=(const G4EventGenerator &right) const;

  public:
    virtual G4double GetWidth(G4VertexCode &theCode) = 0;
    virtual G4double GetCrossSection(G4InteractionCode &theCode) = 0;
    virtual G4ReactionProduct * GetFinalState(G4InteractionCode &theCode) = 0;
    // please also define ApplyYourself for scattering off Hydrogen
};

#endif

