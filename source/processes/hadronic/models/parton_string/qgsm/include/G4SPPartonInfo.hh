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
#ifndef G4SPPartonInfo_h
#define G4SPPartonInfo_h

#include "G4Types.hh"

class G4SPPartonInfo
{
  public:
    G4SPPartonInfo(G4int diq, G4int q, G4double prob)
    { diQuarkPDGCode = diq; quarkPDGCode = q; probability = prob; }
    G4int GetQuark() const {return quarkPDGCode;}
    G4int GetDiQuark() const {return diQuarkPDGCode;}
    G4double GetProbability() const {return probability;}
  private:
    G4bool operator == (const G4SPPartonInfo & aInfo) const;

  private:
    G4int quarkPDGCode;
    G4int diQuarkPDGCode;
    G4double probability;
};

#endif

