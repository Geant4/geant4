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
// G4FermiBreakUpAN alternative de-excitation model
// by A. Novikov (January 2025)
//

#ifndef G4FERMIPHASEDECAY_HH
#define G4FERMIPHASEDECAY_HH

#include "G4HadPhaseSpaceKopylov.hh"

class G4FermiPhaseDecay
{
  public:
    std::vector<G4LorentzVector> CalculateDecay(const G4LorentzVector& totalMomentum,
                                                const std::vector<G4double>& fragmentsMass) const
    {
      std::vector<G4LorentzVector> results;
      KopylovDecay().Generate(totalMomentum.m(), fragmentsMass, results);
      return results;
    }

  private:
    static G4HadPhaseSpaceKopylov& KopylovDecay()
    {
      static G4HadPhaseSpaceKopylov phaseDecay;
      return phaseDecay;
    }
};

#endif  // G4FERMIPHASEDECAY_HH
