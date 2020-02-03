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
// J.M. Quesada (August2008). Based on:
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Oct 1998)
//
// Modified:
// 17-11-2010 V.Ivanchenko integer Z and A
//

#ifndef G4ProtonEvaporationProbability_h
#define G4ProtonEvaporationProbability_h 1

#include "G4EvaporationProbability.hh"

class G4ProtonEvaporationProbability : public G4EvaporationProbability
{
public:

  explicit G4ProtonEvaporationProbability();

  ~G4ProtonEvaporationProbability() override;

protected:

  G4double CalcAlphaParam(const G4Fragment & fragment) override;
 
  G4double CalcBetaParam(const G4Fragment & fragment) override;
 
private:  

  G4ProtonEvaporationProbability(const G4ProtonEvaporationProbability &right);
  const G4ProtonEvaporationProbability & operator=
  (const G4ProtonEvaporationProbability &right);
  G4bool operator==(const G4ProtonEvaporationProbability &right) const;
  G4bool operator!=(const G4ProtonEvaporationProbability &right) const;

};

#endif
