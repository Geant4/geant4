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
// G4DNAQuadrupleIonisation.hh
//
//  Created at 2024/04/03 (Thu.)
//  Author: Shogo OKADA @KEK-CRC (shogo.okada@kek.jp)
//

#ifndef G4DNA_QUADRUPLE_IONISATION_HH_
#define G4DNA_QUADRUPLE_IONISATION_HH_

#include "G4VEmProcess.hh"
#include "G4DNAGenericIonsManager.hh"

class G4DNAQuadrupleIonisation : public G4VEmProcess {
public:
  G4DNAQuadrupleIonisation(const G4String& pname ="DNAQuadrupleIonisation",
                           G4ProcessType type = fElectromagnetic);

  virtual ~G4DNAQuadrupleIonisation() = default;

  virtual G4bool IsApplicable(const G4ParticleDefinition&);

  virtual void PrintInfo();

protected:

  virtual void InitialiseProcess(const G4ParticleDefinition*);

private:
  G4bool is_initialized_;
};

#endif // G4DNA_QUADRUPLE_IONISATION_HH_
