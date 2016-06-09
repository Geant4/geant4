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
// $Id: TiaraPrimaryGeneratorAction.hh,v 1.1.1.2 2006/06/29 15:44:07 gunter Exp $
// GEANT4 tag $Name: geant4-09-01 $
//

#ifndef TiaraPrimaryGeneratorAction_hh
#define TiaraPrimaryGeneratorAction_hh TiaraPrimaryGeneratorAction_hh 

#include "G4VUserPrimaryGeneratorAction.hh"

#include "TiaraTally.hh"


class G4ParticleGun;
class G4Event;
class TiaraVSourceEnergyGenerator;
class TiaraVDirectionGenerator;
class TiaraDimensions;

class TiaraPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
  TiaraPrimaryGeneratorAction(const TiaraVSourceEnergyGenerator& eG,
			      const TiaraVDirectionGenerator& dG,
			      const TiaraTally &tally,
			      const TiaraDimensions &tiaraDimensions);
  TiaraPrimaryGeneratorAction(const TiaraPrimaryGeneratorAction &rhs);

  ~TiaraPrimaryGeneratorAction();
  
  void GeneratePrimaries(G4Event* anEvent);
  const TiaraVSourceEnergyGenerator *GetEnergyGenerator() const;
  const TiaraTally &GetTally() const;


  TiaraPrimaryGeneratorAction &operator=
  (const TiaraPrimaryGeneratorAction &rhs);
  
private:
  TiaraVSourceEnergyGenerator *fEnergyGenerator;
  TiaraVDirectionGenerator *fDirectionGenerator;
  G4ParticleGun* particleGun;
  TiaraTally fTally;
};

#endif
