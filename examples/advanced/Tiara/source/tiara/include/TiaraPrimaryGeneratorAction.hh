//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: TiaraPrimaryGeneratorAction.hh,v 1.1.1.1 2003/06/12 13:08:24 dressel Exp $
// GEANT4 tag $Name: geant4-07-00-cand-01 $
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
