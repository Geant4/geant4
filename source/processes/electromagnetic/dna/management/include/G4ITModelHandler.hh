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
// $Id: G4ITModelHandler.hh 100802 2016-11-02 14:55:27Z gcosmo $
//
// Author: Mathieu Karamitros

// The code is developed in the framework of the ESA AO7146
//
// We would be very happy hearing from you, send us your feedback! :)
//
// In order for Geant4-DNA to be maintained and still open-source,
// article citations are crucial. 
// If you use Geant4-DNA chemistry and you publish papers about your software, 
// in addition to the general paper on Geant4-DNA:
//
// Int. J. Model. Simul. Sci. Comput. 1 (2010) 157–178
//
// we would be very happy if you could please also cite the following
// reference papers on chemistry:
//
// J. Comput. Phys. 274 (2014) 841-882
// Prog. Nucl. Sci. Tec. 2 (2011) 503-508 

#ifndef G4ITModelHandler_H
#define G4ITModelHandler_H

#include "G4Types.hh"
#include "G4ITType.hh"
#include "G4ITModelManager.hh"

/**
 * G4ITModelHandler holds for two IT types the corresponding model manager
 */
class G4ITModelHandler
{
public:
  G4ITModelHandler();
  G4ITModelHandler(const G4ITModelHandler& other);
  G4ITModelHandler& operator=(const G4ITModelHandler& rhs);

  /** Default destructor */
  ~G4ITModelHandler();

  void Initialize();

  // Register a model at a starting time (time1)
  // if a second model is registered at a later time (time2);
  // the second model will be considered from
  // time2 to the end of simulation
  void RegisterModel(G4VITStepModel* aModel, const G4double globalTime);

  // Model applying for type 1 and type 2
  inline G4ITModelManager* GetModelManager(G4ITType, G4ITType);
  void SetModel(G4ITType,
                G4ITType,
                G4VITStepModel* aModel,
                G4double startingTime);
  G4VITStepModel* GetModel(G4ITType, G4ITType, const G4double globalTime);

  //
  inline const std::vector<std::vector<G4ITModelManager*> >* GetAllModelManager()
  {
    return &fModelManager;
  }

  inline bool GetTimeStepComputerFlag()
  {
    return fTimeStepComputerFlag;
  }
  inline bool GetReactionProcessFlag()
  {
    return fReactionProcessFlag;
  }

protected:
  G4bool fIsInitialized;
  std::vector<std::vector<G4ITModelManager*> > fModelManager;

  G4bool fTimeStepComputerFlag; // Set true if a computer is registered
  G4bool fReactionProcessFlag; // Set true if a reaction process is registered
};

inline G4ITModelManager* G4ITModelHandler::GetModelManager(G4ITType type1,
                                                           G4ITType type2)
{
  if (fModelManager.empty())
  {
    return 0;
  }

  if ((int) fModelManager.size() < type1) return 0;

  std::vector<G4ITModelManager*>* v = &(fModelManager.at(type1));

  if ((int) v->size() < type2) return 0;

  return v->at(type2);
}

#endif // G4ITModelHandler_H
