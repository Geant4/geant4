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
// $Id: G4ITModelHandler.hh 64057 2012-10-30 15:04:49Z gcosmo $
//
// Author: Mathieu Karamitros (kara (AT) cenbg . in2p3 . fr) 
//
// WARNING : This class is released as a prototype.
// It might strongly evolve or even disapear in the next releases.
//
// History:
// -----------
// 10 Oct 2011 M.Karamitros created
//
// -------------------------------------------------------------------

#ifndef G4ITModelHandler_H
#define G4ITModelHandler_H

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
    void RegisterModel(G4VITModel* aModel, const G4double globalTime);

    // Model applying for type 1 and type 2
    inline G4ITModelManager* GetModelManager(G4ITType, G4ITType);
    void              SetModel(G4ITType, G4ITType, G4VITModel* aModel, G4double startingTime);
    G4VITModel*       GetModel(G4ITType, G4ITType, const G4double globalTime);

    //
    inline const std::vector<std::vector<G4ITModelManager*> >* GetAllModelManager()
    {
        return &fModelManager;
    }

    inline bool GetTimeStepComputerFlag() {return fTimeStepComputerFlag;}
    inline bool GetReactionProcessFlag() {return fReactionProcessFlag;}

protected:
    G4bool fIsInitialized;
    std::vector<std::vector<G4ITModelManager*> > fModelManager ;

    G4bool fTimeStepComputerFlag; // Set true if a computer is registered
    G4bool fReactionProcessFlag; // Set true if a reaction process is registered
};

inline
G4ITModelManager* G4ITModelHandler::GetModelManager(G4ITType type1, G4ITType type2)
{
    if(fModelManager.empty())
    {
        return 0;
    }

    if((int) fModelManager.size() < type1) return 0;

    std::vector<G4ITModelManager*>* v = &(fModelManager.at(type1));

    if((int) v->size() < type2) return 0;

    return v->at(type2);
}


#endif // G4ITModelHandler_H
