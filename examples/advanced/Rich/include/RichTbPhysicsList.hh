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
// Rich advanced example for Geant4
// RichTbPhysicsList.hh for Rich of LHCb
// History:
// Created: Sajan Easo (Sajan.Easo@cern.ch)
// Revision and changes: Patricia Mendez (Patricia.Mendez@cern.ch)
/////////////////////////////////////////////////////////////////////////////
#ifndef RichTbPhysicsList_h
#define RichTbPhysicsList_h 1

#include "globals.hh"
#include "G4VUserPhysicsList.hh"
#include "RichTbRunConfig.hh"
#include "RichTbAnalysisManager.hh"
#include "G4ParticleTable.hh"

class RichTbPhysicsList : public G4VUserPhysicsList
{
  public:
  RichTbPhysicsList();
  RichTbPhysicsList(RichTbRunConfig*);
    virtual ~RichTbPhysicsList();

  protected:
    // Construct particles and processes
    virtual void ConstructParticle();
    virtual void ConstructProcess();

    //
    virtual void SetCuts();

  protected:
    // these methods Construct particles
    virtual void ConstructBosons();
    virtual void ConstructLeptons();
    virtual void ConstructMesons();
    virtual void ConstructBaryons();

  protected:
  // these methods Construct physics processes and register them
  virtual void ConstructGeneral();
  virtual void ConstructEM();
  virtual void ConstructOp();

private:
  RichTbRunConfig* rConfigPh;
  RichTbAnalysisManager* rAnalysisPh;

    // the particle table has the complete List of existing particle types
    G4ParticleTable* theParticleTable;
    G4ParticleTable::G4PTblDicIterator* theParticleIterator;

};

#endif 




