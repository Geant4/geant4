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




