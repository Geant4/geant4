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
//---------------------------------------------------------------------------
//
// ClassName:   G4QGSPLundStrFragmProtonBuilder
//
// Author: Julia Yarba, FNAL/CD (2014)
// 12.04.2017 A.Dotti move to new design with base class
//
//----------------------------------------------------------------------------
//
#ifndef G4QGSPLundStrFragmProtonBuilder_h
#define G4QGSPLundStrFragmProtonBuilder_h 

#include "globals.hh"

#include "G4HadronElasticProcess.hh"
#include "G4HadronFissionProcess.hh"
#include "G4HadronCaptureProcess.hh"
#include "G4ProtonInelasticProcess.hh"
#include "G4VProtonBuilder.hh"

#include "G4TheoFSGenerator.hh"
#include "G4GeneratorPrecompoundInterface.hh"
#include "G4QGSModel.hh"
#include "G4QGSParticipants.hh"
#include "G4LundStringFragmentation.hh"
#include "G4ExcitedStringDecay.hh"
#include "G4QuasiElasticChannel.hh"

class G4QGSPLundStrFragmProtonBuilder : public G4VProtonBuilder
{
  public: 
  
    // ctor & dtor
    G4QGSPLundStrFragmProtonBuilder( G4bool quasiElastic=false ); 
    virtual ~G4QGSPLundStrFragmProtonBuilder();

  public: 
    virtual void Build(G4HadronElasticProcess *) final override {}
    virtual void Build(G4ProtonInelasticProcess * aP) final override;
    
    virtual void SetMinEnergy(G4double aM) final override {theMin = aM;}

    using G4VProtonBuilder::Build; //Prevent compiler warning
  private:

    G4TheoFSGenerator*               theModel;
    G4GeneratorPrecompoundInterface* theCascade;
    G4QGSModel< G4QGSParticipants >* theStringModel;
    G4ExcitedStringDecay*            theStringDecay;
    G4QuasiElasticChannel*           theQuasiElastic;
    G4LundStringFragmentation*       theStrFragm;
    G4double                         theMin;

};

#endif

