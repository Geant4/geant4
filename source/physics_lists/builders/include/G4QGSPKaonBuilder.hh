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
// $Id$
//
//---------------------------------------------------------------------------
//
// ClassName:   G4QGSPKaonBuilder
//
// Author: 23-May-2014 D. Mancusi
//
// Description: Modified version of G4QGSBinaryKaonBuilder to use PreCompound
//
// Modified:
// 12.04.2017 A.Dotti move to new design with base class
//
//----------------------------------------------------------------------------
//
#ifndef G4QGSPKaonBuilder_h
#define G4QGSPKaonBuilder_h 1

#include "globals.hh"

#include "G4HadronElasticProcess.hh"
#include "G4VKaonBuilder.hh"

#include "G4TheoFSGenerator.hh"
#include "G4GeneratorPrecompoundInterface.hh"
#include "G4QGSModel.hh"
#include "G4QGSParticipants.hh"
#include "G4QGSMFragmentation.hh"
#include "G4ExcitedStringDecay.hh"
#include "G4QuasiElasticChannel.hh"


class G4QGSPKaonBuilder : public G4VKaonBuilder
{
  public: 
    G4QGSPKaonBuilder(G4bool quasiElastic=false);
    virtual ~G4QGSPKaonBuilder();

  public: 
    virtual void Build(G4HadronElasticProcess * aP) final override;
    virtual void Build(G4KaonPlusInelasticProcess * aP) final override;
    virtual void Build(G4KaonMinusInelasticProcess * aP) final override;
    virtual void Build(G4KaonZeroLInelasticProcess * aP) final override;
    virtual void Build(G4KaonZeroSInelasticProcess * aP) final override;
    
    virtual void SetMinEnergy(G4double aM) final override {theMin = aM;}

    using G4VKaonBuilder::Build; //Prevent compiler warning
  private:
    G4TheoFSGenerator * theModel;
    G4GeneratorPrecompoundInterface * theCascade;
    G4QGSModel< G4QGSParticipants > * theStringModel;
    G4ExcitedStringDecay * theStringDecay;
    G4QuasiElasticChannel * theQuasiElastic;
    G4double theMin;
    G4double theMax;

};

#endif

