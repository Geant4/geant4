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
// ClassName:   G4QGSBinaryKaonBuilder
//
// Author: 14-Mar-2013 A. Ribon
//
// Description: Modified version of G4QGSBinaryPiKBuilder to include on kaons.
//
// Modified:
//
//----------------------------------------------------------------------------
//
#ifndef G4QGSBinaryKaonBuilder_h
#define G4QGSBinaryKaonBuilder_h 1

#include "globals.hh"

#include "G4HadronElasticProcess.hh"
#include "G4VKaonBuilder.hh"

#include "G4TheoFSGenerator.hh"
#include "G4BinaryCascade.hh"
#include "G4QGSModel.hh"
#include "G4QGSParticipants.hh"
#include "G4QGSMFragmentation.hh"
#include "G4ExcitedStringDecay.hh"
#include "G4QuasiElasticChannel.hh"


class G4QGSBinaryKaonBuilder : public G4VKaonBuilder
{
  public: 
    G4QGSBinaryKaonBuilder(G4bool quasiElastic=false);
    virtual ~G4QGSBinaryKaonBuilder();

  public: 
    virtual void Build(G4HadronElasticProcess *) final override {};
    virtual void Build(G4HadronInelasticProcess * aP) final override;
    
    void SetMinEnergy(G4double aM) {theMin = aM;}

  private:
    G4TheoFSGenerator * theModel;
    G4double theMin;

};

#endif

