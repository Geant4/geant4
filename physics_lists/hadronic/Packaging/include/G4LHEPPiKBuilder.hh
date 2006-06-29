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
#ifndef G4LHEPPiKBuilder_h
#define G4LHEPPiKBuilder_h 1

#include "globals.hh"
#include "G4ios.hh"

#include "G4VPiKBuilder.hh"

#include "G4LEPionPlusInelastic.hh"
#include "G4LEPionMinusInelastic.hh"
#include "G4LEKaonPlusInelastic.hh"
#include "G4LEKaonZeroSInelastic.hh"
#include "G4LEKaonZeroLInelastic.hh"
#include "G4LEKaonMinusInelastic.hh"

#include "G4HEPionPlusInelastic.hh"
#include "G4HEPionMinusInelastic.hh"
#include "G4HEKaonPlusInelastic.hh"
#include "G4HEKaonZeroInelastic.hh"
#include "G4HEKaonZeroInelastic.hh"
#include "G4HEKaonMinusInelastic.hh"

class G4LHEPPiKBuilder : public G4VPiKBuilder
{
  public: 
    G4LHEPPiKBuilder();
    virtual ~G4LHEPPiKBuilder();
    
  public: 
    virtual void Build(G4HadronElasticProcess *);
    virtual void Build(G4PionPlusInelasticProcess * aP);
    virtual void Build(G4PionMinusInelasticProcess * aP);
    virtual void Build(G4KaonPlusInelasticProcess * aP);
    virtual void Build(G4KaonMinusInelasticProcess * aP);
    virtual void Build(G4KaonZeroLInelasticProcess * aP);
    virtual void Build(G4KaonZeroSInelasticProcess * aP);
    
    void SetMinEnergy(G4double aM) 
    {
      theM=aM;
      theMinPion=theM;
    }
    void SetMinPionEnergy(G4double aM)
    {
      theMinPion = aM;
    }
    
  private:
    G4double theM;
    G4double theMinPion;
    
    // Pi + 
    G4LEPionPlusInelastic* theLEPionPlusModel;
    G4HEPionPlusInelastic* theHEPionPlusModel;

    // Pi -
    G4LEPionMinusInelastic* theLEPionMinusModel;
    G4HEPionMinusInelastic* theHEPionMinusModel;

    // K + 
    G4LEKaonPlusInelastic* theLEKaonPlusModel;
    G4HEKaonPlusInelastic* theHEKaonPlusModel;

    // K -
    G4LEKaonMinusInelastic* theLEKaonMinusModel;
    G4HEKaonMinusInelastic* theHEKaonMinusModel;

    // K0L
    G4LEKaonZeroLInelastic* theLEKaonZeroLModel;
    G4HEKaonZeroInelastic* theHEKaonZeroLModel;

    // K0S
    G4LEKaonZeroSInelastic* theLEKaonZeroSModel;
    G4HEKaonZeroInelastic* theHEKaonZeroSModel;

};
// 2002 by J.P. Wellisch

#endif
