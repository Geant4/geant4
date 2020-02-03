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
// ClassName:   G4BertiniPiKBuilder
//
// Author: 2002 H.P. Wellisch
//
// Modified:
// 30.03.2009 V.Ivanchenko create cross section by new
// 12.04.2017 A.Dotti move to new design with base class
//
//----------------------------------------------------------------------------
//
#ifndef G4BertiniPiKBuilder_h
#define G4BertiniPiKBuilder_h 1

#include "globals.hh"

#include "G4HadronElasticProcess.hh"
#include "G4HadronFissionProcess.hh"
#include "G4HadronCaptureProcess.hh"
#include "G4NeutronInelasticProcess.hh"
#include "G4VPiKBuilder.hh"
#include "G4CascadeInterface.hh"   


class G4BertiniPiKBuilder : public G4VPiKBuilder
{
  public: 
    G4BertiniPiKBuilder();
    virtual ~G4BertiniPiKBuilder() {}

    virtual void Build(G4HadronElasticProcess *) final override {}
    virtual void Build(G4PionPlusInelasticProcess * aP) final override;
    virtual void Build(G4PionMinusInelasticProcess * aP) final override;
    virtual void Build(G4KaonPlusInelasticProcess * aP) final override;
    virtual void Build(G4KaonMinusInelasticProcess * aP) final override;
    virtual void Build(G4KaonZeroLInelasticProcess * aP) final override;
    virtual void Build(G4KaonZeroSInelasticProcess * aP) final override;
    
    virtual void SetMinEnergy(G4double aM) final override {theMin = aM;}
    virtual void SetMaxEnergy(G4double aM) final override {theMax = aM;}

    using G4VPiKBuilder::Build; //Prevent Compilation warning

  private:
    G4VCrossSectionDataSet * kaonxs;
    G4CascadeInterface * theModel;    
    G4double theMin;
    G4double theMax;
};

#endif

