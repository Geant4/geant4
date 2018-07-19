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
// $Id: G4BinaryDeuteronBuilder.hh 66892 2013-01-17 10:57:59Z gunter $
//
//---------------------------------------------------------------------------
//
// ClassName:   G4BinaryDeuteronBuilder
//
// Author: 2002 H.P. Wellisch
//
// Modified:
// 30.03.2009 V.Ivanchenko create cross section by new
// 12.04.2017 A.Dotti move to new design with base class
//
//----------------------------------------------------------------------------
//
#ifndef G4BinaryDeuteronBuilder_h
#define G4BinaryDeuteronBuilder_h 

#include "globals.hh"

#include "G4HadronElasticProcess.hh"
#include "G4HadronFissionProcess.hh"
#include "G4HadronCaptureProcess.hh"
#include "G4DeuteronInelasticProcess.hh"
#include "G4VDeuteronBuilder.hh"

#include "G4BinaryCascade.hh"   

class G4BinaryDeuteronBuilder : public G4VDeuteronBuilder
{
  public: 
    G4BinaryDeuteronBuilder();
    virtual ~G4BinaryDeuteronBuilder() {}

    virtual void Build(G4HadronElasticProcess *) final override {}
    virtual void Build(G4DeuteronInelasticProcess * aP) final override;
    
    virtual void SetMinEnergy(G4double aM) final override {theMin = aM;}
    virtual void SetMaxEnergy(G4double aM) final override {theMax = aM;}

    using G4VDeuteronBuilder::Build; //Prevent compiler warning

  private:

    G4BinaryCascade * theModel;    
    G4double theMin;
    G4double theMax;
};

// 2002 by J.P. Wellisch

#endif

