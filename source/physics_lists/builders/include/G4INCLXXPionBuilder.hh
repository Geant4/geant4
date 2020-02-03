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
#ifndef G4INCLXXPionBuilder_h
#define G4INCLXXPionBuilder_h 1

#include "globals.hh"

#include "G4HadronElasticProcess.hh"
#include "G4HadronFissionProcess.hh"
#include "G4HadronCaptureProcess.hh"
#include "G4NeutronInelasticProcess.hh"
#include "G4VPionBuilder.hh"

#include "G4INCLXXInterface.hh"   
#include "G4CascadeInterface.hh"

/**
 * Builder for pion processes using the INCL++ intra-nuclear
 * cascade model.
 *
 * By default the INCL++ model is used for projectile energies 0 - 20 GeV.
 *
 * The builder uses INCL++ cascade model with G4ExcitationHandler
 * de-excitation. This is implemented in interface
 * G4INCLXXInterface.
 *
 * @see G4INCLXXInterface
 * @see G4INCLXXProtonBuilder
 * @see G4INCLXXNeutronBuilder
 */


class G4INCLXXPionBuilder : public G4VPionBuilder
{
  public: 
    G4INCLXXPionBuilder();
    virtual ~G4INCLXXPionBuilder();

  public: 
    virtual void Build(G4HadronElasticProcess * aP);
    virtual void Build(G4PionPlusInelasticProcess * aP);
    virtual void Build(G4PionMinusInelasticProcess * aP);

    void SetMinEnergy(G4double aM) {theMin = aM;}
    void SetMaxEnergy(G4double aM) {theMax = aM;}

  private:
    G4INCLXXInterface * theModel;    
    G4double theMin;
    G4double theMax;
};

#endif

