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
#ifndef G4LHEPProtonBuilder_h
#define G4LHEPProtonBuilder_h 1

#include "globals.hh"

#include "G4HadronElasticProcess.hh"
#include "G4ProtonInelasticProcess.hh"
#include "G4VProtonBuilder.hh"

#include "G4LEProtonInelastic.hh"
#include "G4HEProtonInelastic.hh"

class G4LHEPProtonBuilder : public G4VProtonBuilder
{
  public: 
    G4LHEPProtonBuilder();
    virtual ~G4LHEPProtonBuilder();

  public: 
    virtual void Build(G4ProtonInelasticProcess * aP);
    virtual void Build(G4HadronElasticProcess *);
    
    void SetMinEnergy(G4double aM) 
    {
      theMin = aM;
    }

  private:
    G4LEProtonInelastic * theLEProtonModel;
    G4HEProtonInelastic * theHEProtonModel;
    
    G4double theMin;

};

// 2002 by J.P. Wellisch

#endif

