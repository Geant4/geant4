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
// GEANT4 tag $Name: not supported by cvs2svn $
//
//---------------------------------------------------------------------------
//
// ClassName:   G4LHEPAntiBarionBuilder
//
// Author: 2011 J. Apostolakis
//
// Modified:
// 22.02.2011 J.Apostolakis Started from G4MiscLHEPBuilder
//
//----------------------------------------------------------------------------
//
#ifndef G4LHEPAntiBarionBuilder_h
#define G4LHEPAntiBarionBuilder_h 1

#include "globals.hh"

#include "G4AntiProtonInelasticProcess.hh"
#include "G4AntiNeutronInelasticProcess.hh"

#include "G4LEAntiProtonInelastic.hh"
#include "G4LEAntiNeutronInelastic.hh"

// High-energy Models
#include "G4HEAntiProtonInelastic.hh"
#include "G4HEAntiNeutronInelastic.hh"

class G4LHEPAntiBarionBuilder 
{
  public: 
    G4LHEPAntiBarionBuilder();
    virtual ~G4LHEPAntiBarionBuilder();

  public: 
    void Build();

  private:
 
    // anti-proton
    G4AntiProtonInelasticProcess* theAntiProtonInelastic;
    G4LEAntiProtonInelastic* theLEAntiProtonModel;
    G4HEAntiProtonInelastic* theHEAntiProtonModel;

    // anti-neutron
    G4AntiNeutronInelasticProcess*  theAntiNeutronInelastic;
    G4LEAntiNeutronInelastic* theLEAntiNeutronModel;
    G4HEAntiNeutronInelastic* theHEAntiNeutronModel;

    G4bool wasActivated;
};
// 2011 J. Apostolakis

#endif
