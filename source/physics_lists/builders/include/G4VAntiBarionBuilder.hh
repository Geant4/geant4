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
// $Id: G4VAntiBarionBuilder.hh 103801 2017-04-27 13:59:03Z gcosmo $
// GEANT4 tag $Name: $
//
//---------------------------------------------------------------------------
//
// ClassName:   G4VAntiBarionBuilder
//
// Author: 2011 J. Apostolakis
//
// Modified:
// 12.04.2017 A.Dotti move to new design with base class
//
//----------------------------------------------------------------------------
//
#ifndef G4VAntiBarionBuilder_h
#define G4VAntiBarionBuilder_h

#include "G4PhysicsBuilderInterface.hh"
#include "G4HadronElasticProcess.hh"
#include "G4AntiProtonInelasticProcess.hh"
#include "G4AntiNeutronInelasticProcess.hh"
#include "G4AntiDeuteronInelasticProcess.hh"
#include "G4AntiTritonInelasticProcess.hh"
#include "G4AntiHe3InelasticProcess.hh"
#include "G4AntiAlphaInelasticProcess.hh"

class G4VAntiBarionBuilder : public G4PhysicsBuilderInterface
{
public:
  G4VAntiBarionBuilder() = default;
  virtual ~G4VAntiBarionBuilder() {} 
  virtual void Build(G4HadronElasticProcess * aP) = 0;
  virtual void Build(G4AntiProtonInelasticProcess * aP) = 0;
  virtual void Build(G4AntiNeutronInelasticProcess * aP) = 0;
  virtual void Build(G4AntiDeuteronInelasticProcess * aP) = 0;
  virtual void Build(G4AntiTritonInelasticProcess * aP) = 0;
  virtual void Build(G4AntiHe3InelasticProcess * aP) = 0;
  virtual void Build(G4AntiAlphaInelasticProcess * aP) = 0;
  using G4PhysicsBuilderInterface::Build; //Prevent compiler warning
};

#endif
