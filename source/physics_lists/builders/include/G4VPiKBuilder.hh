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
// ClassName:   G4VPiKBuilder
//
// Author: 2002 J.P. Wellisch
//
// Modified:
// 16.11.2005 G.Folger: don't  keep processes as data members, but new these
// 30.03.2009 V.Ivanchenko move constructor and destructor to source
// 12.04.2017 A.Dotti move to new design with base class
//
//----------------------------------------------------------------------------
//
#ifndef G4VPiKBuilder_h
#define G4VPiKBuilder_h

#include "G4PhysicsBuilderInterface.hh"
#include "G4HadronElasticProcess.hh"
#include "G4PionPlusInelasticProcess.hh"
#include "G4PionMinusInelasticProcess.hh"
#include "G4KaonPlusInelasticProcess.hh"
#include "G4KaonMinusInelasticProcess.hh"
#include "G4KaonZeroLInelasticProcess.hh"
#include "G4KaonZeroSInelasticProcess.hh"

class G4VPiKBuilder : public G4PhysicsBuilderInterface
{
public:
  G4VPiKBuilder() = default;
  virtual ~G4VPiKBuilder() {}
  virtual void Build(G4HadronElasticProcess * aP) = 0;
  virtual void Build(G4PionPlusInelasticProcess * aP) = 0;
  virtual void Build(G4PionMinusInelasticProcess * aP) = 0;
  virtual void Build(G4KaonPlusInelasticProcess * aP) = 0;
  virtual void Build(G4KaonMinusInelasticProcess * aP) = 0;
  virtual void Build(G4KaonZeroLInelasticProcess * aP) = 0;
  virtual void Build(G4KaonZeroSInelasticProcess * aP) = 0;
  using G4PhysicsBuilderInterface::Build; //Prevent compiler warning;
};

#endif

