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
// $Id: G4VKaonBuilder.hh 103801 2017-04-27 13:59:03Z gcosmo $
//
//---------------------------------------------------------------------------
//
// ClassName:   G4VKaonBuilder
//
// Author: 2010 G.Folger
//  devired from G4VPiKBuilder
//
// Modified:
// 12.04.2017 A.Dotti move to new design with base class
//
//----------------------------------------------------------------------------
//
#ifndef G4VKaonBuilder_h
#define G4VKaonBuilder_h

#include "G4PhysicsBuilderInterface.hh"
#include "G4HadronElasticProcess.hh"
#include "G4KaonPlusInelasticProcess.hh"
#include "G4KaonMinusInelasticProcess.hh"
#include "G4KaonZeroLInelasticProcess.hh"
#include "G4KaonZeroSInelasticProcess.hh"

class G4VKaonBuilder : public G4PhysicsBuilderInterface
{
public:
  G4VKaonBuilder() = default;
  virtual ~G4VKaonBuilder() {} 
  virtual void Build(G4HadronElasticProcess * aP) = 0;
  virtual void Build(G4KaonPlusInelasticProcess * aP) = 0;
  virtual void Build(G4KaonMinusInelasticProcess * aP) = 0;
  virtual void Build(G4KaonZeroLInelasticProcess * aP) = 0;
  virtual void Build(G4KaonZeroSInelasticProcess * aP) = 0;
  using G4PhysicsBuilderInterface::Build; //Prevent compiler warning
};
#endif
