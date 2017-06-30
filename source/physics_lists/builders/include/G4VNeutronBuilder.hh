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
// $Id: G4VNeutronBuilder.hh 103801 2017-04-27 13:59:03Z gcosmo $
//
//---------------------------------------------------------------------------
//
// ClassName:   G4VNeutronBuilder
//
// Author: 2002 J.P. Wellisch
//
// Modified:
// 21.11.2005 G.Folger: don't  keep processes as data members, but new these
// 30.03.2009 V.Ivanchenko move constructor and destructor to source
// 12.04.2017 A.Dotti move to new design with base class
//
//----------------------------------------------------------------------------
//
#ifndef G4VNeutronBuilder_h
#define G4VNeutronBuilder_h

#include "G4PhysicsBuilderInterface.hh"

class G4HadronElasticProcess;
class G4HadronFissionProcess;
class G4HadronCaptureProcess;
class G4NeutronInelasticProcess;

class G4VNeutronBuilder : public G4PhysicsBuilderInterface
{
public:
  G4VNeutronBuilder() = default;
  virtual ~G4VNeutronBuilder() {} 
  virtual void Build(G4HadronElasticProcess * aP) = 0;
  virtual void Build(G4HadronFissionProcess * aP) = 0;
  virtual void Build(G4HadronCaptureProcess * aP) = 0;
  virtual void Build(G4NeutronInelasticProcess * aP) = 0;
  using G4PhysicsBuilderInterface::Build; //Prevent compiler warning

};
// 2002 by J.P. Wellisch

#endif
