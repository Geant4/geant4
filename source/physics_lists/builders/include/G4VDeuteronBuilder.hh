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
//---------------------------------------------------------------------------
//
// ClassName:   G4VDeuteronBuilder
//
// Author: 2013 P. Arce
//
// Modified
// 12.04.2017 A.Dotti move to new design with base class
//
//----------------------------------------------------------------------------
//
#ifndef G4VDeuteronBuilder_h
#define G4VDeuteronBuilder_h

#include "G4PhysicsBuilderInterface.hh"

class G4DeuteronInelasticProcess;
class G4HadronElasticProcess;

class G4VDeuteronBuilder : public G4PhysicsBuilderInterface
{
public:
  G4VDeuteronBuilder() = default;
  virtual ~G4VDeuteronBuilder() {} 
  virtual void Build(G4HadronElasticProcess * aP) = 0;
  virtual void Build(G4DeuteronInelasticProcess * aP) = 0;
  using G4PhysicsBuilderInterface::Build; //Prevent compilation warning
};

#endif
