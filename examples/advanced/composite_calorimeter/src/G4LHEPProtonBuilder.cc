//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
#include "G4ProcessManager.hh"

G4LHEPProtonBuilder::
G4LHEPProtonBuilder() 
{
  theMin = 0;
}

G4LHEPProtonBuilder::
~G4LHEPProtonBuilder() {}

void G4LHEPProtonBuilder::
Build(G4HadronElasticProcess & aP)
{
  theElasticModel = new G4LElastic();
  aP.RegisterMe(theElasticModel);
}

void G4LHEPProtonBuilder::
Build(G4ProtonInelasticProcess & aP)
{
  theLEProtonModel = new G4LEProtonInelastic();
  theHEProtonModel = new G4HEProtonInelastic();
  theLEProtonModel->SetMinEnergy(theMin);
  theLEProtonModel->SetMaxEnergy(55*GeV);
  theHEProtonModel->SetMinEnergy(25*GeV);
  aP.RegisterMe(theLEProtonModel);
  aP.RegisterMe(theHEProtonModel);
}

// 2002 by J.P. Wellisch
