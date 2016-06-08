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
//
#include "G4ProtonBuilder.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4ProcessManager.hh"

G4ProtonBuilder::
G4ProtonBuilder() {}

G4ProtonBuilder::
~G4ProtonBuilder() {}

void G4ProtonBuilder::
Build()
{
  G4std::vector<G4VProtonBuilder *>::iterator i;
  for(i=theModelCollections.begin(); i!=theModelCollections.end(); i++)
  {
    (*i)->Build(theProtonElasticProcess);
    (*i)->Build(theProtonInelastic);
  }
  G4ProcessManager * theProcMan = G4Proton::Proton()->GetProcessManager();
  theProcMan->AddDiscreteProcess(&theProtonElasticProcess);
  theProcMan->AddDiscreteProcess(&theProtonInelastic);
}
// 2002 by J.P. Wellisch
