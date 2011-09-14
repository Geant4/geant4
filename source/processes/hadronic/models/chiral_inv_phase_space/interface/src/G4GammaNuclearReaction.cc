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
// Short description: The CHIPS model provides the G4QHadronVector
// output, which is converted to the G4 particle-singletons
// --------------------------------------------------------------------
//
// Created: J.P. Wellisch, 2000/08/18 
// 01.09.2008 V.Ivanchenko 
//

#include "G4GammaNuclearReaction.hh"
#include "G4Gamma.hh"
#include "G4Nucleus.hh"
#include "G4HadFinalState.hh"
#include "G4HadProjectile.hh"
#include <iostream>


G4GammaNuclearReaction::G4GammaNuclearReaction(const G4String& name): 
  G4HadronicInteraction(name)
{
  Description();
}


G4GammaNuclearReaction::~G4GammaNuclearReaction()
{}


void G4GammaNuclearReaction::Description() const
{
  char* dirName = getenv("G4PhysListDocDir");
  if (dirName) {
    std::ofstream outFile;
    G4String outFileName = GetModelName() + ".html";
    G4String pathName = G4String(dirName) + "/" + outFileName;

    outFile.open(pathName);
    outFile << "<html>\n";
    outFile << "<head>\n";

    outFile << "<title>Description of CHIPS Gamma Nuclear Model</title>\n";
    outFile << "</head>\n";
    outFile << "<body>\n";

    outFile << "G4GammaNuclearReaction handles inelastic gamma scattering\n"
            << "using the Chiral Invariant Phase Space (CHIPS) model of\n"
            << "M. Kossov.  The incident gamma is interacted directly with\n"
            << "nucleons and clusters of nucleons within the target by\n"
            << "forming quasmons (or generalized excited hadrons) which then\n"
            << "decay into final state hadrons.  The model is valid for\n"
            << "incident gamma energies up to 3.5 GeV.\n";

    outFile << "</body>\n";
    outFile << "</html>\n";
    outFile.close();
  }
}


G4HadFinalState * G4GammaNuclearReaction::ApplyYourself(
	   const G4HadProjectile& aTrack, 
	   G4Nucleus& aTargetNucleus)
{
  if(aTrack.GetDefinition() != G4Gamma::GammaDefinition())
  {
    throw G4HadronicException(__FILE__, __LINE__, 
			      "Called G4GammaNuclearReaction for particle other than gamma");
  }
  return theModel.ApplyYourself(aTrack, aTargetNucleus);
}

