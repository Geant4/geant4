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
// Hadronic Inelastic Process Class
// J.L. Chuma, TRIUMF, 24-Mar-1997
// Last modified: 27-Mar-1997
// J.P. Wellisch: Bug hunting, 23-Apr-97
// Modified by J.L.Chuma 8-Jul-97 to eliminate possible division by zero for sigma
//
// 14-APR-98 F.W.Jones: variant G4HadronInelastic process for
// G4CrossSectionDataSet/DataStore class design.
//
// 17-JUN-98 F.W.Jones: removed extraneous code causing core dump.
// 01-SEP-2008 V.Ivanchenko: use methods from the base class
// 14-Sep-12 M.Kelsey -- Pass subType code to base ctor
//
 
#include "G4HadronInelasticProcess.hh"
#include "G4HadronInelasticDataSet.hh"
#include "G4GenericIon.hh"
#include "G4ParticleDefinition.hh"
  
G4HadronInelasticProcess::G4HadronInelasticProcess(const G4String& processName,
                                                   G4ParticleDefinition* aParticle):
  G4HadronicProcess(processName,fHadronInelastic)
{
  AddDataSet(new G4HadronInelasticDataSet());
  theParticle = aParticle;
}

G4HadronInelasticProcess::~G4HadronInelasticProcess() 
{}

G4bool G4HadronInelasticProcess::IsApplicable(const G4ParticleDefinition& aP)
{
  return  theParticle == &aP || theParticle == G4GenericIon::GenericIon();
}
