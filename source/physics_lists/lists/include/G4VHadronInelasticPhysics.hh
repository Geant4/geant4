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
// $Id: G4VHadronInelasticPhysics.hh,v 1.1 2009-07-02 09:32:05 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//---------------------------------------------------------------------------
//
// ClassName:  G4VHadronInelasticPhysics
//
// Author: 28 June 2009 V.Ivanchenko
//
// Modified:
//
//----------------------------------------------------------------------------
//
#ifndef G4VHadronInelasticPhysics_h
#define G4VHadronInelasticPhysics_h 1

#include "G4VPhysicsConstructor.hh"
#include "G4VHadronModelBuilder.hh"
#include "globals.hh"
#include <vector>

class G4ParticleDefinition;
class G4VCrossSectionDataSet;

class G4VHadronInelasticPhysics : public G4VPhysicsConstructor
{
public: 

  G4VHadronInelasticPhysics(const G4String& name ="hInelastic", G4int verbose = 0);

  virtual ~G4VHadronInelasticPhysics();

  virtual void ConstructParticle();

  //  virtual void ConstructProcess();

  G4HadronicInteraction* BuildModel(G4VHadronModelBuilder*,
				    G4double emin, 
				    G4double emax);

  void AddInelasticCrossSection(const G4String&, 
				G4VCrossSectionDataSet*);

  void AddInelasticCrossSection(const G4ParticleDefinition*, 
				G4VCrossSectionDataSet*);

  void AddCaptureCrossSection(G4VCrossSectionDataSet*);

  void AddFissionCrossSection(G4VCrossSectionDataSet*);

private:

  // copy constructor and hide assignment operator
  G4VHadronInelasticPhysics(G4VHadronInelasticPhysics &);
  G4VHadronInelasticPhysics & operator=(const G4VHadronInelasticPhysics &right);

  std::vector<G4VHadronModelBuilder*> builders;
 
};

#endif

