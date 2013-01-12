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
// *                                                                  *
// * Parts of this code which have been  developed by QinetiQ Ltd     *
// * under contract to the European Space Agency (ESA) are the        *
// * intellectual property of ESA. Rights to use, copy, modify and    *
// * redistribute this software for general public use are granted    *
// * in compliance with any licensing, distribution and development   *
// * policy adopted by the Geant4 Collaboration. This code has been   *
// * written by QinetiQ Ltd for the European Space Agency, under ESA  *
// * contract 17191/03/NL/LvH (Aurora Programme).                     *
// *                                                                  *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
#ifndef G4EMDissociationCrossSection_h
#define G4EMDissociationCrossSection_h 1
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// MODULE:		G4EMDissociationCrossSection.hh
//
// Version:		B.1
// Date:		15/04/04
// Author:		P R Truscott
// Organisation:	QinetiQ Ltd, UK
// Customer:		ESA/ESTEC, NOORDWIJK
// Contract:		17191/03/NL/LvH
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// CHANGE HISTORY
// --------------
//
// 17 October 2003, P R Truscott, QinetiQ Ltd, UK
// Created.
//
// 15 March 2004, P R Truscott, QinetiQ Ltd, UK
// Beta release
//
// 17 August 2011, V.Ivanchenko, provide migration to new design of cross 
//                 sections considering this cross section as element-wise
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
///////////////////////////////////////////////////////////////////////////////
//
#include "G4VCrossSectionDataSet.hh"
#include "G4DynamicParticle.hh"
#include "G4Element.hh"
#include "G4ParticleDefinition.hh"
#include "G4EMDissociationSpectrum.hh"
#include "G4PhysicsFreeVector.hh"
#include "globals.hh"

///////////////////////////////////////////////////////////////////////////////
//

class G4Material;

class G4EMDissociationCrossSection : public G4VCrossSectionDataSet
{
  public:
    G4EMDissociationCrossSection ();
    ~G4EMDissociationCrossSection ();

    virtual G4bool IsElementApplicable (const G4DynamicParticle*, G4int Z,
					const G4Material*);
 
    virtual G4double GetElementCrossSection (const G4DynamicParticle *,
					     G4int Z, const G4Material *);
      
    G4PhysicsFreeVector * GetCrossSectionForProjectile
      (G4double, G4double, G4double, G4double, G4double, G4double);
    G4PhysicsFreeVector * GetCrossSectionForTarget
      (G4double, G4double, G4double, G4double, G4double, G4double);
    G4double GetWilsonProbabilityForProtonDissociation
      (G4double, G4double);

  private:
    G4EMDissociationSpectrum *thePhotonSpectrum;
    G4double                 r0;
    G4double                 J;
    G4double                 Qprime;
    G4double                 epsilon;
    G4double                 xd;
};

#endif
