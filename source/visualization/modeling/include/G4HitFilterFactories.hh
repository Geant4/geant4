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
/// $Id: G4HitFilterFactories.hh 68043 2013-03-13 14:27:49Z gcosmo $
//
//
// Hit filter model factories creating filters
// and associated messengers.
//
// Jane Tinslay March 2006
//
#ifndef G4HITFILTERFACTORIES_HH
#define G4HITFILTERFACTORIES_HH

#include "G4VFilter.hh"
#include "G4VModelFactory.hh"
#include "G4VHit.hh"

// Attribute filter
class G4HitAttributeFilterFactory : public G4VModelFactory< G4VFilter<G4VHit>  > {

public: // With description

  typedef std::vector<G4UImessenger*> Messengers;
  typedef std::pair< G4VFilter<G4VHit> *, Messengers > ModelAndMessengers;

  G4HitAttributeFilterFactory();

  virtual ~G4HitAttributeFilterFactory();
  
  ModelAndMessengers Create(const G4String& placement, const G4String& name);
    
};

#endif

