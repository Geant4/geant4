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
// $Id: G4EvaporationGEMFactoryVI.hh 66241 2012-12-13 18:34:42Z gunter $
//
// GEM de-excitation model
// by V. Ivanchenko (July 2016)
//


#ifndef G4EvaporationGEMFactoryVI_hh
#define G4EvaporationGEMFactoryVI_hh 1

#include "G4VEvaporationFactory.hh"

class G4EvaporationGEMFactoryVI : public G4VEvaporationFactory
{
public:

  explicit G4EvaporationGEMFactoryVI(G4VEvaporationChannel* ptotoEvaporation);

  virtual ~G4EvaporationGEMFactoryVI(); 

  virtual std::vector<G4VEvaporationChannel*>* GetChannel();

private:

  G4EvaporationGEMFactoryVI(const G4EvaporationGEMFactoryVI & ) = delete;
  const G4EvaporationGEMFactoryVI & operator=
  (const G4EvaporationGEMFactoryVI & val) = delete;
  G4bool operator==(const G4EvaporationGEMFactoryVI & val) const = delete;
  G4bool operator!=(const G4EvaporationGEMFactoryVI & val) const = delete;

};

#endif
