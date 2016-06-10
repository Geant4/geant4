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
// $Id: G4EvaporationGEMFactory.hh 67983 2013-03-13 10:42:03Z gcosmo $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara
//
// Modifications:
// 23 January 2012 by V.Ivanchenko added pointer of G4VPhotonEvaporation to 
//    the constructor


#ifndef G4EvaporationGEMFactory_hh
#define G4EvaporationGEMFactory_hh 1

#include "G4VEvaporationFactory.hh"

class G4EvaporationGEMFactory : public G4VEvaporationFactory
{
public:

  G4EvaporationGEMFactory(G4VEvaporationChannel* ptotoEvaporation);

  virtual ~G4EvaporationGEMFactory(); 

  virtual std::vector<G4VEvaporationChannel*>* GetChannel();

private:

  G4EvaporationGEMFactory(const G4EvaporationGEMFactory & );
  const G4EvaporationGEMFactory & operator=(const G4EvaporationGEMFactory & val);
  G4bool operator==(const G4EvaporationGEMFactory & val) const;
  G4bool operator!=(const G4EvaporationGEMFactory & val) const;

};

#endif
