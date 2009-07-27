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
// $Id: G4EvaporationDefaultGEMFactory.hh,v 1.1 2009-07-27 10:20:13 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Hadronic Process: Nuclear De-excitations
// by J. M. Quesada (July 2009) on base of V. Lara code
// V.Ivanchenko cleanup
//
// new hybrid Default-GEM evaoration model:
//      - default evaporation for n,p,d,t and alpha particles
//      - GEM evaporation for light nuclei evaporation (2<Z<13,4<A<29) 


#ifndef G4EvaporationDefaultGEMFactory_h
#define G4EvaporationDefaultGEMFactory_h 1

#include "G4VEvaporationFactory.hh"

class G4EvaporationDefaultGEMFactory : public G4VEvaporationFactory
{
public:

  G4EvaporationDefaultGEMFactory();
  virtual ~G4EvaporationDefaultGEMFactory();

private:

  G4EvaporationDefaultGEMFactory(const G4EvaporationDefaultGEMFactory & );
  const G4EvaporationDefaultGEMFactory & operator=(const G4EvaporationDefaultGEMFactory & val);
  G4bool operator==(const G4EvaporationDefaultGEMFactory & val) const;
  G4bool operator!=(const G4EvaporationDefaultGEMFactory & val) const;

  std::vector<G4VEvaporationChannel*> * CreateChannel();

};

#endif
