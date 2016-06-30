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
// $Id: G4VEvaporationFactory.hh 96402 2016-04-12 07:30:48Z gcosmo $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara
//
// Modifications:
// 23 January 2012 by V.Ivanchenko added pointer of G4VPhotonEvaporation to 
//    the constructor

#ifndef G4VEvaporationFactory_hh
#define G4VEvaporationFactory_hh 1

#include "G4VEvaporationChannel.hh"
#include <vector>

class G4VEvaporationFactory
{
public:

  explicit G4VEvaporationFactory(G4VEvaporationChannel* photoEvaporation);

  virtual ~G4VEvaporationFactory();

  virtual std::vector<G4VEvaporationChannel*> * GetChannel() = 0;

private:

  G4VEvaporationFactory(const G4VEvaporationFactory & ) = delete;
  const G4VEvaporationFactory & operator=(const G4VEvaporationFactory & val) = delete;
  G4bool operator==(const G4VEvaporationFactory & val) const = delete;
  G4bool operator!=(const G4VEvaporationFactory & val) const = delete;

protected:

  G4VEvaporationChannel* thePhotonEvaporation;

};

#endif
