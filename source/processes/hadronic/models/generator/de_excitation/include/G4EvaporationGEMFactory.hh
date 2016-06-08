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
// * authors in the GEANT4 collaboration.                             *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4EvaporationGEMFactory.hh,v 1.1 2002/06/06 17:41:24 larazb Exp $
// GEANT4 tag $Name: geant4-04-01-patch-01 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara


#ifndef G4EvaporationGEMFactory_hh
#define G4EvaporationGEMFactory_hh


#include "G4VEvaporationFactory.hh"

class G4EvaporationGEMFactory : public G4VEvaporationFactory
{
public:
  G4EvaporationGEMFactory() {};
  virtual ~G4EvaporationGEMFactory() {};

private:
  G4EvaporationGEMFactory(const G4EvaporationGEMFactory & val) {};
  const G4EvaporationGEMFactory & operator=(const G4EvaporationGEMFactory & val);
  G4bool operator==(const G4EvaporationGEMFactory & val) const;
  G4bool operator!=(const G4EvaporationGEMFactory & val) const;

private:
  G4std::vector<G4VEvaporationChannel*> * CreateChannel();


};
#endif
