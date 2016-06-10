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
// $Id: G4PreCompoundEmissionFactory.hh 68028 2013-03-13 13:48:15Z gcosmo $
//
//

#ifndef G4PreCompoundEmissionFactory_hh
#define G4PreCompoundEmissionFactory_hh 

#include "G4VPreCompoundEmissionFactory.hh"

class G4PreCompoundEmissionFactory : public G4VPreCompoundEmissionFactory
{
public:

  G4PreCompoundEmissionFactory();
  virtual ~G4PreCompoundEmissionFactory();

protected:

  virtual std::vector<G4VPreCompoundFragment*> *  CreateFragmentVector();

private:

  G4PreCompoundEmissionFactory(const G4PreCompoundEmissionFactory & );
  const G4PreCompoundEmissionFactory & operator=
  (const G4PreCompoundEmissionFactory & val);
  G4bool operator==(const G4PreCompoundEmissionFactory & val) const;
  G4bool operator!=(const G4PreCompoundEmissionFactory & val) const;

};

#endif
