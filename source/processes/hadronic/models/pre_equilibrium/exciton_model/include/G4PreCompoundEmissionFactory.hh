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
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
#ifndef G4PreCompoundEmissionFactory_hh
#define G4PreCompoundEmissionFactory_hh 

#include "G4VPreCompoundEmissionFactory.hh"


class G4PreCompoundEmissionFactory : public G4VPreCompoundEmissionFactory
{
public:

  G4PreCompoundEmissionFactory() {};
  virtual ~G4PreCompoundEmissionFactory() {};

private:

  G4PreCompoundEmissionFactory(const G4PreCompoundEmissionFactory & ) : G4VPreCompoundEmissionFactory() {};
  const G4PreCompoundEmissionFactory & operator=(const G4PreCompoundEmissionFactory & val);
  G4bool operator==(const G4PreCompoundEmissionFactory & val) const;
  G4bool operator!=(const G4PreCompoundEmissionFactory & val) const;

private:

  virtual std::vector<G4VPreCompoundFragment*> *  CreateFragmentVector();

};

#endif
