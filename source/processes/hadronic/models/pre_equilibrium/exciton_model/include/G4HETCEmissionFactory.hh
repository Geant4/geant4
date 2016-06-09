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
#ifndef G4HETCEmissionFactory_hh
#define G4HETCEmissionFactory_hh 

#include "G4VPreCompoundEmissionFactory.hh"


class G4HETCEmissionFactory : public G4VPreCompoundEmissionFactory
{
public:

  G4HETCEmissionFactory() {};
  virtual ~G4HETCEmissionFactory() {};

private:

  G4HETCEmissionFactory(const G4HETCEmissionFactory & ) : G4VPreCompoundEmissionFactory() {};
  const G4HETCEmissionFactory & operator=(const G4HETCEmissionFactory & val);
  G4bool operator==(const G4HETCEmissionFactory & val) const;
  G4bool operator!=(const G4HETCEmissionFactory & val) const;

private:

  virtual std::vector<G4VPreCompoundFragment*> *  CreateFragmentVector();

};

#endif
