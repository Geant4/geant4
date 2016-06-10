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
// $Id: G4VPreCompoundEmissionFactory.hh 68028 2013-03-13 13:48:15Z gcosmo $
//
// by V. Lara
//

#ifndef G4VPreCompoundEmissionFactory_hh
#define G4VPreCompoundEmissionFactory_hh

#include "G4VPreCompoundFragment.hh"
#include <vector>

class G4VPreCompoundEmissionFactory
{
public:

  G4VPreCompoundEmissionFactory();

  virtual ~G4VPreCompoundEmissionFactory();
  
  inline std::vector<G4VPreCompoundFragment*> * GetFragmentVector();

protected:

  virtual std::vector<G4VPreCompoundFragment*> * CreateFragmentVector() = 0;

private:

  G4VPreCompoundEmissionFactory(const G4VPreCompoundEmissionFactory & );
  const G4VPreCompoundEmissionFactory & operator=
  (const G4VPreCompoundEmissionFactory & val);
  G4bool operator==(const G4VPreCompoundEmissionFactory & val) const;
  G4bool operator!=(const G4VPreCompoundEmissionFactory & val) const;

  std::vector<G4VPreCompoundFragment*> * fragvector;

  struct DeleteFragment 
  {
    template<typename T>
    void operator()(const T* ptr) const
    {
      delete ptr;
    }
  };
  
};

inline std::vector<G4VPreCompoundFragment*> * 
G4VPreCompoundEmissionFactory::GetFragmentVector()
{
  if (fragvector == 0) { fragvector = CreateFragmentVector(); }
  return fragvector;
}

#endif
