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

// by V. Lara


#ifndef G4VPreCompoundEmissionFactory_hh
#define G4VPreCompoundEmissionFactory_hh


#include "G4VPreCompoundFragment.hh"
#include <vector>

class G4VPreCompoundEmissionFactory
{
public:
  G4VPreCompoundEmissionFactory() : _fragvector(0) {};
  virtual ~G4VPreCompoundEmissionFactory();
  
private:
  G4VPreCompoundEmissionFactory(const G4VPreCompoundEmissionFactory & ) {};
  const G4VPreCompoundEmissionFactory & operator=(const G4VPreCompoundEmissionFactory & val);
  G4bool operator==(const G4VPreCompoundEmissionFactory & val) const;
  G4bool operator!=(const G4VPreCompoundEmissionFactory & val) const;

public:
  
  std::vector<G4VPreCompoundFragment*> * GetFragmentVector();

protected:

  virtual std::vector<G4VPreCompoundFragment*> * CreateFragmentVector() = 0;

private:
  std::vector<G4VPreCompoundFragment*> * _fragvector;

  struct DeleteFragment 
  {
    template<typename T>
    void operator()(const T* ptr) const
    {
      delete ptr;
    }
  };
  
};
#endif
