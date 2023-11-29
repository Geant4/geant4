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
// G4tgrSolid implementation
//
// Author: P.Arce, CIEMAT (November 2007)
// --------------------------------------------------------------------

#include <map>
#include <set>

#include "G4tgrSolid.hh"

#include "G4SystemOfUnits.hh"
#include "G4tgrUtils.hh"
#include "G4tgrMessenger.hh"
#include "G4tgrVolumeMgr.hh"

// --------------------------------------------------------------------
G4tgrSolid::G4tgrSolid()
{
}

// --------------------------------------------------------------------
G4tgrSolid::~G4tgrSolid()
{
}

// --------------------------------------------------------------------
G4tgrSolid::G4tgrSolid(const std::vector<G4String>& wl)
{
  //---------- set name
  theName = G4tgrUtils::GetString(wl[1]);

  //---------- set solid type
  theType = G4tgrUtils::GetString(wl[2]);

  //---------- create only vector<double> of theSolidParams
  FillSolidParams(wl);

  G4tgrVolumeMgr::GetInstance()->RegisterMe(this);

#ifdef G4VERBOSE
  if(G4tgrMessenger::GetVerboseLevel() >= 1)
  {
    G4cout << " Created " << *this << G4endl;
  }
#endif
}

// --------------------------------------------------------------------
const std::vector<std::vector<G4double>*> G4tgrSolid::GetSolidParams() const
{
  return theSolidParams;
}

// --------------------------------------------------------------------
const G4String& G4tgrSolid::GetRelativeRotMatName() const
{
  return theName;  // Dummy ...
}

// --------------------------------------------------------------------
G4ThreeVector G4tgrSolid::GetRelativePlace() const
{
  return G4ThreeVector(0, 0, 0);  // Dummy...
}

// --------------------------------------------------------------------
void G4tgrSolid::FillSolidParams(const std::vector<G4String>& wl)
{
  //---- Setting which are angle parameters (for dimensions...)
  std::map<G4String, std::set<G4int>> angleParams;
  std::set<G4int> apar;
  apar.clear();
  apar.insert(3);
  apar.insert(4);
  angleParams["TUBS"] = apar;
  apar.clear();
  apar.insert(5);
  apar.insert(6);
  angleParams["CONS"] = apar;
  apar.clear();
  apar.insert(3);
  apar.insert(4);
  apar.insert(5);
  angleParams["PARA"] = apar;
  apar.clear();
  apar.insert(1);
  apar.insert(2);
  apar.insert(6);
  apar.insert(10);
  angleParams["TRAP"] = apar;
  apar.clear();
  apar.insert(2);
  apar.insert(3);
  apar.insert(4);
  apar.insert(5);
  angleParams["SPHERE"] = apar;
  apar.clear();
  apar.insert(3);
  apar.insert(4);
  angleParams["TORUS"] = apar;
  apar.clear();
  apar.insert(0);
  apar.insert(1);
  angleParams["POLYCONE"] = apar;
  apar.clear();
  apar.insert(0);
  apar.insert(1);
  angleParams["POLYHEDRA"] = apar;
  apar.clear();
  apar.insert(2);
  apar.insert(3);
  angleParams["HYPE"] = apar;
  apar.clear();
  apar.insert(0);
  angleParams["TWISTED_BOX"] = apar;
  apar.clear();
  apar.insert(0);
  apar.insert(2);
  apar.insert(3);
  apar.insert(10);
  angleParams["TWISTED_TRAP"] = apar;
  apar.clear();
  apar.insert(5);
  angleParams["TWISTED_TRD"] = apar;
  apar.clear();
  apar.insert(0);
  apar.insert(4);
  angleParams["TWISTED_TUBS"] = apar;

  std::vector<G4double>* vd = new std::vector<G4double>;
  theSolidParams.push_back(vd);
  std::size_t noParRead = wl.size() - 3;

  G4String solidType = wl[2];
  //--- Default unit (mm) if length, deg if angle
  for(G4int ii = 0; ii < (G4int)noParRead; ++ii)
  {
    G4bool isAngle = 0;
    std::map<G4String, std::set<G4int>>::const_iterator ite =
      angleParams.find(solidType);
    if(ite != angleParams.cend())
    {
      std::set<G4int> apar2 = (*ite).second;
      if(apar2.find(ii) != apar2.cend())
      {
        isAngle = 1;
        vd->push_back(G4tgrUtils::GetDouble(wl[3 + ii], deg));
#ifdef G4VERBOSE
        if(G4tgrMessenger::GetVerboseLevel() >= 3)
        {
          G4cout << " G4tgrSolid::FillSolidParams() - Angle param found "
                 << solidType << " " << ii << G4endl;
        }
#endif
      }
    }
    if(!isAngle)
    {
      vd->push_back(G4tgrUtils::GetDouble(wl[3 + ii]));
    }
  }
}

// --------------------------------------------------------------------
std::ostream& operator<<(std::ostream& os, const G4tgrSolid& sol)
{
  os << "G4tgrSolid= " << sol.theName << " of type " << sol.theType
     << " PARAMS: ";
  if(sol.theSolidParams.size() != 0)
  {
    std::vector<G4double> solpar = *(sol.theSolidParams[0]);
    for(std::size_t ii = 0; ii < solpar.size(); ++ii)
    {
      os << solpar[ii] << " ";
    }
  }
  os << G4endl;

  return os;
}
