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
// $Id: G4GPRAssocT.hh,v 1.1 2007-07-27 22:13:08 tinslay Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// 
// J. Tinslay, July 2007
//
#ifndef G4GPRASSOCT_HH
#define G4GPRASSOCT_HH

#include <map>
#include <vector>
#include <algorithm>	

template <typename A, typename B>
class G4GPRAssocT {

public:

  void Register(const A& a, const B& b)
  {
    fAssoc[a] = b;
  }

  G4bool Retrieve(const A& a, B& result)
  {
    typename std::map<A, B>::iterator iter;

    if ((iter = fAssoc.find(a)) == fAssoc.end()) return false;
   
    result = iter->second;

    return true;
  }

private:

  std::map<A, B> fAssoc;
};

#endif
