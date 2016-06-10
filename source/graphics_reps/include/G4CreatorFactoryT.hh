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
// $Id: G4CreatorFactoryT.hh 78955 2014-02-05 09:45:46Z gcosmo $
//
// Generic identifier-creator based factory. Based on 
// factory presented in "Modern C++ Design, Andrei Alexandrescu"
//
// Jane Tinslay, September 2006
//
#ifndef G4CREATORFACTORYT_HH
#define G4CREATORFACTORYT_HH

#include "globals.hh"
#include <map>

template <typename T, typename Identifier, typename Creator>
class G4CreatorFactoryT {

public:

  // Constructor
  G4CreatorFactoryT();

  // Destructor
  virtual ~G4CreatorFactoryT();

  // Register identifier<->creator pairs
  G4bool Register(const Identifier& id, Creator creator);

  // Create product with given identifier
  T* Create(const Identifier& id) const; 

private:

  typedef std::map<Identifier, Creator> Map;

  // Data member
  Map fMap;

};

template <typename T, typename Identifier, typename Creator>
G4CreatorFactoryT<T, Identifier, Creator>::G4CreatorFactoryT() {}

template <typename T, typename Identifier, typename Creator>
G4CreatorFactoryT<T, Identifier, Creator>::~G4CreatorFactoryT() {}

template <typename T, typename Identifier, typename Creator>
G4bool
G4CreatorFactoryT<T, Identifier, Creator>::Register(const Identifier& id, 
                                                    Creator creator)
{
  if (fMap.find(id) != fMap.end()) {
    G4ExceptionDescription ed;
    ed << "Creator with identifier "<<id<<" already exists."<<G4endl;
    G4Exception
      ("G4CreatorFactoryT::Register(const Identifier& id, Creator creator)",
       "greps0102", JustWarning, ed,
       "Creator exists");
    return false;
  }

  // Insert identifier<->creator pair into map
  std::pair<Identifier, Creator> myPair(id, creator);
  return fMap.insert(myPair).second;
}

template <typename T, typename Identifier, typename Creator>
T* 
G4CreatorFactoryT<T, Identifier, Creator>::Create(const Identifier& id) const
{
  typename Map::const_iterator iter = fMap.find(id);
  
  if (iter == fMap.end()) {
    G4ExceptionDescription ed;
    ed << "Identifier "<<id<<" does not exist."<<G4endl;
    G4Exception("G4CreatorFactoryT::Create(const Identifier& id)",
                "greps0103", JustWarning, ed,
                "Non-existent identifier");
    return 0;
  }
  
  return iter->second();
}

#endif
