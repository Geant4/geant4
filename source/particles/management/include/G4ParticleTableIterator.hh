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
// $Id$
//
// 
// ------------------------------------------------------------
//	GEANT 4 class header file 
//
//	History: first implementation, based on object model of
//	28 Oct 1999, H.Kurashige
// ------------------------------------------------------------

#ifndef G4ParticleTableIterator_h
#define G4ParticleTableIterator_h 1

#include <map>

template < class K, class V > class G4ParticleTableIterator
{
 public:
  typedef  std::map<K, V, std::less<K> > Map;
  G4ParticleTableIterator( Map &adict):
    it(adict.begin()),
    mydict(&adict),
    defined(false)
     {}

  G4bool operator++ ()
    {
      if(!defined) return false;
      it++;
      return it!=mydict->end() ? true : false;
    }
      
  G4bool operator()()
    {
      if(defined) {
	    return operator++();
      } else {
	    defined=true;
	    it=mydict->begin();
	    return it!=mydict->end() ? true : false;
      }
    }

  void reset (){defined=false;}
  K* key() const { return &((*it).first); }
  V  value() const { return (*it).second; }

 private:
  typename Map::iterator it;  
  Map * mydict;
  G4bool defined;
};

#endif
