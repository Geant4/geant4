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
// INCL++ intra-nuclear cascade model
// Alain Boudard, CEA-Saclay, France
// Joseph Cugnon, University of Liege, Belgium
// Jean-Christophe David, CEA-Saclay, France
// Pekka Kaitaniemi, CEA-Saclay, France, and Helsinki Institute of Physics, Finland
// Sylvie Leray, CEA-Saclay, France
// Davide Mancusi, CEA-Saclay, France
//
#define INCLXX_IN_GEANT4_MODE 1

#include "globals.hh"

/*
 * \file G4INCLUnorderedVector.hh
 *
 * \date 2nd October 2014
 * \author Davide Mancusi
 */

#ifndef G4INCLUNORDEREDVECTOR_HH_
#define G4INCLUNORDEREDVECTOR_HH_

#include <vector>
#include <algorithm>

#if !defined(NDEBUG) && !defined(INCLXX_IN_GEANT4_MODE)
// Force instantiation of all the std::vector<Particle*> methods for debugging
// purposes
namespace G4INCL {
  class Particle;
}
template class std::vector<G4INCL::Particle*>;
#endif

namespace G4INCL {

  template<class T>
    class UnorderedVector : private std::vector<T> {
      public:
        UnorderedVector() {}
        using std::vector<T>::push_back;
        using std::vector<T>::pop_back;
        using std::vector<T>::size;
        using std::vector<T>::begin;
        using std::vector<T>::end;
        using std::vector<T>::rbegin;
        using std::vector<T>::rend;
        using std::vector<T>::front;
        using std::vector<T>::back;
        using std::vector<T>::clear;
        using std::vector<T>::empty;
        using std::vector<T>::insert;
        using std::vector<T>::erase;
        using std::vector<T>::operator[];
        using std::vector<T>::reserve;
        using std::vector<T>::resize;
        using std::vector<T>::at;
        using typename std::vector<T>::iterator;
        using typename std::vector<T>::reverse_iterator;
        using typename std::vector<T>::const_iterator;
        using typename std::vector<T>::const_reverse_iterator;

        void remove(const T &t) {
          const typename std::vector<T>::iterator removeMe = std::find(begin(), end(), t);
// assert(removeMe!=end());
          *removeMe = back();
          pop_back();
        }

        G4bool contains(const T &t) const {
          return (std::find(begin(), end(), t)!=end());
        }
    };

}

#endif // G4INCLUNORDEREDVECTOR_HH_
