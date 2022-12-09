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
// Author: Mathieu Karamitros

// The code is developed in the framework of the ESA AO7146
//
// We would be very happy hearing from you, send us your feedback! :)
//
// In order for Geant4-DNA to be maintained and still open-source,
// article citations are crucial. 
// If you use Geant4-DNA chemistry and you publish papers about your software, 
// in addition to the general paper on Geant4-DNA:
//
// Int. J. Model. Simul. Sci. Comput. 1 (2010) 157–178
//
// we would be very happy if you could please also cite the following
// reference papers on chemistry:
//
// J. Comput. Phys. 274 (2014) 841-882
// Prog. Nucl. Sci. Tec. 2 (2011) 503-508 

#ifndef G4ITMAPROOM_HH
#define G4ITMAPROOM_HH

#include <map>
//#include <unordered_map>
#include <vector>
#include <deque>
#include <cmath>
#include <iostream>

#include "G4Types.hh"

class G4KDNode_Base;

class __1DSortOut
{
public :
  __1DSortOut(std::size_t dimension);
  __1DSortOut(const __1DSortOut& right);
  G4int GetDimension();
  G4KDNode_Base* GetMidle(std::size_t& /*G4KDNode_deque*/);

  std::deque<G4KDNode_Base*>::iterator Insert(G4KDNode_Base*);
  G4KDNode_Base* PopOutMiddle();
  void Sort();
  void Erase(std::deque<G4KDNode_Base*>::iterator &);
  std::size_t Size()
  {
    return fContainer.size();
  }

protected :
  struct sortOutNDim
  {
     sortOutNDim(std::size_t dimension)
     {
  	   fDimension = dimension;
     }
     G4bool operator() (G4KDNode_Base* const& lhs, G4KDNode_Base* const& rhs);
     std::size_t fDimension;
  };

  std::deque<G4KDNode_Base*> fContainer;
  sortOutNDim fSortOutNDim;
};

class G4KDMap
{
public:
  G4KDMap(std::size_t dimensions): fSortOut(dimensions, __1DSortOut(dimensions))
  {
        fIsSorted = false;
//        for(std::size_t i = 0 ; i < dimensions ; ++i)
//        {
//            fSortOut[i] = new __1DSortOut(i);
//        }
  }

  void Insert(G4KDNode_Base* pos);
  void Sort();

  G4KDNode_Base* PopOutMiddle(std::size_t dimension);
  std::size_t GetDimension()
  {
      return fSortOut.size();
  }

  std::size_t GetSize()
  {
      return fMap.size();
  }

private:
  G4bool fIsSorted;
  std::vector<__1DSortOut> fSortOut;
  std::map<G4KDNode_Base*, std::vector<std::deque<G4KDNode_Base*>::iterator>> fMap;

  // A mettre directement dans G4KDNode
};


#endif // G4ITMAPROOM_HH
