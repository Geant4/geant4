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
#include "G4KDMap.hh"
#include "globals.hh"
#include "G4KDNode.hh"
#include <algorithm>

using namespace std;

typedef std::deque<G4KDNode_Base*>::iterator _deq_iterator;

bool __1DSortOut::sortOutNDim::operator()(G4KDNode_Base* const & lhs,
                                          G4KDNode_Base* const & rhs) //const
{
  return (*lhs)[fDimension] < (*rhs)[fDimension];
}

__1DSortOut::__1DSortOut(size_t dimension) :
    fSortOutNDim(dimension)
{
}

__1DSortOut::__1DSortOut(const __1DSortOut& right) :
    fContainer(right.fContainer), fSortOutNDim(right.fSortOutNDim)
{
}

int __1DSortOut::GetDimension()
{
  return fSortOutNDim.fDimension;
}

G4KDNode_Base* __1DSortOut::GetMidle(size_t& main_middle)
{
  size_t contSize = fContainer.size();
  main_middle = (size_t) ceil(contSize / 2.); // ceil = round up
  return fContainer[main_middle];
}

_deq_iterator __1DSortOut::Insert(G4KDNode_Base* pos)
{
  return fContainer.insert(fContainer.end(), pos);
}

G4KDNode_Base* __1DSortOut::PopOutMiddle()
{
  size_t middle;
  G4KDNode_Base* pos = GetMidle(middle);
  _deq_iterator deq_pos = fContainer.begin() + middle;

  if(deq_pos == fContainer.end()) return 0; // this is a double check

  fContainer.erase(deq_pos);
  return pos;
}

void __1DSortOut::Sort()
{
  sort(fContainer.begin(), fContainer.end(), fSortOutNDim);
}

void __1DSortOut::Erase(_deq_iterator& deq_pos)
{
  fContainer.erase(deq_pos);
}

void G4KDMap::Insert(G4KDNode_Base* pos)
{
  vector<_deq_iterator>& vit = fMap[pos];

  size_t maxSize = fSortOut.size();

  G4cout << "G4KDMap::Insert : " << maxSize << G4endl;

  vit.reserve(maxSize);

  for (size_t i = 0; i < fSortOut.size(); ++i)
  {
    vit[i] = fSortOut[i].Insert(pos);

//		if(*(vit[i]) != pos)
//		{
//			G4cout << "insert wrong iterator" << G4endl;
//			abort();
//		}
  }
  /*
   std::map<G4KDNode*, std::vector<_deq_iterator> >::iterator fMap_it
   = fMap.begin();

   for( ; fMap_it != fMap.end() ; fMap_it++)
   {
   std::vector<_deq_iterator>&  vit = fMap_it->second;

   G4KDNode* tmpNode = fMap_it->first;

   for(size_t i = 0 ; i < fSortOut.size() ; i++)
   {
   G4cout << "i = " << i << G4endl;
   G4cout << "vit[i] = " << *(vit[i]) << G4endl;
   if(*(vit[i]) != tmpNode)
   {
   G4cout << "!!!! Wrong iterator" << G4endl;
   abort();
   }
   }

   }
   */

  fIsSorted = false;
}

G4KDNode_Base* G4KDMap::PopOutMiddle(size_t dimension)
{
  G4cout << "_____________" << G4endl;
  G4cout << "G4KDMap::PopOutMiddle ( "<< dimension << " )" << G4endl;

  if(fIsSorted == false) Sort();
  G4KDNode_Base* output_node = fSortOut[dimension].PopOutMiddle();

  if(output_node == 0) return 0;

  G4cout << "output_node : " << output_node << G4endl;
  G4cout << "output_node : " << output_node->GetAxis() << G4endl;

  std::map<G4KDNode_Base*, std::vector<_deq_iterator> >::iterator fMap_it
  = fMap.find(output_node);


   if(fMap_it == fMap.end())
   {
     G4cout << "fMap_it == fMap.end()" << G4endl;
     G4cout << "output_node = " << output_node << G4endl;
     return output_node;
   }

  std::vector<_deq_iterator>& vit = fMap_it->second;

  /*
   if(fMap_it->first != output_node)
   {
   G4cout << "fMap_it->first ! output_node"<< G4endl;
   G4cout << "fMap_it->first = " << fMap_it->first << G4endl;
   abort();
   }
   */

  for(size_t i = 0; i < fSortOut.size(); i++)
  {
    if(i != dimension)
    {
      G4cout << "i = " << i << G4endl;

      /*
       // G4cout << "i = " << i << G4endl;
       // G4cout << "vit[i] = " << *(vit[i]) << G4endl;
       if(*(vit[i]) != output_node)
       {
       G4cout << "deleting wrong iterator" << G4endl;
       abort();
       }
       */
      fSortOut[i].Erase(vit[i]);
    }
  }

  fMap.erase(fMap_it);

  return output_node;
}

void G4KDMap::Sort()
{
  for (size_t i = 0; i < fSortOut.size(); ++i)
  {
    fSortOut[i].Sort();
  }

  fIsSorted = true;
}
