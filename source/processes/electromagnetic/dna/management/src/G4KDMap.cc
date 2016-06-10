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
#include <algorithm>

using namespace std;

//typedef std::_Deque_iterator<G4KDNode*,G4KDNode*&,G4KDNode**> _deq_iterator ;
typedef std::deque<G4KDNode*>::iterator _deq_iterator ;

__1DSortOut::__1DSortOut(int dimension) : fSortOutNDim(dimension)
{}

int __1DSortOut::GetDimension()
{
    return fSortOutNDim.fDimension;
}

G4KDNode* __1DSortOut::GetMidle(int& main_middle)
{
    int contSize = fContainer.size();
    main_middle = (int) ceil(contSize/2.); // ceil = round up
    return fContainer[main_middle];
}

_deq_iterator __1DSortOut::Insert(G4KDNode* pos)
{
    return fContainer.insert(fContainer.end(),pos);
}

G4KDNode* __1DSortOut::PopOutMiddle()
{
    int middle;
    G4KDNode* pos = GetMidle(middle);
    _deq_iterator deq_pos = fContainer.begin()+middle;
    fContainer.erase(deq_pos);
    return pos;
}

void __1DSortOut::Sort()
{
    sort(fContainer.begin(),fContainer.end(),fSortOutNDim);
}

void __1DSortOut::Erase(_deq_iterator& deq_pos)
{
	fContainer.erase(deq_pos);
}

void G4KDMap::Insert(G4KDNode* pos)
{
    vector<_deq_iterator>& vit = fMap[pos];

    size_t maxSize = fSortOut.size();

    cout <<  maxSize << endl;
    cout <<  fSortOut.capacity() << endl;

    vit.reserve(maxSize);

    for (size_t i = 0; i < fSortOut.size(); ++i)
    {
        vit[i]=fSortOut[i]->Insert(pos);
    }

    fIsSorted = false;
}

G4KDNode* G4KDMap::PopOutMiddle(int dimension)
{
    if(fIsSorted == false) Sort();
    G4KDNode* output_node = fSortOut[dimension]->PopOutMiddle();

    std::map<G4KDNode*, std::vector<_deq_iterator> >::iterator fMap_it
    = fMap.find(output_node);

    std::vector<_deq_iterator>&  vit = fMap_it->second;

    for(int i = 0 ; i < (int) fSortOut.size() ; i++)
    {
    	if(i != dimension) fSortOut[i]->Erase(vit[i]);
    }

    fMap.erase(fMap_it);

    return output_node;
}

void G4KDMap::Sort()
{
    for (size_t i = 0; i < fSortOut.size(); ++i)
    {
        fSortOut[i]->Sort();
    }

    fIsSorted = true;
}
