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
#ifndef G4ITMAPROOM_HH
#define G4ITMAPROOM_HH

#include <map>
//#include <unordered_map>
#include <vector>
#include <deque>
#include <cmath>
#include <iostream>

#include "G4KDNode.hh"

class __1DSortOut
{
public :
    __1DSortOut(int dimension);
    int GetDimension();
    G4KDNode *GetMidle(int& G4KDNode_deque);

    //std::_Deque_iterator<G4KDNode*,G4KDNode*&,G4KDNode**> Insert(G4KDNode*);

    std::deque<G4KDNode*>::iterator Insert(G4KDNode*);
    G4KDNode* PopOutMiddle();
    void Sort();
    //void Erase(std::_Deque_iterator<G4KDNode*,G4KDNode*&,G4KDNode**>&);
    void Erase(std::deque<G4KDNode*>::iterator &);

protected :
    struct sortOutNDim
    {
        sortOutNDim( int dimension)
        {
            fDimension = dimension;
        }

        bool operator() (G4KDNode* const& lhs, G4KDNode* const& rhs) //const
        {
            return lhs->GetPosition()[fDimension] < rhs->GetPosition()[fDimension];
        }

        int fDimension;
    };

    std::deque<G4KDNode*> fContainer;
    sortOutNDim fSortOutNDim;
};

class G4KDMap
{
public:
    G4KDMap(int dimensions) : fSortOut(dimensions)
    {
        fIsSorted = false;
        for(int i = 0 ; i < dimensions ; i++)
        {
            fSortOut[i] = new __1DSortOut(i);
        }
    }

    void Insert(G4KDNode* pos);
    void Sort();

    G4KDNode* PopOutMiddle(int dimension);
    int GetDimension()
    {
        return fSortOut.size();
    }

    size_t GetSize()
    {
        return fMap.size();
    }

private:
    bool fIsSorted;
    std::vector<__1DSortOut*> fSortOut;
    // std::map<G4KDNode*, std::vector<std::_Deque_iterator<G4KDNode*,G4KDNode*&,G4KDNode**> > > fMap;

    std::map<G4KDNode*, std::vector<std::deque<G4KDNode*>::iterator > > fMap;
};


#endif // G4ITMAPROOM_HH
