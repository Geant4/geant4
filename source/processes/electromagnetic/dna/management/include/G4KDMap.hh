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
