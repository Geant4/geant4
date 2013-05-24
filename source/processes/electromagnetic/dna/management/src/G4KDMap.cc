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
