#ifndef __PERMUT__
#define __PERMUT__

#include "g4std/vector"

template<class t>
vector<t>& operator>>(vector<t>& v,int n); 

template<class t>
vector<t> operator>>(const vector<t>& v,int n);

template<class t>
vector<t> swap(const vector<t>& w,int i,int j);

template<class t>
vector< vector<t> > Permutations(const vector<t>& start);

template<class t>
bool compare(const vector<t>&,const vector<t>&);

#ifndef IS_GCC
#include "Permutations.tcc"
#endif

#endif
