#ifndef __PERMUT__
#define __PERMUT__

#include "g4std/vector"

template<class t>
G4std::vector<t>& operator>>(G4std::vector<t>& v,int n); 

template<class t>
G4std::vector<t> operator>>(const G4std::vector<t>& v,int n);

template<class t>
G4std::vector<t> swap(const G4std::vector<t>& w,int i,int j);

template<class t>
G4std::vector< G4std::vector<t> > Permutations(const G4std::vector<t>& start);

template<class t>
bool compare(const G4std::vector<t>&,const G4std::vector<t>&);

// -----------------------------------------------------
// implementation from Permutations.tcc

template<class t>
G4std::vector<t>& operator>>(G4std::vector<t>& v,int n) 
{
  return shift(v,0,n);
}

template<class t>
G4std::vector<t>& shift(G4std::vector<t>& v,int p,int n) 
{
  int m = v.size()-p;
  if ( m > 1 ) {
    n = n % m;
    t x = v[v.size()-1];
    for (int i=v.size()-1; i>p; i--) {
      signed int k = i-n;
      if ( k<p )
        k += m;
      v[i] = v[k];
    }
    v[p+n-1] = x;
  }
  return v;
}

template<class t>
G4std::vector<t> operator>>(const G4std::vector<t>& v,int n) 
{
  G4std::vector<t> x = v;
  x >> n;
  return x;
}

template<class t>
G4std::vector<t> swap(const G4std::vector<t>& w,int i,int j) 
{
  G4std::vector<t> v = w;
  t x = v[i];
  v[i] = v[j];
  v[j] = x;
  return v;
}

template<class t>
void Permutate(G4std::vector< G4std::vector<t> >& P,const G4std::vector<t> start,int k)
{
  int n = start.size();
  if ( k == n-1 ) {
    P.insert(P.end(),start);
  }
  else {
    G4std::vector<t> copy(start);
    for (int i=k; i<n; i++) {
      Permutate(P,copy,k+1);
      shift(copy,k,1);
    }
  }
}

template<class t>
G4std::vector< G4std::vector<t> > Permutations(const G4std::vector<t>& start) 
{
  G4std::vector< G4std::vector<t> > P;
  Permutate(P,start,0);
  return P;
}

template<class t>
bool compare(const G4std::vector<t>& V1,const G4std::vector<t>& V2)
{
  G4std::vector< G4std::vector<t> > Perm = Permutations(V2);
  int best = -1,val = -1;
  for (unsigned int i=0; i<Perm.size(); i++) {
    bool y = false;
    for (unsigned int j=0; j<V1.size(); j++) {
      y = (*Perm[i][j] == *V1[j]);
      if ( !y ) 
        break;
    }
    if ( y ) return true;
  }
  return false;
}

// -----------------------------------------------------

#endif
