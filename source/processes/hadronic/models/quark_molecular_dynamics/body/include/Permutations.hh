#ifndef __PERMUT__
#define __PERMUT__

#include <vector>

template<class t>
std::vector<t>& operator>>(std::vector<t>& v,int n); 

template<class t>
std::vector<t> operator>>(const std::vector<t>& v,int n);

template<class t>
std::vector<t> swap(const std::vector<t>& w,int i,int j);

template<class t>
std::vector< std::vector<t> > Permutations(const std::vector<t>& start);

template<class t>
bool compare(const std::vector<t>&,const std::vector<t>&);

// -----------------------------------------------------
// implementation from Permutations.tcc

template<class t>
std::vector<t>& operator>>(std::vector<t>& v,int n) 
{
  return shift(v,0,n);
}

template<class t>
std::vector<t>& shift(std::vector<t>& v,int p,int n) 
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
std::vector<t> operator>>(const std::vector<t>& v,int n) 
{
  std::vector<t> x = v;
  x >> n;
  return x;
}

template<class t>
std::vector<t> swap(const std::vector<t>& w,int i,int j) 
{
  std::vector<t> v = w;
  t x = v[i];
  v[i] = v[j];
  v[j] = x;
  return v;
}

template<class t>
void Permutate(std::vector< std::vector<t> >& P,const std::vector<t> start,int k)
{
  int n = start.size();
  if ( k == n-1 ) {
    P.insert(P.end(),start);
  }
  else {
    std::vector<t> copy(start);
    for (int i=k; i<n; i++) {
      Permutate(P,copy,k+1);
      shift(copy,k,1);
    }
  }
}

template<class t>
std::vector< std::vector<t> > Permutations(const std::vector<t>& start) 
{
  std::vector< std::vector<t> > P;
  Permutate(P,start,0);
  return P;
}

template<class t>
bool compare(const std::vector<t>& V1,const std::vector<t>& V2)
{
  std::vector< std::vector<t> > Perm = Permutations(V2);
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
