#include "MultiIndex.hh"
#include "g4std/iostream"

G4std::ostream& operator<<(G4std::ostream& o,MultiIndex& i) {
  for (int k=0; k<i.N; k++) {
    o << i.array[k];
    if ( k<i.N-1 ) 
      o << ",";
  }
  return o;
}

MultiIndex::MultiIndex(int n,int b) 
  : N(n),array(new int[n]),n_min(new int[n]),n_max(new int[n]),valid(true)
{
  for (int i=0; i<n; i++) {
    n_min[i] = 0;
    n_max[i] = b;
    array[i] = n_min[i];
  }
}

MultiIndex::MultiIndex(int n,int* b) 
  : N(n),array(new int[n]),n_min(new int[n]),n_max(new int[n]),valid(true)
{
  for (int i=0; i<n; i++) {
    n_min[i] = 0;
    n_max[i] = b[i];
    array[i] = n_min[i];
  }
}

MultiIndex::MultiIndex(int n,int* a,int* b) 
  : N(n),array(new int[n]),n_min(new int[n]),n_max(new int[n]),valid(true)
{
  for (int i=0; i<n; i++) {
    n_min[i] = a[i];
    n_max[i] = b[i];
    array[i] = n_min[i];
  }
}

MultiIndex::MultiIndex(const MultiIndex& x) 
  : N(x.N),array(new int[N]),n_min(x.n_min),n_max(x.n_max),valid(x.valid)
{
  for (int i=0; i<N; i++) {
    n_min[i] = x.n_min[i];
    n_max[i] = x.n_max[i];
    array[i] = x.array[i];
  }
}

MultiIndex::~MultiIndex()
{
  delete [] n_min;
  delete [] n_max;
  delete [] array;
}

void MultiIndex::reset() 
{
  for (int i=0;i<N; i++)
    array[i] = n_min[i];
  valid = true;
}

void MultiIndex::increase()
{
  int k=0;
  bool go_on = true;
  do { 
    ++array[k];
    if ( array[k]>=n_max[k] ) {
      array[k] = n_min[k];
    }
    else
      go_on = false;
  }
  while ( ++k<N && go_on );
  valid = !go_on;
}

