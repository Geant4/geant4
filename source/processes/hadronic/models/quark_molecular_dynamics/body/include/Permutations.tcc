template<class t>
vector<t>& operator>>(vector<t>& v,int n) 
{
  return shift(v,0,n);
}

template<class t>
vector<t>& shift(vector<t>& v,int p,int n) 
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
vector<t> operator>>(const vector<t>& v,int n) 
{
  vector<t> x = v;
  x >> n;
  return x;
}

template<class t>
vector<t> swap(const vector<t>& w,int i,int j) 
{
  vector<t> v = w;
  t x = v[i];
  v[i] = v[j];
  v[j] = x;
  return v;
}

template<class t>
void Permutate(vector< vector<t> >& P,const vector<t> start,int k)
{
  int n = start.size();
  if ( k == n-1 ) {
    P.insert(P.end(),start);
  }
  else {
    vector<t> copy(start);
    for (int i=k; i<n; i++) {
      Permutate(P,copy,k+1);
      shift(copy,k,1);
    }
  }
}

template<class t>
vector< vector<t> > Permutations(const vector<t>& start) 
{
  vector< vector<t> > P;
  Permutate(P,start,0);
  return P;
}

template<class t>
bool compare(const vector<t>& V1,const vector<t>& V2)
{
  vector< vector<t> > Perm = Permutations(V2);
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
