#include "newvector.hh"
#include "Random.hh"

Vektor3 Vektor3::isotropy(double r)
{
  double phi = 2*M_PI*rand_gen::Random();
  double cos_theta = 2*rand_gen::Random()-1;
  double sin_theta = sqrt(1-cos_theta*cos_theta);
  return Vektor3(r*sin_theta*cos(phi),r*sin_theta*sin(phi),r*cos_theta);
}

G4std::istream& operator>>(G4std::istream& in,Vektor& v)
{
  String vs;
  in >> vs;
  int n=1; 
  for (int i=0; i<length(vs); i++) {
    if ( vs[i] == ',' ) ++n;
  }
  v = Vektor(n,vs);
  return in;
}

Vektor::Vektor(int n,const String& s) : N(n),a(new double[N])
{
  int l=0;
  if ( s[0] != '(' || s[length(s)-1] != ')' ) 
    throw "Cannot form string!";
  String vs = s.subString(1,length(s)-2);
  int last = 0;
  for (int i=0; i<length(vs); i++) {
    if ( vs[i] == ',' ) {
      ++l;
      a[l-1] = atof(vs.subString(last,i-1));
      last = i+1;
    }
  }
  a[l] = atof(vs.subString(last));
}


void Vektor::reset(int dN,int where)
{
  int N_new = N+dN;
  if ( where<0 ) where = N;
  double* newArray = new double[N_new]; 
  if ( a ) {
    if ( where<N )
      for (int i=where; i<N; i++)
	newArray[i+dN] = a[i];
    if ( where>0 )
      for (int i=0; i<where; i++)
	newArray[i] = a[i];
    delete [] a;
  }
  for (int i=0; i<dN; i++)
    newArray[where+i] = 0.0;
  a = newArray;
  N = N_new;
}

void Vektor::remove(int pos,int dN)
{
  int N_new = N-dN;
  double* newArray = new double[N_new]; 
  if ( a ) {
    if (pos>1)
      memcpy(newArray,a,(pos-1)*sizeof(double));
    memcpy(newArray+pos-1,a+pos+dN-1,(N_new-pos+1)*sizeof(double));
    delete [] a;
  }
  a = newArray;
  N = N_new;
}

Matrize::Matrize(int n) : N(n),M(n),a(new Vektor[M])
{
  for (int i=0; i<M; i++)
    a[i].reset(N);
}

Matrize::Matrize(int n,int m) : N(n),M(m),a(new Vektor[M])
{
  for (int i=0; i<M; i++)
    a[i].reset(N);
}

void Matrize::set(int n,int m)
{
  if ( a ) 
    throw "Matrix size alread set!";
  N = n; M = m;
  a = new Vektor[M];
  for (int i=0; i<M; i++)
    a[i].reset(N);
}

Matrize operator*(const Matrize& A,const Matrize& B)
{
  Matrize C(Rows(A),Cols(B));
  for (int i=1; i<=Rows(A); i++)
    for (int j=1; j<=Cols(B); j++) {
      double s = 0.0;
      for (int k=1; k<=Cols(A); k++)
	s += A(i,k)*B(k,j);
      C(i,j) = s;
    }
  return C;
}

Vektor operator*(const Matrize& A,const Vektor& B)
{
  Vektor C(Rows(A));
  for (int i=1; i<=Rows(A); i++) {
    double s = 0.0;
    for (int k=1; k<=Cols(A); k++)
      s += A(i,k)*B[k];
    C[i] = s;
  }
  return C;
}

Matrize transpose(const Matrize& A)
{
  Matrize C(A.M,A.N);
  
  for (int i=1; i<=A.M; i++)
    for (int j=1; j<=A.N; j++)
      C(i,j) = A(j,i);
  return C;
}

G4std::ostream& operator<<(G4std::ostream& o,const Matrize& A)
{
  o << "[";
  for (int i=0; i<A.M-1; i++)
    o << A.a[i] << ",";
  o << A.a[A.M-1] << "]";
  return o;
}
