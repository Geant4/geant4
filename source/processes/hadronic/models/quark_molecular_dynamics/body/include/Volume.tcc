#include "Nbody.hh"

template<class t>
Volume<t>::Volume(Geometry& G_,vector<ParticleType*>& L,REAL T_,REAL mu_,REAL mus_,int num) 
  : List(L),T(T_),mu(mu_),mu_s(mus_),Ensemble(NEW t,L),G(G_)
{
}

template<class t>
void Volume<t>::createParticles(double frac)
{
  for ( int i=0; i<List.size(); i++) {
    int N = rand_gen::Random.randInt(Ensemble.N(*List[i],parameter(T,G.getVolume(),mu,mu_s),1));
    int Nbar = -1;
    if ( List[i]->B() ) 
      Nbar = rand_gen::Random.randInt(Ensemble.N(*List[i],parameter(T,G.getVolume(),mu,mu_s),-1));
    setSpecies(i,N,Nbar,frac);
  }
}

template<class t>
void Volume<t>::setSpecies(int r,int N,int Nbar,double frac)
{
  ParticleBase *p;
  for (int i=0; i<N; i++) {
    p = makeParticle(*List[r]);
  }
  {for (int i=0; i<Nbar; i++) {
    p = makeParticle(*List[r]);
    p->Conjugate();
  }}
}

template<class t>
void QuarkVolume<t>::setSpecies(int r,int Nq,int Nqq,double frac)
{
  int N_both = min(Nq,Nqq);
  int N_add = Nq-Nqq;
  N_add -= N_add % 3; // ganze Baryonenzahl!!
  for (int i=0; i<N_both; i++) {
    int c = rand_gen::Random(1,3);
    if ( List[r]->isDiquark() ) 
      c = -c;
    Vektor4 xx;
    G.homogeneous(xx);
    for (int k=-1; k<=1; k+=2) {
      QuantumState pp(*List[r],QuantumProjections(1,c,List[r]->getIso3(),List[r]->getSpin3()));
      ParticleBase* p;
      if (k<0) 
	p = makeParticle(*List[r],anti(pp));
      else
	p = makeParticle(*List[r],pp);
      p->SetCoordinates(Vektor3(xx[1],xx[2],xx[3]));
    }
  }
  for (int j=0; j<N_add/3; j++) {
    Vektor4 xx;
    G.homogeneous(xx);
    for (int c=1; c<=3; c++) {
      double iso = ( rand_gen::Random()<frac ) ? 0.5 : -0.5;
      ParticleBase* p = makeParticle(*List[r],QuantumProjections(1,c,iso,List[r]->getSpin3()));
      p->SetCoordinates(Vektor3(xx[1],xx[2],xx[3]));
    }
  }
}

