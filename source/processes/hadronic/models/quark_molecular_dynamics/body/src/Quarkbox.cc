#include "Quarkbox.hh"
#include "G4ios.hh"

QuarkBox::QuarkBox(Geometry& G,double T,Colour& C_) 
  : Metropolis(G,C_.List.size(),T,G.size()/20,0.05),C(C_),cnt(0)
{
  //  for (int i=0; i<C.List.size(); i++)
  //    setParticle(C.List[i]->Coordinates(),C.List[i]->Momentum());
}

double QuarkBox::Eqgp() const
{
  double E = 0;
  for (int i=0; i<C.Nquark; i++) 
    if ( C.List[i]->Color() ) {
      E += C.List[i]->E();
      for (int j=i+1; j<C.Nquark; j++) 
	E += C.E_int(i,j);
    }
  return E;
}

double QuarkBox::Ehad() const
{
  double E = 0;
  for (int i=C.Nquark; i<C.Npart; i++) 
      E += C.List[i]->E();
  return E;
}

double QuarkBox::pressure() const
{
  double x = 0;
  for (int i=0; i<C.Nquark; i++) {
    x+=(C.List[i]->Momentum()*C.List[i]->Momentum())/C.List[i]->E();
    Vektor3 F = -C.dHdx(i);
    x+=C.List[i]->Coordinates()*F;
  }
  return x/3.0/G.getVolume();
}

double QuarkBox::phad() const
{
  double x = 0;
  for (int i=0; i<C.Nquark; i++) {
    Vektor3 F = -C.dHdx(i);
    x+=C.List[i]->Coordinates()*F;
  }
  return x/3.0/G.getVolume();
}

void QuarkBox::doSomething()
{
  double e = Eqgp();
  Eq.addEntry(e);
  Eh.addEntry(Ehad());
  Pq.addEntry(pressure());
  Ph.addEntry(phad());
  double p = 0;
  double f=0;
  for (int i=0; i<C.Nquark; i++) {
    p+=length(C.List[i]->Momentum());
    f+=length(C.List[i]->Force());
  }
  Force.addEntry(f/C.Nquark);
  G4cerr << ++cnt << "  " << f/C.Nquark << "  " << e/G.getVolume() << G4endl;
}
  //  C.print(cout);

