#include <strstream.h>
#include "globals.hh"
#include "Random.hh"
#include "reactionChannels.hh"
#include "ParticleProperties.hh"
#include "ParticleKinematics.hh"
#include "array.hh"
#include "iso.hh"
#include "Propagation.hh"
#include "Memory.hh"
#include <ctype.h>
#include "Tree.hh"
#include "EventHandling.hh"
#include "Permutations.hh"

isotop::isotop(const ParticleType& pp_,double c,int anti_) 
  : pp(pp_),charge(c),mass(-1),anti(anti_),keep(false)
{
}

isotop::isotop(const ParticleType& pp_) 
  : pp(pp_),charge(-1000),mass(-1),anti(0),keep(true)
{
}

isotop::isotop(int C,const ParticleBase& p_ ) 
  : pp(p_.getType()),charge(C*p_.Charge()),mass(-1),anti(C*p_.C()),keep(true)
{
}

bool operator==(const isotop& a,const isotop& b)
{
  return (a.pp.isEqualTo(b.pp)>0  && (a.charge == b.charge || a.charge==-1000 || b.charge==-1000) && (a.anti == b.anti || a.anti*b.anti == 0) );
}

double FunctionType::getValue(double x) 
{
  if ( x != x0 ) {
    x0 = x;
    value = f(x);
  }
  return value;
}

decayMode::decayMode(int C_,const vector<ParticleBase*>& P)
  : SumMass(0),isoSet(false),elastic(false) 
{
  for (int i=0; i<P.size(); i++) {
    products.insert(products.end(),new isotop(C_,*P[i]));
  }
}

decayMode::decayMode(istream& in,const vector<ParticleType*>& incoming) 
  : SumMass(0),isoSet(false),elastic(false)
{
  String prods,Mode,Type;
  int k = -1;
  in >> prods >> Type >> Mode;
  if ( Type[1] != ':' ) {
    istrstream in1((char*)Type);
    double f;
    in1 >> f;
    sigmaPart = new ConstantFunction(f);
  }
  if ( length(Mode) == 0 || Mode == "NORMAL" )
    select = ALL;
  else if ( Mode == "DECOMPOSITION" )
    select = DECOMPOSITION;
  else if ( Mode == "FLUXTUBE" ) 
    select = FLUXTUBE;
  else
    throw "decayMode: reaction mode undefined...";
  Array<String> parts = prods.divide(',');
  for (int i=0; i<parts.size(); i++) {
    if ( parts[i] == "in" ) {
      SumMass += incoming[++k]->PeakMass();
      products.insert(products.end(),new isotop(*incoming[k]));
      isoSet = false;
    }
    else {
      int bracketOpen = 0;
      int j;
      for (j=0; j<parts[i].getLength() && (isalpha(parts[i][j]) || parts[i][j]=='(' || parts[i][j]==')' || parts[i][j]=='_'|| bracketOpen); j++) {
	switch ( parts[i][j] ) {
	case '(' : ++bracketOpen; break;
	case ')' : --bracketOpen; break;
	}
      }
      String s = parts[i].subString(0,j-1);
      int c = 0;
      int anti = 1;
      while ( j<parts[i].getLength() ) {
	switch ( parts[i][j] ) {
	case '+' : ++c; isoSet = true; break;
	case '-' : --c; isoSet = true; break;
	case '0' : isoSet = true; break;
	case '\'' : anti = -1; break;
	default : throw "Wrong decay code!!";
	}
	++j;
      }
      ParticleType& Y = Knot<ParticleType>::FindKnot(s);
      if ( !isoSet ) c = -1000;
      isotop* p = NEW isotop(Y,c,anti);
      SumMass += p->pp.PeakMass();
      products.insert(products.end(),p);
    }
  }
}
/*
bool decayMode::compareProducts(const vector<ParticleBase*>& P) const
{
  vector<isotop*> L;
  for (int i=0; i<P.size(); i++)
    L.insert(L.end(),new isotop(*P[i]));
  return compare(L,products);
}
*/
ostream& operator<<(ostream& o,decayMode& d)
{
  for (int i=0; i<d.products.size(); i++) {
    o << d.products[i]->pp.Name();
    if ( i < d.products.size()-1 )
      o << ", ";
  }
  return o;
}

class NBodyDecay
{
  vector<Vektor4> MomList;
  const vector<double>&  MassList;
  double Decay(int n,int k,double& mmin,double mmax);
public:
  NBodyDecay(const vector<double>& M) : MassList(M) {}
  vector<Vektor4> getMomenta(double Emax);
};

double NBodyDecay::Decay(int n,int k,double& mmin,double mmax)
{
  double m_k,m_kplus1;

  if ( k>1 ) {
    if ( k == n-1 ) 
      m_kplus1 = MassList[k];
    else 
      m_kplus1 = Decay(n,k+1,mmin,mmax-MassList[k-1]);
    mmin = MassList[k-1]+m_kplus1;
    m_k = rand_gen::Random(mmin,mmax);
  }
  else {
    if ( k == n-1 ) 
      m_kplus1 = MassList[k];
    else 
      m_kplus1 = Decay(n,k+1,mmin,mmax-MassList[k-1]);
    m_k = mmax;
  }
  double p_k = CMmomentum(m_k,MassList[k-1],m_kplus1);
  // Kinematics: 
  double P_k_star_2 = sqr(p_k);
  Vektor3 P_k_star = Vektor3::isotropy(p_k);
  if ( k == n-1 ) 
     MomList.insert(MomList.begin(),Vektor4(-P_k_star,sqrt(sqr(double(MassList[k]))+P_k_star_2)));
  else {
    // beta would normally be -P_k_star/E, but the sign cancels with the
    // minus from the BACKtransformation!!
    Vektor3 beta_star = (1.0/sqrt(sqr(m_kplus1)+P_k_star_2))*P_k_star;
    for (int i=0; i<MomList.size(); i++)
      MomList[i] = ::Lorentz(beta_star,MomList[i]);
  }

MomList.insert(MomList.begin(),Vektor4(P_k_star,sqrt(sqr(double(MassList[k-1]))+P_k_star_2)));
  // End of Kinematics.
  return m_k;
}

vector<Vektor4> NBodyDecay::getMomenta(double Emax)
{
  MomList.erase(MomList.begin(),MomList.end());
  double mmin = 0;
  Decay(MassList.size(),1,mmin,Emax);
  return MomList;
}

void decayMode::performDecay(const vector<ParticleBase*>& p_list,double Etot,const Vektor3& beta,const Vektor3& x,bool force)
{
  Array<double> jk(N());
  Array<double> iso3(N());
  Array<double> spin3(N());
  Array<bool> isSet(N());
  vector<double> MassList;
  int no_in = -1;
  {for (int i=0; i<N(); i++) {
    if ( products[i]->keep ) {
      iso3[i] = p_list[++no_in]->Iso3();
      spin3[i] = p_list[no_in]->Spin3();
      isSet[i] = true;
    }
    else {
      isSet[i] = false;
    }
  }}
  QuantumState p = p_list[0]->getProperties();
  //  Now calculate the iso3-projections
  if ( !isoSet ) {
    for (int j=0; j<N(); j++) {
      jk[j] = products[j]->pp.Isospin();
    }
    Iso::Projections(N(),p.Isospin(),p.Iso3(),jk,iso3,isSet);
    //	cerr << iso3 << endl;
  }
  else {
    for (int j=0; j<N(); j++) 
      iso3[j] = products[j]->charge-0.5*(products[j]->pp.B()+products[j]->pp.S());
  }
  // and now the same thing with spin projections...
  for (int j=0; j<N(); j++) {
    jk[j] = products[j]->pp.Spin();
  }
  Iso::Projections(N(),p.Spin(),p.Spin3(),jk,spin3,isSet);
  //      cerr << spin3 << endl;
  Array<RGB> colors(N());
  if ( isDecomposition() ) {
    if ( N() == 2 ) {  // qq->q,q or M->q,q'
      if ( p.Color() != RGB::WHITE ) { // qq->q,q
	colors[0] = -p.Color() >> 1;
	colors[1] = colors[0] >> 1;
      }
      else {  // M -> q,q'
	colors[0] = int(products[0]->anti*RGB::pickColor());
	colors[1] = -colors[0];
      }
    }
    else if ( N() == 3 ) {  // B->q,q,q
      colors[0] = int(products[0]->anti*RGB::pickColor());
      colors[1] = colors[0] >> 1;
      colors[2] = colors[1] >> 1;
    }
    else
      throw "FATAL ERROR: unexcpected color combination...";
  } 
  Vektor3 Ptot;
  double E = Etot;
  for ( int i=0; i<N(); i++) {
    MassList.insert(MassList.end(),products[i]->pp.getMass(-1));
    E -= MassList[i];
  }
  //  Calculate momenta of outgoing particles:
  bool energyOk = ( E>=0 || !force );
  vector<Vektor4> MomList;
  if ( energyOk ) {
    NBodyDecay Decay(MassList);
    MomList = Decay.getMomenta(Etot);
  }
  else {
    cerr << "ATTENTION: Forcing decay. Energy not conserved!!!\n";
    for (int i=0; i<N(); i++)
      MomList.insert(MomList.end(),Vektor4(MassList[i],0,0,0));
  }
  {for (int i=0; i<N(); i++) {
    MomList[i] = ::Lorentz(beta,MomList[i]);
    const ParticleType& h = products[i]->pp;
    QuantumProjections qp(products[i]->anti*p.C(),p.C()*colors[i],iso3[i],spin3[i]);
    ParticleBase* Prod = makeParticle(h,qp);
    Prod->SetMass(MassList[i]);
    Prod->SetMomentum(spatial(MomList[i]));
    Prod->SetCoordinates(x);
    cerr << Prod->Name() << " (" << MomList[i][0] << "), ";
    Ptot += spatial(MomList[i]);
  }}
  cerr << endl << Ptot << "  " << p_list[0]->Momentum() << "  " << Etot << endl;
  cerr << endl;
}
