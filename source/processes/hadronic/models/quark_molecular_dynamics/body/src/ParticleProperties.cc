#include "ParticleKinematics.hh"
#include "ParticleProperties.hh"
#include "G4ios.hh"
#include "math.hh"
#include "Random.hh"
#include "iso.hh"
#include "Collision.hh"
#include "EventHandling.hh"


QuantumProjections::QuantumProjections(int c_,RGB col,double i3,double s3)
  : c(c_),color(col),iso3(i3),spin3(s3) {}

QuantumProjections::QuantumProjections(const ParticleType& h)
  : c(1),color(h.getColor()),iso3(h.getIso3()),spin3(h.getSpin3()) {}

QuantumProjections& QuantumProjections::operator=(const QuantumProjections& Q)
{
  c = Q.c; color = Q.color; iso3 = Q.iso3; spin3 = Q.spin3; 
  return *this;
}

QuantumState& QuantumState::operator=(const QuantumState& x)
{
  QuantumProjections::operator=(x);
  Pattern<double>::operator=(x);
  return *this;
}

QuantumState anti(const QuantumState& p0) 
{
  QuantumState p = p0;
  p.c = -p0.c;
  p.color = -p0.color;
  p.iso3 = -p0.iso3;
  p.spin3 = -p0.spin3;
  return p;
}

QuantumState operator+(const QuantumState& p1,const QuantumState& p2)
{
  QuantumState p;
  if ( p1.C() == p2.C() ) 
    p.c = p1.C();
  else
    p.c = 1;
  p.color = p1.color+p2.color;
  p.iso3 = p1.iso3+p2.iso3;
  p.spin3 = p1.spin3+p2.spin3;

  if ( fabs(p1.iso3+p2.iso3) == p1.Isospin()+p2.Isospin() ) 
    p.setIsospin(p1.Isospin()+p2.Isospin());
  if ( fabs(p1.spin3+p2.spin3) == p1.Spin()+p2.Spin() ) 
    p.setSpin(p1.Spin()+p2.Spin());
  p.setB(p.C()*(p1.C()*p1.B()+p2.C()*p2.B()));
  p.setS(p.C()*(p1.C()*p1.S()+p2.C()*p2.S()));
  p.setC(p.C()*(p1.C()*p1.Charm()+p2.C()*p2.Charm()));
  p.setMaxColor(p.color ? 3 : 1);
  return p;
}


QuantumState addToLowest(const QuantumState& p1,const QuantumState& p2)
{
  QuantumState p;
  if ( p1.C() == p2.C() ) 
    p.c = p1.C();
  else
    p.c = 1;
  p.color = p1.color+p2.color;
  p.iso3 = p1.iso3+p2.iso3;
  p.spin3 = p1.spin3+p2.spin3;

  p.setIsospin(Iso::chooseMultiplett(p1.Isospin(),p1.iso3,p2.Isospin(),p2.iso3,Iso::LOWEST));
  p.setSpin(Iso::chooseMultiplett(p1.Spin(),p1.spin3,p2.Spin(),p2.spin3,Iso::LOWEST));
  p.setB(p.C()*(p1.C()*p1.B()+p2.C()*p2.B()));
  p.setS(p.C()*(p1.C()*p1.S()+p2.C()*p2.S()));
  p.setC(p.C()*(p1.C()*p1.Charm()+p2.C()*p2.Charm()));
  p.setMaxColor(p.color ? 3 : 1);
  return p;
}

ParticleProperties::ParticleProperties(const ParticleType& h,double Emax)
  : QuantumProjections(1,h.getColor(),h.getIso3(),h.getSpin3()),Type(h),
    mass(h.getMass(Emax)),flag(0)
{}

ParticleProperties::ParticleProperties(const ParticleType& h,const QuantumProjections& q,double Emax)
  : QuantumProjections(q),Type(h),
    mass(h.getMass(Emax)),flag(0)
{
}

ParticleProperties::ParticleProperties(const ParticleProperties& p)
  : QuantumProjections(p),Type(p.Type),
    mass(p.mass),lifetime(p.lifetime),flag(p.flag)
{
}

void ParticleProperties::writeOut(G4std::ostream& o) const
{
  Type.writeOut(o);
  o << c << "  " << color << "  " << iso3 << "  " << spin3 << "  ";
}

void ParticleProperties::ChargeProjection(double ch)
{ 
  iso3 = ch-(Type.B()+Type.S())/2.0; 
}

String ParticleProperties::Name() const 
{ 
  String s;
  if ( C()<0 && B() ) 
    s = "anti-";
  s += Type.Name()+printCharge();
  return s;
}

int ParticleProperties::PDGCode() const 
{ 
  String pdgc;
	int PDGC = 111;
  if ( C()<0 && B() ) 
    pdgc = "anti-";
  pdgc += Type.Name()+printCharge();

	if (pdgc == "Delta(1620)++") PDGC=2222 ;
	if (pdgc == "Delta(1620)+") PDGC=2122 ;
	if (pdgc == "Delta(1620)0") PDGC=1212 ;
	if (pdgc == "Delta(1620)-") PDGC=1112;
	if (pdgc == "Delta(1700)++") PDGC=12224 ;
	if (pdgc == "Delta(1700)+") PDGC=12214 ;
	if (pdgc == "Delta(1700)0") PDGC=12114 ;
	if (pdgc == "Delta(1700)-") PDGC=11114 ;
	if (pdgc == "Delta(1900)++") PDGC=2226 ;
	if (pdgc == "Delta(1900)+") PDGC=2126 ;
	if (pdgc == "Delta(1900)0") PDGC=1216 ;
	if (pdgc == "Delta(1900)-") PDGC=1116 ;
	if (pdgc == "Delta(1905)++") PDGC=2226 ;
	if (pdgc == "Delta(1905)+") PDGC=2126 ;
	if (pdgc == "Delta(1905)0") PDGC=1216 ;
	if (pdgc == "Delta(1905)-") PDGC=1116 ;
	if (pdgc == "Delta(1910)++") PDGC=22222 ;
	if (pdgc == "Delta(1910)+") PDGC=22122 ;
	if (pdgc == "Delta(1910)0") PDGC=21212 ;
	if (pdgc == "Delta(1910)-") PDGC=21112 ;
	if (pdgc == "Delta(1920)++") PDGC=22224 ;
	if (pdgc == "Delta(1920)+") PDGC=22214 ;
	if (pdgc == "Delta(1920)0") PDGC=22114 ;
	if (pdgc == "Delta(1920)-") PDGC=21114 ;
	if (pdgc == "Delta(1930)++") PDGC=12226 ;
	if (pdgc == "Delta(1930)+") PDGC=12126 ;
	if (pdgc == "Delta(1930)0") PDGC=11216 ;
	if (pdgc == "Delta(1930)-") PDGC=11116 ;
	if (pdgc == "Delta(1950)++") PDGC=2228 ;
	if (pdgc == "Delta(1950)+") PDGC=2218 ;
	if (pdgc == "Delta(1950)0") PDGC=2118 ;
	if (pdgc == "Delta(1950)-") PDGC=1118 ;
	if (pdgc == "Delta++") PDGC=2224 ;
	if (pdgc == "Delta+") PDGC=2214 ;
	if (pdgc == "Delta0") PDGC=2114 ;
	if (pdgc == "Delta-") PDGC=1114 ;

	if (pdgc == "Lambda(1405)0") PDGC=13122 ;
	if (pdgc == "Lambda(1520)0") PDGC=3124 ;
	if (pdgc == "Lambda(1600)0") PDGC=23122 ;
	if (pdgc == "Lambda(1670)0") PDGC=33122 ;
	if (pdgc == "Lambda(1690)0") PDGC=13124 ;
	if (pdgc == "Lambda(1800)0") PDGC=43122 ;
	if (pdgc == "Lambda(1810)0") PDGC=53122 ;
	if (pdgc == "Lambda(1820)0") PDGC=3126 ;
	if (pdgc == "Lambda(1830)0") PDGC=13126 ;
	if (pdgc == "Lambda(1890)0") PDGC=23124 ;
	if (pdgc == "Lambda(2100)0") PDGC=3128 ;
	if (pdgc == "Lambda(2110)0") PDGC=23126 ;

	if (pdgc == "Sigma(1385)+") PDGC=3224 ;
	if (pdgc == "Sigma(1385)0") PDGC=3214 ;
	if (pdgc == "Sigma(1385)-") PDGC=3114 ;
	if (pdgc == "Sigma(1660)+") PDGC=13222 ;
	if (pdgc == "Sigma(1660)0") PDGC=13212 ;
	if (pdgc == "Sigma(1660)-") PDGC=13112 ;
	if (pdgc == "Sigma(1670)+") PDGC=13224 ;
	if (pdgc == "Sigma(1670)0") PDGC=13214 ;
	if (pdgc == "Sigma(1670)-") PDGC=13114 ;
	if (pdgc == "Sigma(1750)+") PDGC=23222 ;
	if (pdgc == "Sigma(1750)0") PDGC=23212 ;
	if (pdgc == "Sigma(1750)-") PDGC=23112 ;
	if (pdgc == "Sigma+") PDGC=3222 ;
	if (pdgc == "Sigma0") PDGC=3212 ;
	if (pdgc == "Sigma-") PDGC=3112 ;

	if (pdgc == "Xi(1530)0") PDGC=3324 ;
	if (pdgc == "Xi(1530)-") PDGC=3314 ;

	if (pdgc == "N(1440)+") PDGC=12212 ;
	if (pdgc == "N(1440)0") PDGC=12112 ;
	if (pdgc == "N(1520)+") PDGC=2124 ;
	if (pdgc == "N(1520)0") PDGC=1214 ;
	if (pdgc == "N(1535)+") PDGC=32212 ;
	if (pdgc == "N(1535)0") PDGC=32112 ;
	if (pdgc == "N(1650)+") PDGC=32212 ;
	if (pdgc == "N(1650)0") PDGC=32112 ;
	if (pdgc == "N(1675)+") PDGC=2216 ;
	if (pdgc == "N(1675)0") PDGC=2116 ;
	if (pdgc == "N(1680)+") PDGC=12216 ;
	if (pdgc == "N(1680)0") PDGC=12116 ;
	if (pdgc == "N(1700)+") PDGC=22124 ;
	if (pdgc == "N(1700)0") PDGC=21214 ;
	if (pdgc == "N(1710)+") PDGC=42212 ;
	if (pdgc == "N(1710)0") PDGC=42112 ;
	if (pdgc == "N(1720)+") PDGC=32124 ;
	if (pdgc == "N(1720)0") PDGC=31214 ;
	if (pdgc == "N(2190)+") PDGC=2128 ;
	if (pdgc == "N(2190)0") PDGC=1218 ;

	if (pdgc == "N+") PDGC=2212 ;
	if (pdgc == "N0") PDGC=2112 ;


	if (pdgc == "anti-Delta(1620)++") PDGC=-2222 ;
	if (pdgc == "anti-Delta(1620)+") PDGC=-2122 ;
	if (pdgc == "anti-Delta(1620)0") PDGC=-1212 ;
	if (pdgc == "anti-Delta(1620)-") PDGC=-1112;
	if (pdgc == "anti-Delta(1700)++") PDGC=-12224 ;
	if (pdgc == "anti-Delta(1700)+") PDGC=-12214 ;
	if (pdgc == "anti-Delta(1700)0") PDGC=-12114 ;
	if (pdgc == "anti-Delta(1700)-") PDGC=-11114 ;
	if (pdgc == "anti-Delta(1900)++") PDGC=-2226 ;
	if (pdgc == "anti-Delta(1900)+") PDGC=-2126 ;
	if (pdgc == "anti-Delta(1900)0") PDGC=-1216 ;
	if (pdgc == "anti-Delta(1900)-") PDGC=-1116 ;
	if (pdgc == "anti-Delta(1905)++") PDGC=-2226 ;
	if (pdgc == "anti-Delta(1905)+") PDGC=-2126 ;
	if (pdgc == "anti-Delta(1905)0") PDGC=-1216 ;
	if (pdgc == "anti-Delta(1905)-") PDGC=-1116 ;
	if (pdgc == "anti-Delta(1910)++") PDGC=-22222 ;
	if (pdgc == "anti-Delta(1910)+") PDGC=-22122 ;
	if (pdgc == "anti-Delta(1910)0") PDGC=-21212 ;
	if (pdgc == "anti-Delta(1910)-") PDGC=-21112 ;
	if (pdgc == "anti-Delta(1920)++") PDGC=-22224 ;
	if (pdgc == "anti-Delta(1920)+") PDGC=-22214 ;
	if (pdgc == "anti-Delta(1920)0") PDGC=-22114 ;
	if (pdgc == "anti-Delta(1920)-") PDGC=-21114 ;
	if (pdgc == "anti-Delta(1930)++") PDGC=-12226 ;
	if (pdgc == "anti-Delta(1930)+") PDGC=-12126 ;
	if (pdgc == "anti-Delta(1930)0") PDGC=-11216 ;
	if (pdgc == "anti-Delta(1930)-") PDGC=-11116 ;
	if (pdgc == "anti-Delta(1950)++") PDGC=-2228 ;
	if (pdgc == "anti-Delta(1950)+") PDGC=-2218 ;
	if (pdgc == "anti-Delta(1950)0") PDGC=-2118 ;
	if (pdgc == "anti-Delta(1950)-") PDGC=-1118 ;
	if (pdgc == "anti-Delta++") PDGC=-2224 ;
	if (pdgc == "anti-Delta+") PDGC=-2214 ;
	if (pdgc == "anti-Delta0") PDGC=-2114 ;
	if (pdgc == "anti-Delta-") PDGC=-1114 ;

	if (pdgc == "anti-Lambda(1405)0") PDGC=-13122 ;
	if (pdgc == "anti-Lambda(1520)0") PDGC=-3124 ;
	if (pdgc == "anti-Lambda(1600)0") PDGC=-23122 ;
	if (pdgc == "anti-Lambda(1670)0") PDGC=-33122 ;
	if (pdgc == "anti-Lambda(1690)0") PDGC=-13124 ;
	if (pdgc == "anti-Lambda(1800)0") PDGC=-43122 ;
	if (pdgc == "anti-Lambda(1810)0") PDGC=-53122 ;
	if (pdgc == "anti-Lambda(1820)0") PDGC=-3126 ;
	if (pdgc == "anti-Lambda(1830)0") PDGC=-13126 ;
	if (pdgc == "anti-Lambda(1890)0") PDGC=-23124 ;
	if (pdgc == "anti-Lambda(2100)0") PDGC=-3128 ;
	if (pdgc == "anti-Lambda(2110)0") PDGC=-23126 ;

	if (pdgc == "anti-Sigma(1385)+") PDGC=-3224 ;
	if (pdgc == "anti-Sigma(1385)0") PDGC=-3214 ;
	if (pdgc == "anti-Sigma(1385)-") PDGC=-3114 ;
	if (pdgc == "anti-Sigma(1660)+") PDGC=-13222 ;
	if (pdgc == "anti-Sigma(1660)0") PDGC=-13212 ;
	if (pdgc == "anti-Sigma(1660)-") PDGC=-13112 ;
	if (pdgc == "anti-Sigma(1670)+") PDGC=-13224 ;
	if (pdgc == "anti-Sigma(1670)0") PDGC=-13214 ;
	if (pdgc == "anti-Sigma(1670)-") PDGC=-13114 ;
	if (pdgc == "anti-Sigma(1750)+") PDGC=-23222 ;
	if (pdgc == "anti-Sigma(1750)0") PDGC=-23212 ;
	if (pdgc == "anti-Sigma(1750)-") PDGC=-23112 ;
	if (pdgc == "anti-Sigma+") PDGC=-3222 ;
	if (pdgc == "anti-Sigma0") PDGC=-3212 ;
	if (pdgc == "anti-Sigma-") PDGC=-3112 ;

	if (pdgc == "anti-Xi(1530)0") PDGC=-3324 ;
	if (pdgc == "anti-Xi(1530)-") PDGC=-3314 ;

	if (pdgc == "anti-N(1440)+") PDGC=-12212 ;
	if (pdgc == "anti-N(1440)0") PDGC=-12112 ;
	if (pdgc == "anti-N(1520)+") PDGC=-2124 ;
	if (pdgc == "anti-N(1520)0") PDGC=-1214 ;
	if (pdgc == "anti-N(1535)+") PDGC=-32212 ;
	if (pdgc == "anti-N(1535)0") PDGC=-32112 ;
	if (pdgc == "anti-N(1650)+") PDGC=-32212 ;
	if (pdgc == "anti-N(1650)0") PDGC=-32112 ;
	if (pdgc == "anti-N(1675)+") PDGC=-2216 ;
	if (pdgc == "anti-N(1675)0") PDGC=-2116 ;
	if (pdgc == "anti-N(1680)+") PDGC=-12216 ;
	if (pdgc == "anti-N(1680)0") PDGC=-12116 ;
	if (pdgc == "anti-N(1700)+") PDGC=-22124 ;
	if (pdgc == "anti-N(1700)0") PDGC=-21214 ;
	if (pdgc == "anti-N(1710)+") PDGC=-42212 ;
	if (pdgc == "anti-N(1710)0") PDGC=-42112 ;
	if (pdgc == "anti-N(1720)+") PDGC=-32124 ;
	if (pdgc == "anti-N(1720)0") PDGC=-31214 ;
	if (pdgc == "anti-N(2190)+") PDGC=-2128 ;
	if (pdgc == "anti-N(2190)0") PDGC=-1218 ;

	if (pdgc == "anti-N+") PDGC=-2212 ;
	if (pdgc == "anti-N0") PDGC=-2112 ;

	if (pdgc == "Phi0") PDGC=333 ;
	if (pdgc == "K(892)+") PDGC=323 ;
	if (pdgc == "K(892)0") PDGC=313 ;
	if (pdgc == "K_bar(892)-") PDGC=-323 ;
	if (pdgc == "K_bar(892)0") PDGC=-313 ;
	if (pdgc == "K+") PDGC=321 ;
	if (pdgc == "K0") PDGC=311 ;
	if (pdgc == "K_bar-") PDGC=-321 ;
	if (pdgc == "K_bar0") PDGC=-311 ;
	if (pdgc == "a_0+") PDGC=9000211 ;
	if (pdgc == "a_0-") PDGC=-9000211 ;
	if (pdgc == "a_00") PDGC=9000111 ;
	if (pdgc == "a_1+") PDGC=20213 ;
	if (pdgc == "a_1-") PDGC=-20213 ;
	if (pdgc == "a_10") PDGC=20113 ;
	if (pdgc == "b_1+") PDGC=10213 ;
	if (pdgc == "b_1-") PDGC=-10213 ;
	if (pdgc == "b_10") PDGC= 10113;
	if (pdgc == "eta(1295)0") PDGC=100221 ;
	if (pdgc == "eta0") PDGC=221 ;
	if (pdgc == "eta_s0") PDGC=311 ;
	if (pdgc == "f_00") PDGC=9000221 ;
	if (pdgc == "f_10") PDGC=20223;
	if (pdgc == "h_10") PDGC=10223 ;
	if (pdgc == "omega(1390)0") PDGC=100223 ;
	if (pdgc == "omega0") PDGC=223 ;
	if (pdgc == "pi(1300)+") PDGC=100211 ;
	if (pdgc == "pi(1300)-") PDGC=-100211 ;
	if (pdgc == "pi(1300)0") PDGC=100111 ;
	if (pdgc == "pi+") PDGC=211 ;
	if (pdgc == "pi-") PDGC=-211 ;
	if (pdgc == "pi0") PDGC=111 ;
	if (pdgc == "rho+") PDGC=213 ;
	if (pdgc == "rho-") PDGC=-213 ;
	if (pdgc == "rho0") PDGC=113 ;

	return PDGC;
}


String ParticleProperties::printCharge() const
{
  String s;
  //  int ch = int((B() ? sign(B()) : 1)*Charge());
  int ch = (int)Charge();
  char c;
  if (ch>0) c = '+';
  if (ch<0) c = '-';
  if (ch==0) c = '0';
  int j=0;
  do {
    s+=c;
    ++j;
  }
  while (j<abs(ch));
  if ( j>2 ) {
    int ch1 = (int)Charge();
    G4cerr << "ATTENTION: " << ch << "  " << c << G4endl;
  }
  return s;
}

void ParticleProperties::setLifetime(double lt) 
{ 
  if ( lt < 0 ) 
    lifetime = Time()+Type.getLifetime(); 
  else if ( Type.getWidth() > 0 ) {
    vector<CollisionTab*>::iterator Y = CollisionTab::exists(this);
    if ( Y != CollisionTab::Root.end() ) 
      CollisionTab::remove(Y);
    lifetime = lt;
  }
  else 
    lifetime = 1e30;
  if ( lifetime<1e30 ) {
    CollisionTab::addEntry(lifetime,*this,ALL); 
  }
}
