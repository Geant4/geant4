#include "g4std/fstream"
#include <algo.h>
#include "newvector.H"
#include "Random.H"
#include "newBinning.H"
#include "Geometry.H"
#include "Arguments.H"
#include "HadronGas.H"
#include "MathTools.H"
#include "Error.H"
//#include "outputList.H"
#include "HeatBath.H"
#include "Propagation.H"
#include "ParticleBase.H"
#include "reactionChannels.H"
#include "ParticleKinematics.H"
#include "iso.H"
#include "array.H"
#include "Memory.H"
#include "Collision.H"
#include "genericRead.H"
#include "StreamBuffer.H"
#include "RunningVariable.H"
#include "InputVariable.H"
#include "InputReader.H"
#include "Volume.H"
#include "Potentials.H"
#include "Quarkbox.H"
#include "String.H"
#include "output.H"

Array<String> xx;
 
G4std::ostream* Output::fileout = 0;

class q_Quark : public QuantumNumbers
{
public:
  q_Quark() { setB(1.0/3.0); setS(0); }
};

class Pion : public QuantumNumbers
{
public:
  Pion() { setB(0.0); setS(0); setIsospin(1.0); setSpin(0.0); setPeakMass(0.138); }
};

class Proton : public QuantumNumbers
{
public:
  Proton() { setB(1.0); setS(0); setIsospin(0.5); setSpin(0.5); setPeakMass(0.938); }
};

class transverseMomentum : private InverseFunction
{
  double N,T;
  virtual REAL ToBeInverted(REAL x) const { return 1-exp(-x)*(1+x)/N; }
public:
  transverseMomentum(double m,double t) 
    : InverseFunction(m/t,20),T(t),N(exp(-m/t)*(1+m/t)) {}
  double getValue(double x) const { return T*Inverse(x); }
  double getValue() const { return T*Inverse(rand_gen::Random()); }
};

void setPhaseSpace(vector<Particle*>& L,Geometry& G,REAL T) 
{
  Vektor4 x;
  HeatBath<Relativistic> Distribution;
  for (int i=0; i<L.size(); i++) {
    L[i]->SetMomentum(Distribution(T,L[i]->Mass()));
    //    p->SetMomentum(Distribution(Relativistic(T,(*X).mass)));
    G.homogeneous(x);
    L[i]->SetCoordinates(Vektor3(x[1],x[2],x[3]));
  }
}

void v_profile(vector<Particle*>& L,double tau,double alpha = 1.0) 
{
  for (int i=0; i<L.size(); i++) {
    double z = L[i]->Coordinates(3);
    L[i]->Momentum(3) = z/tau/sqrt(1.0-(z*z/tau/tau))*L[i]->TransverseMass();
  }
}

void trans_profile(vector<Particle*>& L,double R,double alpha = 1.0) 
{
  for (int i=0; i<L.size(); i++) {
    double mz=L[i]->Mass();
    double r = sqrt(sqr(L[i]->Coordinates(1))+sqr(L[i]->Coordinates(2)));
    double br = 0.6*pow(r/R,alpha);
    double pt = br/sqrt(1-sqr(br))*mz;
    double x = 2*M_PI*rand_gen::Random();
    Vektor3 rvek = L[i]->Coordinates();
    rvek[3] = 0;
    double pz = L[i]->Momentum(3);
    L[i]->SetMomentum(pt*rvek/length(rvek));
    L[i]->Momentum(3) = pz;
  }
}

int N_c;

void setPhaseSpace_stream(vector<Particle*>& L,Geometry& G,REAL T,double ymax) 
{
  Vektor4 x;
  for (int i=0; i<L.size(); i++) {
    G.homogeneous(x);
    L[i]->SetCoordinates(Vektor3(x[1],x[2],x[3]));
    double betaz = 2.0*x[3]/((Tube&)G).L;
    double gammaz = 1.0/sqrt(1.0-betaz*betaz);
    transverseMomentum Distribution(L[i]->Mass(),T);
    Vektor3 p;
    double mt = Distribution.getValue();
    double pt = sqrt(sqr(mt)-sqr(L[i]->Mass()));
    double r = 2*mathConstants::Pi*rand_gen::Random();
    p[1] = pt*sin(r);
    p[2] = pt*cos(r);
    p[3] = gammaz*betaz*mt;
    L[i]->SetMomentum(p);
    if ( L[i]->Charm() != 0 ) 
      ++N_c;
  }
}

void waitfor(G4std::istream& in,const String &s) {
  char c;
  while ( in ) {
    while ( in.get(c) && c != s[0] );
    int i;
    for (i=1; in.get(c) && i<length(s) && c == s[i]; i++);
    if ( i == length(s) ) {
      s.putToStream(in);
      break;
    }
  }
}

void skipline(G4std::istream& in,int n=1) {
  char c;
  for (int i=0; i<n; i++) {
    //        G4cerr << "skip: ";
    while ( in.get(c) && c != '\n' ) /*G4cerr << c*/ ;
    //        G4cerr << G4endl;
  }
}

double readEvent(G4std::istream& in,bool final = true) 
{
  bool readit = false;
  int S,C,col;
  double B,m,i,i3,s,s3,lt,time;
  String name,fs;
  Vektor3 xs,ps;
  bool over = false;
  bool data = false;
  while ( in ) {
    in >> name;
    if ( name == "#" ) {
      if ( data ) {
	in.putback('#');
	break;
      }
      in >> fs;
      if ( final && fs.subString(0,4) == "final" ) 
	readit = true;
      else if ( fs.subString(0,3) == "time" ) {
	time = atof(fs.subString(5,length(fs)-2));
	G4cerr << "t = " << time << G4endl;
	skipline(in);
      }
      else 
	skipline(in);
    }
    else if ( readit ) {
      data = true;
      in >> B >> S >> C >> m >> i >> i3 >> s >> s3;
      in >> xs >> ps >> col >> fs >> lt;
      skipline(in);
      if ( fabs(i3) > i ) { i3 = i; }
      if ( fabs(s3) > s ) { s3 = s; }
      int last;
      if ( name[length(name)-2] == '0' ) 
	last = length(name)-3;
      else 
	for (last=length(name)-2; last>0 && ( name[last] == '-' || name[last] == '+'); last--);
      name = name.subString(0,last);
      double sig = sign(B);
      if ( sig == 0 ) sig = 1.0;
      if ( name.subString(0,4) == "anti-" ) 
	name = name.subString(5);

      ParticleType& h = Knot<ParticleType>::FindKnot(name);
      ParticleBase* p = makeParticle(h,QuantumProjections(sig,col,i3,s3),ps,xs,m);
      p->setLifetime(lt);
    }
    else 
      skipline(in);
  }
  if ( data ) 
    return time;
  else
    return 0.0;
}
/*
double readUQMD(G4std::istream& in) 
{
  ParticleType* quarks[3];
  const ParticleType& q = Knot<ParticleType>::FindKnot("q");
  const ParticleType& s = Knot<ParticleType>::FindKnot("s");
  quarks[0] = (ParticleType*)&q;
  quarks[1] = (ParticleType*)&q;
  quarks[2] = (ParticleType*)&s;
  int N;
  double t,p0,mass,dummy,iso;
  int i3,Q,q1,q2;
  Vektor3 x,p;
  skipline(in,14);
  in >> N;
  G4cerr << "N=" << N << G4endl;
  skipline(in,2);
  for (int j=0; j<N; j++) {
    in >> t >> x[1] >> x[2] >> x[3] >> p0 >> p[1] >> p[2] >> p[3]
       >> mass >> dummy >> i3 >> Q >> dummy >> dummy >> dummy >> dummy 
       >> dummy >> dummy >> q1 >> q2;
    //    G4cerr << j+1 << "  " << t << "  " << p0 << "  " << q1 << "  " << q2 << G4endl;
    int sig = sign(q2);
    RGB col(rand_gen::Random(1,3));
    switch ( abs(q2) ) {
    case 1 : iso = 0.5; break; 
    case 2 : iso = -0.5; break;
    default : iso = 0; break;
    }
    ParticleBase* p2 = makeParticle(*quarks[abs(q2)-1],QuantumProjections(sig,RGB((int)(sig*col)),sig*iso,-0.5+rand_gen::Random(0,1)));
    if ( abs(q1)<100 ) {
      switch ( abs(q1) ) {
      case 1 : iso = 0.5; break; 
      case 2 : iso = -0.5; break;
      default : iso = 0; break;
      }
      ParticleBase* p1 = makeParticle(*quarks[abs(q1)-1],QuantumProjections(-sig,RGB((int)(-sig*col)),sig*iso,-0.5+rand_gen::Random(0,1)));
      double s = sqr(mass);
      double mom = sqrt(sqr(s-sqr(p1->Mass())-sqr(p2->Mass()))-4*sqr(p1->Mass()*p2->Mass()))/(2*mass);
      p1->SetMomentum(Vektor3::isotropy(mom));
      p2->SetMomentum(-p1->Momentum());
      p1->SetCoordinates(x);
      p2->SetCoordinates(x);
      Vektor3 beta = -p/p0;
      p1->Lorentz(beta);
      p2->Lorentz(beta);
    }
    else {
      int q3 = int(abs(q1)/1000);
      int q4 = int((abs(q1) % 1000)/100);
      switch ( q3 ) {
      case 1 : iso = 0.5; break; 
      case 2 : iso = -0.5; break;
      default : iso = 0; break;
      }
      ParticleBase* p3 = makeParticle(*quarks[q3-1],QuantumProjections(sig,col>>1,sig*iso,-0.5+rand_gen::Random(0,1)));
      switch ( q4 ) {
      case 1 : iso = 0.5; break; 
      case 2 : iso = -0.5; break;
      default : iso = 0; break;
      }
      ParticleBase* p4 = makeParticle(*quarks[q4-1],QuantumProjections(sig,col>>2,sig*iso,-0.5+rand_gen::Random(0,1)));
      double s = sqr(mass);
      double m_di = 0.5;
      double mom = sqrt(sqr(s-sqr(p2->Mass())-sqr(m_di))-4*sqr(p2->Mass()*m_di))/(2*mass);
      double e_di = m_di;
      double mom1 = sqrt(sqr(e_di*e_di-sqr(p3->Mass())-sqr(p4->Mass()))-4*sqr(p3->Mass()*p4->Mass()))/(2*e_di);
      p2->SetMomentum(Vektor3::isotropy(mom));
      p3->SetMomentum(Vektor3::isotropy(mom1));
      p4->SetMomentum(-p3->Momentum());
      Vektor3 beta = -p/p0;
      Vektor3 beta1 = p2->Momentum()/sqrt(square(p2->Momentum())+e_di*e_di);
      p3->Lorentz(beta1);
      p4->Lorentz(beta1);
      p2->Lorentz(beta);
      p3->Lorentz(beta);
      p4->Lorentz(beta);
      p2->SetCoordinates(x);
      p3->SetCoordinates(x);
      p4->SetCoordinates(x);
    }
  }
  return t;
}
*/
double readUQMD(G4std::istream& in) 
{
  int N;
  double t,p0,mass,dummy,iso;
  int i3,Q,q1,q2;
  Vektor3 x,p;
  skipline(in,14);
  in >> N;
  G4cerr << "N=" << N << G4endl;
  skipline(in,2);
  for (int j=0; j<N; j++) {
    in >> t >> x[1] >> x[2] >> x[3] >> p0 >> p[1] >> p[2] >> p[3]
       >> mass >> dummy >> i3 >> Q >> dummy >> dummy >> dummy >> dummy 
       >> dummy >> dummy >> q1 >> q2;
    QuantumNumbers pp;
    double s = (q2 == 3) ? 1.0 : 0.0;
    double c = 1;
    double s3 = 0;
    if ( abs(q1)>1000 ) {
      if ( q1<0 ) c = -1;
      pp.setB(1);
      pp.setSpin(0.5);
      s3 = (rand_gen::Random()<0.5) ? 0.5 : -0.5;
      if ( (q1 % 1000)/100 == 3 ) s += 1.0;
      if ( q1/1000 == 3 ) s += 1.0;
    }
    else {
      pp.setB(0);
      pp.setSpin(0);
    }
    pp.setMaxColor(1);
    pp.setS(-s);
    pp.setIsospin(0.5*abs(i3));
    pp.setC(0);
    G4cerr << q1 << "  " << q2 << "  ";
    pp.writeOut(cerr);
    ParticleType& knot = Knot<ParticleType>::FindKnot((ParticleType&)pp);
    ParticleType& h = knot.selectType(mass);
    ParticleBase* H = makeParticle(h,QuantumProjections(c,0,0.5*i3,s3));
    G4cerr << H->Name() << G4endl;
    H->SetMomentum(p);
    H->SetCoordinates(x);
    H->SetMass(mass);
    CollisionTab::addEntry(0,*H,DECOMPOSITION);
  }
  //  CollisionTab::perform(0.001);
  skipline(in);
  return t;
}

class prof_tube : public Tube
{
  REAL tau,y,z0;
  double f(double z) const { return ( fabs(z)<double(z0) ) ? double(1.0/sqrt(tau*tau+z*z)) : 0.0; }
public:
  prof_tube(REAL y0,REAL t,REAL R) 
    : Tube(R,2*t*sinh(y0/2.0)),tau(t),y(y0),z0(tau*sinh(y0/2.0)) {}
  bool homogeneous(Vektor4&);
  bool isInside(const Vektor3& x) { return double(sqr(x[1])+sqr(x[2]))<=R*R && f(x[3]) > 0.0; }
};

bool prof_tube::homogeneous(Vektor4& x)
{
    REAL r = R*pow(rand_gen::Random(),0.5);
    REAL Phi = 2*mathConstants::Pi*rand_gen::Random();
    x[1] = r*cos(Phi);
    x[2] = r*sin(Phi);
    x[3] = sinh((rand_gen::Random()-0.5)*y)*tau;
  return true;
}

int main(int argc,char* argv[]) {
  StreamBuffer inputStream;
  inputStream.createStream(argc,argv);

  InputVariable<int> Number("n",1);
  InputVariable<int> nc("nc",0);
  InputVariable<String> dec("decay","yes");
  InputVariable<String> findec("final","yes");
  InputVariable<String> clu("cluster","yes");
  InputVariable<String> dir("direct","no");
  InputVariable<String> final_out("dir","/afs/cern.ch/user/s/sscherer/upp/data");
  InputVariable<double> deconf("deconf",0.01);
  InputVariable<double> kapp("kappa");
  InputVariable<double> Temp("T",0.150);
  InputVariable<double> T2("T::max");
  InputVariable<double> dT("T::step");
  InputVariable<double> Mu("Mu",0.0);
  InputVariable<double> Mu2("Mu::max");
  InputVariable<double> dMu("Mu::step");
  InputVariable<double> Mu_s("Mus",0.0);
  InputVariable<double> Rad("R",2.0);
  InputVariable<double> len("L",2.0);
  InputVariable<double> time("time",-1.0);
  InputVariable<double> dt("dt",1.0);
  InputVariable<double> step("h",0.05);
  InputVariable<double> ymax("ymax",0);
  InputVariable<double> Pt("pt",1);
  InputVariable<double> PtoN("p2n",1.0);
  InputVariable<int> M("number",200);
  InputVariable<int> skip("skip",0);
  InputVariable<String> Potential("potential","linear");
  InputVariable<String> Thermalize("thermalize","yes");
  InputVariable<String> Initial("initial","sphere");
  InputVariable<String> input("input");
  InputVariable<String> file("file");
  InputVariable<String> ForceDecay("forcedecay","no");
  InputVariable<double> transprof("transprof",0.0);
  
  InputReader ReadIn;
  ReadIn >> input >> Number >> Temp >> T2 >> dT >> Rad >> Mu >> Mu2 >> dMu 
	 >> Mu_s >> time >> dt >> ymax >> Potential >> Pt >> nc
	 >> step >> dec >> clu >> dir >> final_out >> findec >> Initial
	 >> deconf >> kapp >> len >> Thermalize >> M >> PtoN >> file
	 >> ForceDecay >> transprof >> skip;
  inputStream >> ReadIn;

  if ( input.isValid() ) {
    InputReader ReadIn;
    ReadIn >> input >> Number >> Temp >> T2 >> dT >> Rad >> Mu >> Mu2 >> dMu 
	   >> Mu_s >> time >> dt >> ymax >> Potential >> Pt >> nc
	   >> step >> dec >> clu >> dir >> final_out >> findec >> Initial
	   >> deconf >> kapp >> len >> Thermalize >> M >> PtoN >> file
	   >> ForceDecay >> transprof >> skip;
    inputStream >> ReadIn;
    G4std::istream* in;
    if ( input == "-" ) 
      in = &cin;
    else
      in = new G4std::ifstream((char*)String(input));
    *in >> ReadIn;
  }
  String pot = String(Potential);
  if ( kapp.isValid() ) 
    Colour::kappa = kapp;
  if ( pot == "interquark" ) 
    Colour::Pot = new InterQuark(Colour::kappa);
  else if ( pot == "cornell" ) 
    Colour::Pot = new Cornell;
  else 
    Colour::Pot = new Linear(Colour::kappa);

  //            Colour::Pot->print(cout,0,2,100);
  //            exit(0);
  String Dir = final_out;
  REAL T = REAL(Temp);
  REAL Time = REAL(time);
  REAL Dt = REAL(dt);
  REAL mu = REAL(Mu);
  REAL mus = REAL(Mu_s);
  double L =len;
  double pt = Pt;
  REAL R = REAL(Rad);
  int num = Number;
  bool therm = ( Thermalize == "yes" );
  bool finalDecay = ( findec == "yes" );
  bool forceDecay = ( ForceDecay == "yes" );
  Colour::allowDecay = (dec == "yes");
  Colour::allowClustering = (clu == "yes");
  Colour::directHadrons = (dir == "yes");
  Colour::deconfined = deconf;

  InverseFunction::reducePrecision = false;

  FileRead<ParticleType> Groups("/afs/cern.ch/user/s/sscherer/qmd_geant/test/Groups.dat");
  FileRead<ParticleType> Particles("/afs/cern.ch/user/s/sscherer/qmd_geant/test/Particles_Radiation.dat");
  Knot<ParticleType>::Root->printTree(cerr);
  FileRead<CollisionType> Collisions("/afs/cern.ch/user/s/sscherer/qmd_geant/test/Collisions.dat");

  double tau = 1.0;   // fm, formation time of QGP
  double tau_c = 0.1; // fm, formation tiem of c-cbar

  G4std::istream* readIn = 0;
  if ( file.isValid() ) 
     if ( file == "-" ) 
	readIn = &cin;
     else	
        readIn = new G4std::ifstream((char*)(String)file);
  G4cerr << "point 2\n";
  Geometry* Blob = 0;
  if ( Initial == "bjorken" )
    Blob = new Tube(R,2*tau);
	else if ( Initial == "tube" ) 
    Blob = new Tube(R,2*L);
  else if ( Initial == "sphere" ) 
    Blob = new Sphere(R);
  else if ( Initial == "cylinder" ) 
    Blob = new Tube(R,L);
  else
    throw "Unknown initial condition...";
  G4cerr << "point 1\n";
  Blob->whatAmI(cerr);
  //halfSpace Blob(R);
  //  Box Blob(R); // 3 fm Kantenlänge
  Output::fileout = new G4std::ofstream(Dir+"/output.out");
  *Output::fileout << "! Start\n";
  Output::fileout->flush();
  //  G4std::ostream& final = cerr;
  int n=0;
  G4cerr << "point 3\n";
  vector<ParticleType*> Quarks = Knot<ParticleType>::Find("SQ");
  Colour::setQuarks(Quarks);
  if ( readIn ) {
    for (int i=0; i<skip; i++) {
      skipline(*readIn);
      waitfor(*readIn,"UQMD");
    }
  }
  while ( ++n <= num ) {
    Colour box(step);
    //Radiation box(*Blob,step);
    //    inTheBox box(Blob,step);
    Particle::Soup = &box;
    
    const ParticleType& q = Knot<ParticleType>::FindKnot("q");
    const ParticleType& s = Knot<ParticleType>::FindKnot("s");

    //    QuarkVolume QGP(Quarks,Blob.getVolume(),T,mu,num);

    if ( readIn ) {
      String id;
      (*readIn) >> id;
      id.putToStream(*readIn);
      if ( id == "UQMD" ) 
	box.setTime(readUQMD(*readIn));
      else
	box.setTime(readEvent(*readIn));
    }
    else {
      QuarkVolume<RelBoltzmann> QGP(*Blob,Quarks,T,mu,mus,num);
      double u_frac = (2.0*double(PtoN)+1.0)/(PtoN+1.0)/3.0;
      QGP.createParticles(u_frac);
      G4cerr << "Energy : " << QGP.Etot()/Blob->getVolume() << "  " << R << "  " << box.List.size() << G4endl;
      N_c = 0;
      //    setPhaseSpace_stream(box.List,Blob,T,ymax);
      box.setup();
      setPhaseSpace(box.List,*Blob,T);
      if ( therm ) {
	QuarkBox run(*Blob,T,box);
	run.kstart = 20000;
	
	Colour::allowDecay = false;
	Colour::allowClustering = false; 
	G4cerr << "Thermalizing...";
	run.evaluate(21000);
	G4cerr << "finished\n";
	//      v_profile(box.List);
      }
      if ( Initial == "bjorken" ) {
	v_profile(box.List,tau);
	if ( transprof>0 ) 
	  trans_profile(box.List,R,double(transprof));
      }
      Colour::allowDecay = (dec == "yes");
      Colour::allowClustering = (clu == "yes");
      //    Colour::allowClustering = false;
      Colour::directHadrons = (dir == "yes");

      cout << "# Temperature " << T << G4endl;
      cout << "# mu_q        " << mu << G4endl;
      cout << "# mu_s        " << mus << G4endl;
      cout << "# pt        " << pt << G4endl;
      
      const ParticleType& cq = Knot<ParticleType>::FindKnot("c");
      
      //      for (int cc=0; cc<nc; cc++) {
      ParticleBase* p1 = makeParticle(q,QuantumProjections(1,RGB::RED,0,-0.5));
      ParticleBase* p2 = makeParticle(q,QuantumProjections(-1,-RGB::RED,0.5,0.5));
      ++N_c;
      //      double phi = 2*M_PI*rand_gen::Random();
      p1->SetMomentum(Vektor3(pt,pt,0));
      p2->SetMomentum(Vektor3(-pt,pt,0));
      box.setup();
      //	double rad = R*sqrt(1-sqrt(1-rand_gen::Random()));
      //	Vektor3 drad_c = p1->Momentum()/sqrt(pt*pt+sqr(p1->Mass()))*(tau-tau_c);
      //	phi = 2*M_PI*rand_gen::Random();
      //	Vektor3 e_r = Vektor3(rad*cos(phi),rad*sin(phi),0);
      //p1->SetCoordinates(e_r+drad_c);
      //p2->SetCoordinates(e_r-drad_c);
	//      }
    }
      

    try {
      box.print(cout);
      while ( ( box.Time() < Time || Time<0 ) && ( box.Nquark>0 || Time>0 ) ) {
	double t1 = box.Time()+dt;
	if ( Time>0 ) 
	  t1 = G4std::min(t1,(double)Time);
	//      box.Correlation();
	if ( box.Nquark ) {
	  while ( box.Time() < t1 && ( box.Nquark>0 || Time>0 ) ) {
	    //final.flush();
	    box.one_step();
	    G4cerr << n << " :  " << box.Time() << " : " << "  " 
		 << box.Etot() << "  "
	      //	       << "  " << "  " << length(box.Ptot()) << "  " 
		 << box.Npart << "  " << -box.Nquark+box.Npart << G4endl;
	  }
	}
	else {
	  box.setTime(t1);
	  CollisionTab::perform(t1);
	}
	box.print(cout);
	//      Output[0] << "# time=" << box.Time() << G4endl;
	//      box.writeField(Output[0],0,0,0,0,0,0);
      }
      if ( finalDecay ) {
	CollisionTab::perform(1000,forceDecay);
	cout << "# final\n";
	box.print(cout);
      }
      Output::fileout->flush();
    }
    catch ( char *s ) {
      G4cerr << "ERROR: " << s << G4endl;
    }
    int c = 1;
    //    double E0 = box.Etot(),E1 = E0;
    //  box.writeField(Feld,-10,10,0.4,-10,10,0.4);
    CollisionTab::erase();
  }
	return 1;
}

