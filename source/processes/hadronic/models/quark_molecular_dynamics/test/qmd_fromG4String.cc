#include "g4std/fstream"
#include "G4ios.hh"
#include <algo.h>
#include "newvector.hh"
#include "Random.hh"
#include "newBinning.hh"
#include "Geometry.hh"
#include "Arguments.hh"
#include "HadronGas.hh"
#include "MathTools.hh"
#include "Error.hh"
#include "HeatBath.hh"
#include "Propagation.hh"
#include "ParticleBase.hh"
#include "reactionChannels.hh"
#include "ParticleKinematics.hh"
#include "iso.hh"
#include "array.hh"
#include "Memory.hh"
#include "Collision.hh"
#include "genericRead.hh"
#include "StreamBuffer.hh"
#include "RunningVariable.hh"
#include "InputVariable.hh"
#include "InputReader.hh"
#include "Volume.hh"
#include "Potentials.hh"
#include "Quarkbox.hh"
#include "String.hh"
#include "output.hh"
#include "globals.hh"


ostream* Output::fileout = 0;

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

void skipline(istream& in,int n=1) {
  char c;
  for (int i=0; i<n; i++) {
     while ( in.get(c) && c != '\n' ) ;
  }
}

double readEvent(istream& in,bool final = true) 
{
  int S,C,col;
  double B,m,i,i3,s,s3,lt,time,flag;
  String name,checkstring;
  Vektor3 xs,ps,fs;

  bool readit = false;
  bool over = false;
  bool data = false;
  int newflag = 1;

  while ( in ) {
    in >> name;
// 
// line starts with # ?
//
    if ( name == "#" ) {
//
// if flag "data" was set (particle data have been read in),
// this means new timestep in list is reached -> STOP reading 
//
      if ( data ) {
       in.putback('#');
       break;
      }
//
// otherwise read on line, 
//           extract time, 
//           check whether "final" or "data" occurs and set flag readit,
//           read in rest of line and ignore it
//
      in >> checkstring;
      if ( (final && (checkstring.subString(0,4) == "final")) || (checkstring.subString(0,3) == "data") ) 
        readit = true;
      else if ( checkstring.subString(0,3) == "time" ) {
        time = atof(checkstring.subString(5,length(checkstring)-2));
        skipline(in);
      }
      else 
        skipline(in);
    }
//
// otherwise, line contais particle data, 
// and name is the particle name
// if flag readit is set, reading line of particle data,
// and set flag data (for reading/having read data)
//
    else if ( readit ) {
      data = true;
//
// read in particle data:
// 
// baryon strangenss charm mass isospin (i,i_z) spin (s,s_z) 
//
      in >> B   >> S       >> C  >> m  >> i >> i3     >> s >> s3;
//
// x-vector p-vector color  f-vector lifetime flag
//
      in >> xs    >> ps    >> col >> fs    >> lt    >> flag;
//
      skipline(in);
//
      if ( fabs(i3) > i ) { i3 = i; }
      if ( fabs(s3) > s ) { s3 = s; }
      double sig = sign(B);
      if ( sig == 0 ) sig = 1.0;
//
// strip of redundant anti- and 0 (for uncharged particles)
//
      int last;
      if ( name[length(name)-2] == '0' ) 
        last = length(name)-3;
      else 
        for (last=length(name)-2; last>0 && ( name[last] == '-' || name[last] == '+'); last--);
      name = name.subString(0,last);
      if ( name.subString(0,4) == "anti-" ) 
        name = name.subString(5);

      ParticleType& h = Knot<ParticleType>::FindKnot(name);
      ParticleBase* p = makeParticle(h,QuantumProjections(sig,col,i3,s3),ps,xs,m);
      p->setLifetime(lt);
      p->SetFlag(newflag);
      newflag++;
    }
//
// otherwise read in rest of line and ignore it ...
//
    else 
      skipline(in);
  }
//
// if data have been read in, flag data should have been set, 
//                            time should have been read in
//
  if ( data ) 
    return time;
  else
    return 0.0;
}



int main(int argc,char* argv[]) {

  StreamBuffer inputStream;
  inputStream.createStream(argc,argv);

  InputVariable<int> Number("n",1);
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
         >> Mu_s >> time >> dt >> ymax >> Potential 
         >> step >> dec >> clu >> dir >> final_out >> findec >> Initial
         >> deconf >> kapp >> len >> Thermalize >> M >> PtoN >> file
         >> ForceDecay >> transprof >> skip;
  inputStream >> ReadIn;

  if ( input.isValid() ) {
    InputReader ReadIn;
    ReadIn >> input >> Number >> Temp >> T2 >> dT >> Rad >> Mu >> Mu2 >> dMu 
           >> Mu_s >> time >> dt >> ymax >> Potential 
           >> step >> dec >> clu >> dir >> final_out >> findec >> Initial
           >> deconf >> kapp >> len >> Thermalize >> M >> PtoN >> file
           >> ForceDecay >> transprof >> skip;
    inputStream >> ReadIn;
    istream* in;
    if ( input == "-" ) 
      in = &cin;
    else
      in = new ifstream((char*)String(input));
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

  String Dir = final_out;
  REAL T = REAL(Temp);
  REAL Time = REAL(time);
  REAL Dt = REAL(dt);
  REAL mu = REAL(Mu);
  REAL mus = REAL(Mu_s);
  double L =len;
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
  Knot<ParticleType>::Root->printTree(G4cerr);
  FileRead<CollisionType> Collisions("/afs/cern.ch/user/s/sscherer/qmd_geant/test/Collisions.dat");

  double tau = 1.0;   // fm, formation time of QGP
  double tau_c = 0.1; // fm, formation tiem of c-cbar

  istream* readIn = 0;
  if ( file.isValid() ) 
    if ( file == "-" ) 
      readIn = &cin;
    else	
      readIn = new ifstream((char*)(String)file);



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

  Blob->whatAmI(G4cerr);


  Output::fileout = new ofstream(Dir+"/output.out");
  *Output::fileout << "! Start\n";
  Output::fileout->flush();

  int n=0;

  vector<ParticleType*> Quarks = Knot<ParticleType>::Find("SQ");
  Colour::setQuarks(Quarks);

  while ( ++n <= num ) {

    Colour box(step);
    Particle::Soup = &box;

    const ParticleType& q = Knot<ParticleType>::FindKnot("q");
    const ParticleType& s = Knot<ParticleType>::FindKnot("s");

    if ( readIn ) {

      box.setTime(readEvent(*readIn));
      box.setup();

    }
    else {

      QuarkVolume<RelBoltzmann> QGP(*Blob,Quarks,T,mu,mus,num);
      double u_frac = (2.0*double(PtoN)+1.0)/(PtoN+1.0)/3.0;
      QGP.createParticles(u_frac);
      G4cerr << " - Energy : " << QGP.Etot()/Blob->getVolume() << "  " << R << "  " << box.List.size() << endl;
  
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
      }
  
      if ( Initial == "bjorken" ) {
        v_profile(box.List,tau);
        if ( transprof>0 ) 
          trans_profile(box.List,R,double(transprof));
      }
  
      Colour::allowDecay = (dec == "yes");
      Colour::allowClustering = (clu == "yes");
      Colour::directHadrons = (dir == "yes");
  
      cout << "# Temperature " << T << endl;
      cout << "#        mu_q " << mu << endl;
      cout << "#        mu_s " << mus << endl;
        
    }

    try {
      box.print(cout);
      while ( ((Time < 0) && (box.Nquark > 0)) || (box.Time() < Time) ) {

        double t1 = box.Time()+dt;
        if ( Time>0 ) 
          t1 = G4std::min(t1,(double)Time);

        if ( box.Nquark ) {
          while ( box.Time() < t1 && ( box.Nquark>0 || Time>0 ) ) {
            box.one_step();
            G4cerr << n << " :  " << box.Time() << " : " << "  " 
                   << box.Etot() << "  "
                   << box.Npart << "  " << -box.Nquark+box.Npart << endl;
          }
        }
        else {
          box.setTime(t1);
          CollisionTab::perform(t1);
        }
        box.print(cout);
      }
      if ( finalDecay ) {
        CollisionTab::perform(1000,forceDecay);
        cout << "# final\n";
        box.print(cout);
      }
      Output::fileout->flush();
    }
    catch ( char *s ) {
      G4cerr << "ERROR: " << s << endl;
    }
    int c = 1;
  }

  return 1;
}

