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
  InputVariable<String> ForceDecay("forcedecay","no");
  InputVariable<double> transprof("transprof",0.0);
  
  InputReader ReadIn;
  ReadIn >> input >> Number >> Temp >> T2 >> dT >> Rad >> Mu >> Mu2 >> dMu 
   >> Mu_s >> time >> dt >> ymax >> Potential >> Pt >> nc
   >> step >> dec >> clu >> dir >> final_out >> findec >> Initial
   >> deconf >> kapp >> len >> Thermalize >> M >> PtoN 
   >> ForceDecay >> transprof >> skip;
  inputStream >> ReadIn;

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

  Blob->whatAmI(cerr);

  Output::fileout = new G4std::ofstream(Dir+"/output.out");
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

    QuarkVolume<RelBoltzmann> QGP(*Blob,Quarks,T,mu,mus,num);

    double u_frac = (2.0*double(PtoN)+1.0)/(PtoN+1.0)/3.0;

    Colour::allowDecay = (dec == "yes");
    Colour::allowClustering = (clu == "yes");
    Colour::directHadrons = (dir == "yes");

    cout << "# Temperature " << T << G4endl;
    cout << "#        mu_q " << mu << G4endl;
         

    ParticleBase* p1 = makeParticle(q,QuantumProjections(1,RGB::RED,0,-0.5));
    ParticleBase* p2 = makeParticle(q,QuantumProjections(-1,-RGB::RED,0.5,0.5));
    ParticleBase* p4 = makeParticle(q,QuantumProjections(1,RGB::BLUE,0,-0.5));
    ParticleBase* p3 = makeParticle(q,QuantumProjections(-1,-RGB::BLUE,0.5,0.5));


    Vektor3 e_x = Vektor3(1,0,0);
    Vektor3 e_y = Vektor3(0,1,0);
    Vektor3 e_z = Vektor3(0,0,1);

    p1->SetMomentum(3*e_x+2*e_z);
    p2->SetMomentum(-2*e_x+1*e_z);
    p3->SetMomentum(1*e_z+e_y);
    p4->SetMomentum(3*e_x-1*e_z-e_y);
    
    p1->SetCoordinates(e_y-4*e_z);
    p2->SetCoordinates(-3*e_y-4*e_z);
    p3->SetCoordinates(-2*e_y+4*e_z);
    p4->SetCoordinates(-1*e_z);
	

    box.setup();
 
 
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
  	             << box.Npart << "  " << -box.Nquark+box.Npart << G4endl;
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
      G4cerr << "ERROR: " << s << G4endl;
    }
    int c = 1;
  }

  return 1;
}

