//
// short version oc qmd_fromfile.cc
//
// - read in initial date from qmd-data file
// - process data
//
// Stefan Scherer, 2000-08-02
//
//
//


#include "g4std/fstream"
#include "G4ios.hh"
#include <algo.h>
#include "newvector.hh"
#include "Random.hh"
#include "newBinning.hh"
#include "Geometry.hh"
#include "Arguments.hh"
#include "MathTools.hh"
#include "Error.hh"
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

void skipline(istream& in,int n=1) 
{
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
// otherwise, line contains particle data, 
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

  InputVariable<String> file("file");
  InputVariable<String> dec("decay","yes");
  InputVariable<String> clu("cluster","yes");
  InputVariable<String> dir("direct","no");
  InputVariable<String> findec("final","yes");
  InputVariable<String> ForceDecay("forcedecay","no");
  InputVariable<String> Potential("potential","linear");
  InputVariable<double> deconf("deconf",0.01);
  InputVariable<double> kapp("kappa");
  InputVariable<double> time("time",-1.0);
  InputVariable<double> dt("dt",1.0);
  InputVariable<double> step("h",0.05);
  

  InputReader ReadIn;

  ReadIn >> file
         >> dec >> clu >> dir >> findec >> ForceDecay
         >> Potential >> deconf >> kapp 
         >> time >> dt >> step;
  inputStream >> ReadIn;

  istream* readIn = 0;
  if ( file.isValid() ) 
    if ( file == "-" ) 
      readIn = &cin;
    else	
      readIn = new ifstream((char*)(String)file);


  if ( kapp.isValid() ) 
    Colour::kappa = kapp;

  String pot = String(Potential);
  if ( pot == "interquark" ) 
    Colour::Pot = new InterQuark(Colour::kappa);
  else if ( pot == "cornell" ) 
    Colour::Pot = new Cornell;
  else 
    Colour::Pot = new Linear(Colour::kappa);

  bool finalDecay = ( findec == "yes" );
  bool forceDecay = ( ForceDecay == "yes" );
  Colour::allowDecay = (dec == "yes");
  Colour::allowClustering = (clu == "yes");
  Colour::directHadrons = (dir == "yes");
  Colour::deconfined = deconf;

  REAL Time = REAL(time);
  REAL Dt = REAL(dt);

  InverseFunction::reducePrecision = false;

  FileRead<ParticleType> Groups("/afs/cern.ch/user/s/sscherer/qmd_geant/test/Groups.dat");
  FileRead<ParticleType> Particles("/afs/cern.ch/user/s/sscherer/qmd_geant/test/Particles_Radiation.dat");
  Knot<ParticleType>::Root->printTree(G4cerr);
  FileRead<CollisionType> Collisions("/afs/cern.ch/user/s/sscherer/qmd_geant/test/Collisions.dat");

  vector<ParticleType*> Quarks = Knot<ParticleType>::Find("SQ");
  Colour::setQuarks(Quarks);

  Colour box(step);
  Particle::Soup = &box;

  const ParticleType& q = Knot<ParticleType>::FindKnot("q");
  const ParticleType& s = Knot<ParticleType>::FindKnot("s");

  if ( readIn ) {

    box.setTime(readEvent(*readIn));
    box.setup();

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
          G4cerr << box.Time() << " : " << "  " 
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

  }
  catch ( char *s ) {
    G4cerr << "ERROR: " << s << endl;
  }

  return 1;

}

