// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// G4qmd.cc  2000/08/03  Stefan Scherer
//
#include "G4qmd.hh"
#include "G4ios.hh"
#include "g4std/fstream"
#include "globals.hh"
#include "G4KineticTrack.hh"
#include "G4KineticTrackVector.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"

G4qmd::G4qmd()
{
}

G4qmd::G4qmd(const G4String & anInputFile)
{
  if (anInputFile == "") {
    theInputFile = "/afs/cern.ch/user/s/sscherer/public/qmd/data/line_120.dat";
  }
  else {
    theInputFile = anInputFile;
  }
  theColorStringDecay = "yes";
  theColorCluster = "yes";
  theDirectHadronFromString = "no";
  theFinalHadronDecay = "yes";
  theForcedHadronDecay = "no";
  theColorPotential = "linear";
  theParameterHadronizationCriterium = 0.01;
  theParameterKappa = 1.8;
  theInternalTimestep = 0.05;
  theFinalTime = -1.0;
  theOutputTimestep = 1.0;

  Colour::kappa = theParameterKappa;
  Colour::deconfined = theParameterHadronizationCriterium;
  Colour::allowDecay = (theColorStringDecay == "yes");
  Colour::allowClustering = (theColorCluster == "yes");
  Colour::directHadrons = (theDirectHadronFromString == "yes");

  if ( theColorPotential == "interquark" ) 
    Colour::Pot = new InterQuark(Colour::kappa);
  else if ( theColorPotential == "cornell" ) 
    Colour::Pot = new Cornell;
  else 
    Colour::Pot = new Linear(Colour::kappa);

  InverseFunction::reducePrecision = false;

  FileRead<ParticleType> Groups("/afs/cern.ch/user/s/sscherer/qmd_geant/test/Groups.dat");
  FileRead<ParticleType> Particles("/afs/cern.ch/user/s/sscherer/qmd_geant/test/Particles_Radiation.dat");
  Knot<ParticleType>::Root->printTree(G4cerr);
  FileRead<CollisionType> Collisions("/afs/cern.ch/user/s/sscherer/qmd_geant/test/Collisions.dat");
  vector<ParticleType*> Quarks = Knot<ParticleType>::Find("SQ");
  Colour::setQuarks(Quarks);

  Colour theQuarkSystem(theInternalTimestep);

}

G4qmd::G4qmd(const G4qmd &right)
{
}

G4qmd::~G4qmd()
{
}

const G4qmd & G4qmd::operator=(const G4qmd &right)
{
  G4Exception("G4qmd::operator= meant to not be accessable");
  return *this;
}

int G4qmd::operator==(const G4qmd &right) const
{
  return 0;
}

int G4qmd::operator!=(const G4qmd &right) const
{
  return 1;
}

//
// --------------------------------------------
//

void G4qmd::skipline(istream& in) 
{
  char c;
  while ( in.get(c) && c != '\n' ) ;
}

double G4qmd::readEvent(istream& in) 
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
//           check "data" occurs and set flag readit,
//           read in rest of line and ignore it
//
      in >> checkstring;
      if (checkstring.subString(0,3) == "data") 
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


void G4qmd::SetupFromFile()
{

  G4std::istream* readIn = 0;
  readIn = new G4std::ifstream(theInputFile);


  Particle::Soup = &theQuarkSystem;

  const ParticleType& q = Knot<ParticleType>::FindKnot("q");
  const ParticleType& s = Knot<ParticleType>::FindKnot("s");

  if ( readIn ) {

    theQuarkSystem.setTime(readEvent(*readIn));
    theQuarkSystem.setup();

  }

}



G4KineticTrackVector * G4qmd::TheHadrons()
{
	G4KineticTrackVector * theHadrons = new G4KineticTrackVector();
  try {

    theQuarkSystem.print(cout);
 
    while ( ((theFinalTime < 0) && (theQuarkSystem.Nquark > 0)) || (theQuarkSystem.Time() < theFinalTime) ) {

      double t1 = theQuarkSystem.Time()+theOutputTimestep;
      if ( theFinalTime>0 ) 
        t1 = G4std::min(t1,(double)theFinalTime);

      if ( theQuarkSystem.Nquark ) {
        while ( theQuarkSystem.Time() < t1 && ( theQuarkSystem.Nquark>0 || theFinalTime>0 ) ) {
          theQuarkSystem.one_step();
          G4cerr << theQuarkSystem.Time() << " : " << "  " 
                 << theQuarkSystem.Etot() << "  "
                 << theQuarkSystem.Npart << "  " << -theQuarkSystem.Nquark+theQuarkSystem.Npart << endl;
        }
      }
      else {
        theQuarkSystem.setTime(t1);
        CollisionTab::perform(t1);
      }
      theQuarkSystem.print(cout);
    }
    if ( theFinalHadronDecay == "yes" ) {
      CollisionTab::perform(1000,(theForcedHadronDecay == "yes"));
      cout << "# final\n";
      theQuarkSystem.print(cout);
    }

		theHadrons = theQuarkSystem.GetNewHadrons();

  }
  catch ( char *s ) {
    G4cerr << "ERROR: " << s << G4endl;
  }

	return theHadrons;

}

