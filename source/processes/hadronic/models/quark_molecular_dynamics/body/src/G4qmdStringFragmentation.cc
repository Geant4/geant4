// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// G4qmdStringFragmentation.cc  2001/03/29  Stefan Scherer
//
#include "G4qmdStringFragmentation.hh"
#include "G4ios.hh"
#include "g4std/fstream"
#include "globals.hh"
#include "G4KineticTrack.hh"
#include "G4KineticTrackVector.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"

//
// ----------- constructor ------------------------------
//
G4qmdStringFragmentation::G4qmdStringFragmentation()
{
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

	G4cerr << "### Object of type G4qmdStringFragmentation created successfully ###"  << G4endl;

}

G4qmdStringFragmentation::G4qmdStringFragmentation(const G4qmdStringFragmentation &right)
{
}

G4qmdStringFragmentation::~G4qmdStringFragmentation()
{
}

const G4qmdStringFragmentation & G4qmdStringFragmentation::operator=(const G4qmdStringFragmentation &right)
{
  G4Exception("G4qmdStringFragmentation::operator= meant to not be accessable");
  return *this;
}

int G4qmdStringFragmentation::operator==(const G4qmdStringFragmentation &right) const
{
  return 0;
}

int G4qmdStringFragmentation::operator!=(const G4qmdStringFragmentation &right) const
{
  return 1;
}
//
// --------------------------------------------
//

G4KineticTrackVector * FragmentStrings(const G4ExcitedStringVector * theStrings)
{
	SetupFromExcitedStringVector(theStrings);
	theResult = TheHadrons();
	return theResult;
}

//
// ----------- Setup from Excited String ---------------------------------
//
void G4qmdStringFragmentation::SetupFromExcitedStringVector(const G4ExcitedStringVector * theStrings)
{
  
  Particle::Soup = &theQuarkSystem;

  const ParticleType& q = Knot<ParticleType>::FindKnot("q");    // u|d-quark
  const ParticleType& s = Knot<ParticleType>::FindKnot("s");    // s-quark
  const ParticleType& cq = Knot<ParticleType>::FindKnot("c");   // c-quark

	for (G4int StringCounter=0; StringCounter < theStrings->length(); StringCounter++) {
 
		G4ExcitedString* ThisExcitedString = theStrings->at(StringCounter);
		const G4PartonVector* ThePartons = ThisExcitedString->GetPartonList();

		G4cerr << "Extract partons from G4ExcitedString " << StringCounter << endl;

	  for (G4int QuarkCounter=0; QuarkCounter < ThePartons->length(); QuarkCounter++) {

			G4int flag = StringCounter*1000 + QuarkCounter;

 		  G4Parton* ThisParton = ThePartons->at(QuarkCounter);

      const G4ThreeVector quark_ThreeVector = ThisParton->GetPosition();
      const G4LorentzVector quark_4Momentum = ThisParton->Get4Momentum(); 
			G4int quark_PDGCode = ThisParton->GetPDGcode();
      G4int quark_Colour = ThisParton->GetColour();
      G4double quark_SpinZ = ThisParton->GetSpinZ();
      G4double quark_IsoSpinZ = ThisParton->GetIsoSpinZ();

      Vektor3 quark_momentum = Vektor3(quark_4MomentumVector->x(),quark_4MomentumVector->y(),quark_4MomentumVector->z());
      Vektor3 quark_position = Vektor3(quark_ThreeVector->x(),quark_ThreeVector->y(),quark_ThreeVector->z());

			G4cerr << " Parton " << flag << " with properties " << G4endl;
      G4cerr << "    PDG-CODE:  " <<  quark_PDGCode << G4endl; 
      G4cerr << "  x-3-Vector:  " <<  quark_ThreeVector << G4endl; 
      G4cerr << "  p-4-Vector:  " <<  quark_4Momentum << G4endl; 
      G4cerr << "       color:  " <<  quark_Colour << G4endl; 
      G4cerr << "      Spin_Z:  " <<  quark_SpinZ << G4endl ;
      G4cerr << "   Isospin_Z:  " <<  quark_IsoSpinZ << G4endl;
      G4cerr << "  3_momentum:  " <<  quark_momentum << G4endl;
      G4cerr << "  3_position:  " <<  quark_position << G4endl;

//
//  everything still missing ...
//

//  general scheme: sig=1 - particle, sig=-1 antiparticle
//
//    ParticleBase* p = makeParticle(q|s|cq, QuantumProjections(sig,col,i3,s3));

//      ParticleBase* p1 = makeParticle(cq,QuantumProjections(1,RGB::RED,quark_IsoSpinZ,quark_SpinZ));  // c-quark
//      p1->SetMomentum(quark_momentum);
//      p1->SetCoordinates(quark_position);
//      p1->SetFlag(i*10000+j);

//    ParticleBase* p2 = makeParticle(q,QuantumProjections(-1,-RGB::RED,0.5,0.5)); // u|d-quark
//    p2->SetMomentum(Vektor3(px,py,pz));
//    p2->SetCoordinates(Vektor3(rx,ry,rz));
//    p2->SetFlag(i*10000+j);


		}

	}

}


G4KineticTrackVector * G4qmdStringFragmentation::TheHadrons()
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

	for (G4int i=0; i<theHadrons->length(); i++) {
 		G4KineticTrack * ThisTrack = (*theHadrons)[i];
		G4cerr << "Track "      << i 
		       << ", particle " << ThisTrack->GetDefinition()->GetPDGEncoding()
		       << " originating at " << ThisTrack->GetInitialCoordinates() 
		       << endl;
	}



	return theHadrons;

}

//
// --------------------------------------------
//

void G4qmdStringFragmentation::skipline(istream& in) 
{
  char c;
  while ( in.get(c) && c != '\n' ) ;
}

double G4qmdStringFragmentation::readEvent(istream& in) 
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


void G4qmdStringFragmentation::SetupFromFile(const G4String & anInputFile="")
{
  if (anInputFile == "") {
    theInputFile = "/afs/cern.ch/user/s/sscherer/public/qmd/data/line_120.dat";
  }
  else {
    theInputFile = anInputFile;
  }

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


