//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// G4qmdStringFragmentation.cc  2001/03/29  Stefan Scherer
//
#include "G4qmdStringFragmentation.hh"
#include "G4ios.hh"
#include "g4std/fstream"
#include "globals.hh"

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

  theQuarkSystem = new Colour(theInternalTimestep);

	G4cerr << "### Object of type G4qmdStringFragmentation created successfully ###"  << G4endl;

}

G4qmdStringFragmentation::G4qmdStringFragmentation(const G4qmdStringFragmentation &right)
{
}

G4qmdStringFragmentation::~G4qmdStringFragmentation()
{
  delete theQuarkSystem;
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

G4KineticTrackVector * G4qmdStringFragmentation::FragmentStrings(const G4ExcitedStringVector * theStrings)
{
	G4KineticTrackVector * theResult = new G4KineticTrackVector();
	SetupFromG4ExcitedStringVector(theStrings);
	theResult = TheHadrons();
	return theResult;
}


G4KineticTrackVector * G4qmdStringFragmentation::FragmentStringsFromFile(const G4String & anInputFile)
{
	G4KineticTrackVector * theResult = new G4KineticTrackVector();
	SetupFromFile(anInputFile);
	theResult = TheHadrons();
	return theResult;
}


//
// ----------- Setup from Excited String ---------------------------------
//

void G4qmdStringFragmentation::SetupFromG4ExcitedStringVector(const G4ExcitedStringVector * theStrings)
{
  
  Particle::Soup = theQuarkSystem;

  const ParticleType& q = Knot<ParticleType>::FindKnot("q");    // u|d-quark
  const ParticleType& s = Knot<ParticleType>::FindKnot("s");    // s-quark
  const ParticleType& cq = Knot<ParticleType>::FindKnot("c");   // c-quark

	for (G4int StringCounter=0; StringCounter < theStrings->length(); StringCounter++) {
 
		G4ExcitedString* ThisExcitedString = theStrings->at(StringCounter);
		const G4PartonVector* ThePartons = ThisExcitedString->GetPartonList();

		G4cerr << "Extract partons from G4ExcitedString " << StringCounter << endl;

	  for (G4int PartonCounter=0; PartonCounter < ThePartons->length(); PartonCounter++) {

			G4int flag = (StringCounter*1000 + PartonCounter)*10;

 		  G4Parton* ThisParton = ThePartons->at(PartonCounter);

			const G4String& parton_type = ThisParton->GetDefinition()->GetParticleType();
			G4int parton_PDGCode = ThisParton->GetPDGcode();
      G4int parton_Colour = ThisParton->GetColour();
      G4double parton_Mass = ThisParton->GetDefinition()->GetPDGMass();
      G4double parton_SpinZ = ThisParton->GetSpinZ();
      G4double parton_IsoSpinZ = ThisParton->GetIsoSpinZ();
      const G4ThreeVector  parton_ThreeVector = ThisParton->GetPosition();
      const G4LorentzVector  parton_4Momentum = ThisParton->Get4Momentum(); 

      Vektor3 parton_momentum =     0.0003*Vektor3(parton_4Momentum.x(),parton_4Momentum.y(),parton_4Momentum.z());
      Vektor3 parton_position =    1.0e+12*Vektor3(parton_ThreeVector.x(),parton_ThreeVector.y(),parton_ThreeVector.z());

      Vektor3 parton_fudgediff =    Vektor3(0.5,0.5,0);

			G4cerr << " Parton " << flag << " with properties " << G4endl;
      G4cerr << "        Type:  " <<  parton_type << G4endl; 
      G4cerr << "    PDG-CODE:  " <<  parton_PDGCode << G4endl; 
      G4cerr << "       color:  " <<  parton_Colour << G4endl; 
      G4cerr << "      Spin_Z:  " <<  parton_SpinZ << G4endl ;
      G4cerr << "   Isospin_Z:  " <<  parton_IsoSpinZ << G4endl;
      G4cerr << "  x-3-Vector:  " <<  parton_ThreeVector << G4endl; 
      G4cerr << "  p-4-Vector:  " <<  parton_4Momentum << G4endl; 
      G4cerr << "  3_momentum:  " <<  parton_momentum << G4endl;
      G4cerr << "  3_position:  " <<  parton_position << G4endl;

      G4cerr << "  ... now checking in:  " << G4endl;

      if (parton_type == "quarks") {
        if (abs(parton_PDGCode) < 3) {
 		     	ParticleBase* p = makeParticle(q,QuantumProjections(sign(parton_Colour),parton_Colour,parton_IsoSpinZ,parton_SpinZ));  
	 	      p->SetMomentum(parton_momentum);
    	    p->SetCoordinates(parton_position);
      	  p->SetFlag(flag);
	        G4cerr << "(anti) u/d-quark created!" << G4endl;
			  }
  			else if ((parton_PDGCode = 3) || (parton_PDGCode = -3)) {
        	ParticleBase* p = makeParticle(s,QuantumProjections(sign(parton_Colour),parton_Colour,parton_IsoSpinZ,parton_SpinZ));  
	        p->SetMomentum(parton_momentum);
  	      p->SetCoordinates(parton_position);
    	    p->SetFlag(flag);
  	      G4cerr << "(anti) s-quark created!" << G4endl;
	  		}
		  	else if ((parton_PDGCode = 4) || (parton_PDGCode = -4)) {
        	ParticleBase* p = makeParticle(cq,QuantumProjections(sign(parton_Colour),parton_Colour,parton_IsoSpinZ,parton_SpinZ));  
  	      p->SetMomentum(parton_momentum);
    	    p->SetCoordinates(parton_position);
      	  p->SetFlag(flag);
	        G4cerr << "(anti) c-quark created!" << G4endl;
			  }
  			else {
	        G4cerr << "3th generation (b or t) quark - ignored!" << G4endl;
		  	}
      }
			else if (parton_type == "gluons") {
	      G4cerr << "gluon - ignored!" << G4endl;
			}
			else if (parton_type == "diquarks") {
				if (abs(parton_PDGCode) < 3000) {
					G4int isDiquark = sign(parton_PDGCode) ;
					// 
					// Isospin from PDG-Code (for u/d-quarks only)
					// 
					G4int Quark1_absPDGCode =  abs(parton_PDGCode)/1000;
					G4int Quark2_absPDGCode =  (abs(parton_PDGCode) - Quark1_absPDGCode*1000)/100 ;
					G4double Quark1_IsoSpin3 = Quark1_absPDGCode - 1.5;
					G4double Quark2_IsoSpin3 = Quark2_absPDGCode - 1.5;
					G4int Quark1_absColour, Quark2_absColour;
					// 
					// Colour from colour of diquark
					// 
					if (abs(parton_Colour)==1) { 
						if (G4UniformRand() > 0.5) {
							Quark1_absColour = 2; 
							Quark2_absColour = 3; 
						}
						else {
							Quark1_absColour = 3; 
							Quark2_absColour = 2; 
						} 
					}
					else if (abs(parton_Colour)==2) { 
						if (G4UniformRand() > 0.5) {
							Quark1_absColour = 1; 
							Quark2_absColour = 3; 
						}
						else {
							Quark1_absColour = 3; 
							Quark2_absColour = 1; 
						} 
					}
					else if (abs(parton_Colour)==3) { 
						if (G4UniformRand() > 0.5) {
							Quark1_absColour = 1; 
							Quark2_absColour = 2; 
						}
						else {
							Quark1_absColour = 2; 
							Quark2_absColour = 1; 
						} 
					}
					else {
			      G4cerr << "ERRROR: Colour " << parton_Colour << " is not a valid diquark colour!" << G4endl;
					}
					// 
					// spin-z from spin-z of diquark
					// 
					G4double Quark1_Spin3, Quark2_Spin3;
					if (parton_SpinZ==0.0) { 
						if (G4UniformRand() > 0.5) {
							Quark1_Spin3 =  0.5; 
							Quark2_Spin3 = -0.5;
						}
						else {
							Quark1_Spin3 = -0.5; 
							Quark2_Spin3 =  0.5;
						} 
					}
					else {
						if (G4UniformRand() > 0.5) {
							Quark1_Spin3 =  0.5; 
							Quark2_Spin3 =  0.5;
						}
						else {
							Quark1_Spin3 =  -0.5; 
							Quark2_Spin3 =  -0.5;
						} 
					}
//
// ------ CAUTION: qMD-to-Geant4-fittig --------
//
//  Masses are in GeV in qMD, in MeV in Geant4					
//  Masses of u,d quarks are EQUAL in qMD:             m_u = 5 MeV, m_d =  5 MeV
//  Masses of u,d quarks are are DIFFERENT in Geant4:  m_u = 8 MeV, m_d = 15 MeV
//
//  3momentum-distribution on the two quarks should consider this!
//
// ---------------------------------------------
					G4double Quark1_Mass = 5;
					G4double Quark2_Mass = 5;
					Vektor3 Quark1_3momentum = (Quark1_Mass/parton_Mass)* parton_momentum;
					Vektor3 Quark2_3momentum = (Quark2_Mass/parton_Mass)* parton_momentum;

        	ParticleBase* q1 = makeParticle(q,QuantumProjections(isDiquark,isDiquark*Quark1_absColour,isDiquark*Quark1_IsoSpin3,Quark1_Spin3));  
  	      q1->SetMomentum(Quark1_3momentum+parton_fudgediff);
    	    q1->SetCoordinates(parton_position);
      	  q1->SetFlag(flag);
	        G4cerr << "u/d-quark created from (anti)-diquark ("  << isDiquark << ")"
					       << " colour: " << isDiquark*Quark1_absColour
					       << " isospin-3: " << isDiquark*Quark1_IsoSpin3
					       << " spin-3: " << Quark1_Spin3
					       << G4endl;

        	ParticleBase* q2 = makeParticle(q,QuantumProjections(isDiquark,isDiquark*Quark2_absColour,isDiquark*Quark2_IsoSpin3,Quark2_Spin3));  
  	      q2->SetMomentum(Quark2_3momentum-parton_fudgediff);
    	    q2->SetCoordinates(parton_position);
      	  q2->SetFlag(flag+1);
	        G4cerr << "u/d-quark created from (anti)-diquark ("  << isDiquark << ")"
					       << " colour: " << isDiquark*Quark2_absColour
					       << " isospin-3: " << isDiquark*Quark2_IsoSpin3
					       << " spin-3: " << Quark2_Spin3
					       << G4endl;
				}
	      else {
					G4cerr << "Diquark with strangeness or heavy flavour - ignored!" << G4endl;
				}
			}
  		else {
	      G4cerr << "ERRROR: Particle with PDGCode " << parton_PDGCode << " is not a vild parton!" << G4endl;
      }
		}

	}
  
  theQuarkSystem->setup();

}

//
// ----------- Setup from file ---------------------------------
//

void G4qmdStringFragmentation::SetupFromFile(const G4String & anInputFile)
{
  if (anInputFile == "") {
    theInputFile = "/afs/cern.ch/user/s/sscherer/public/qmd/data/line_120.dat";
  }
  else {
    theInputFile = anInputFile;
  }

  G4std::istream* readIn = 0;
  readIn = new G4std::ifstream(theInputFile);


  Particle::Soup = theQuarkSystem;

  const ParticleType& q = Knot<ParticleType>::FindKnot("q");    // u|d-quark
  const ParticleType& s = Knot<ParticleType>::FindKnot("s");    // s-quark
  const ParticleType& cq = Knot<ParticleType>::FindKnot("c");   // c-quark

  if ( readIn ) {

    theQuarkSystem->setTime(readEvent(*readIn));

  }

	theQuarkSystem->setup();

}


//
// ----------- Run the code ---------------------------------
//

G4KineticTrackVector * G4qmdStringFragmentation::TheHadrons()
{
	G4KineticTrackVector * theHadrons = new G4KineticTrackVector();
  try {

    theQuarkSystem->print(cout);
 
    while ( ((theFinalTime < 0) && (theQuarkSystem->Nquark > 0)) || (theQuarkSystem->Time() < theFinalTime) ) {

      double t1 = theQuarkSystem->Time()+theOutputTimestep;
      if ( theFinalTime>0 ) 
        t1 = G4std::min(t1,(double)theFinalTime);

      if ( theQuarkSystem->Nquark ) {
        while ( theQuarkSystem->Time() < t1 && ( theQuarkSystem->Nquark>0 || theFinalTime>0 ) ) {
          theQuarkSystem->one_step();
          G4cerr << theQuarkSystem->Time() << " : " << "  " 
                 << theQuarkSystem->Etot() << "  "
                 << theQuarkSystem->Npart << "  " << -theQuarkSystem->Nquark+theQuarkSystem->Npart << endl;
        }
      }
      else {
        theQuarkSystem->setTime(t1);
        CollisionTab::perform(t1);
      }
      theQuarkSystem->print(cout);
    }
    if ( theFinalHadronDecay == "yes" ) {
      CollisionTab::perform(1000,(theForcedHadronDecay == "yes"));
      cout << "# final\n";
      theQuarkSystem->print(cout);
    }

		theHadrons = theQuarkSystem->GetNewHadrons();

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


