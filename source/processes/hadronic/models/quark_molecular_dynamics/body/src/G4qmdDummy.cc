// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// G4qmdDummy.cc  2000/08/03  Stefan Scherer
//
#include "G4qmdDummy.hh"
#include "G4ios.hh"
#include "g4std/fstream"
#include "globals.hh"
#include "String.hh"

G4qmdDummy::G4qmdDummy()
{
}

G4qmdDummy::G4qmdDummy(const G4String & anInputFile)
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
  theFinalTime = -1.0;
  theOutputTimestep = 1.0;
  theInternalTimestep = 0.05;

  Colour theBox(theInternalTimestep);

	G4cout << "... Object of class G4qmdDummy initalized with file " << theInputFile << G4endl;

}

G4qmdDummy::G4qmdDummy(const G4qmdDummy &right)
{
}

G4qmdDummy::~G4qmdDummy()
{
}

const G4qmdDummy & G4qmdDummy::operator=(const G4qmdDummy &right)
{
  G4Exception("G4qmdDummy::operator= meant to not be accessable");
  return *this;
}

int G4qmdDummy::operator==(const G4qmdDummy &right) const
{
  return 0;
}

int G4qmdDummy::operator!=(const G4qmdDummy &right) const
{
  return 1;
}

//
// --------------------------------------------
//

void G4qmdDummy::skipline(istream& in) 
{
	G4cout << "... read in line" << G4endl;
  char c;
  while ( in.get(c) && c != '\n' ) ;
}

double G4qmdDummy::readEvent(istream& in) 
{
  G4String name,checkstring;

  while ( in ) {
    in >> name;

    if ( name == "#" ) {
      in >> checkstring;
    }
   G4cout << checkstring;
  }

  G4cout << "... read in input file" << G4endl;

  return 0;
}


void G4qmdDummy::SetupFromFile()
{
  G4cout << "... set up quark system" << G4endl;
}



void G4qmdDummy::justRun()
{
  G4cout << "... run time evolution" << G4endl;
}


