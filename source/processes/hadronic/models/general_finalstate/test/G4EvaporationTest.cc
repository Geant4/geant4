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
// * authors in the GEANT4 collaboration.                             *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// Viktor Krylov, CERN, October 1998

//      The program is  used for testind of Evaporation Model only.
// Output of this program is sent  to file  'G4EvaporationTest.dat'
// and contains list  of events in a definite format. (be carefull,
// if the Quontity=50000, then the file size will be near 60 MByte)

#include "G4ios.hh"
#include <fstream>
#include <iomanip>
 
//#include "G4NucleiProperties.hh"
//#include "G4ExcitationHandler.hh"
#include "G4Evaporation.hh"
#include "G4CompetitiveFission.hh"
//#include "G4StatFermiBreakUp.hh"
//#include "G4DummyMF.hh"

//#include "G4DynamicParticle.hh"
//#include "G4Proton.hh"

//#include "g4templates.hh"
//#include "g4templates_hadronic_generators.hh"
 
//#include "NametoGeant3Number.cc"

void LinePrint() {
   G4cout << "-----------------------------------------------------------------" << G4endl;
}
void TitlePrint() {
    G4cout << "!" << std::setw(3) << "A " <<
              "!" << std::setw(3) << "Z " <<
              "!" << std::setw(13) << "Px     " <<
              "!" << std::setw(13) << "Py     " <<
              "!" << std::setw(13) << "Pz     " <<
              "!" << std::setw(13) << "E      " <<
	      "!" << G4endl;    
}

void FragmentPrint( const G4Fragment& rFragment ) {
    G4cout << "!" << std::setw(3) << G4int(rFragment.GetA()) <<
              "!" << std::setw(3) << G4int(rFragment.GetZ()) <<
              "!" << std::setw(13) << rFragment.GetMomentum().px() <<
              "!" << std::setw(13) << rFragment.GetMomentum().py() <<
              "!" << std::setw(13) << rFragment.GetMomentum().pz() <<
              "!" << std::setw(13) << rFragment.GetMomentum().e() <<
	      "!" << G4endl;    
}
void FragmentWrite( std::ofstream& rStream, const G4Fragment& rFragment ) {
    rStream << " " << G4int(rFragment.GetA()) <<
               " " << G4int(rFragment.GetZ()) <<
               " " << G4NucleiProperties::GetAtomicMass( rFragment.GetA(), rFragment.GetZ() ) <<
               " " << rFragment.GetExcitationEnergy() <<    

               " " << rFragment.GetMomentum().px() <<
               " " << rFragment.GetMomentum().py() <<
               " " << rFragment.GetMomentum().pz() <<
               " " << rFragment.GetMomentum().e() <<    

               " " << rFragment.GetAngularMomentum().x() <<    
               " " << rFragment.GetAngularMomentum().y() <<    
               " " << rFragment.GetAngularMomentum().z() <<    

               " " << rFragment.GetNumberOfExcitons() <<    
               " " << rFragment.GetNumberOfHoles() <<    
               " " << rFragment.GetNumberOfCharged() <<    
               " " << G4endl;    
}
int main( int argc, char* argv[] )
{
    G4cout << G4endl << "The Evaporation Test Program." << G4endl;
    G4cout << "Usage: EvaporationTest [-h] [A Z U [Quantity]]" << G4endl;
    if( argc > 5 || 
	argc > 1 && strcmp( argv[1], "-h" ) == 0 ) {
      G4cout << "Where:" << G4endl;
      G4cout << "\t-h - this message" << G4endl;
      G4cout << "\t A - atom number of a fragment A = 1..300 (defaule 238)" << G4endl;
      G4cout << "\t Z - charge of a fragment Z = 0..110 and Z <= A (default 92)" << G4endl;
      G4cout << "\t U - excitation energy U >= 0.0 * MeV (default 100.0*MeV)" << G4endl;
      G4cout << "\t Quantity - number of events N > 0 (default 1)" << G4endl;
      exit( 1 );
    }
    G4int A = argc > 1 ? atoi( argv[1] ) : 238;
    G4int Z = argc > 2 ? atoi( argv[2] ) : 92;
    G4double U = argc > 3 ? atof( argv[3] )* MeV : 100.0 * MeV;
    G4int N = argc > 4 ? atoi( argv[4] ) : 1;
    if( A == 0 || A > 300 || Z > 110 || Z > A || U < 0.0 || N <= 0 ) {
      G4cerr << "Parameters Error! ([-h] - for help)" << G4endl;
      exit( 2 );
    }

    G4cout.setf( std::ios::scientific, std::ios::floatfield );
    std::ofstream outFile( "G4EvaporationTest.dat", std::ios::out);
    outFile.setf( std::ios::scientific, std::ios::floatfield );
    /*
    // Excess Mass Test
    G4int nA;
    G4int nZ;
    for( nA = 5; nA < 300; nA++ )
      for( nZ = 3; nZ <= nA; nZ++ ) {
	G4double dT = -0.0;
	if( G4NucleiPropertiesTable::IsInTable( nZ, nA ) )
	  dT = G4NucleiPropertiesTable::GetMassExcess( nZ, nA );
	G4double dC = -0.0;
	if( nA >= nZ && nZ < 130 && nA - nZ < 200 )
	  dC = G4NucleiProperties::CameronMassExcess( nA, nZ );
	outFile << "A = " << nA << "\tZ = " << nZ << 
	  "\t Table = " << dT <<
	  "\t Cameron = " << dC << G4endl;
      }
      */
    G4Evaporation theEvaporation;
    G4double M = G4NucleiProperties::GetAtomicMass( A, Z );
    
    G4LorentzVector theLV( 0.0, 0.0, 0.0, M + U );
    G4Fragment theNucleus( A, Z, theLV );
    theNucleus.SetExcitationEnergy( U ); // its not needed for new constructor.

    G4cout << "\t"; LinePrint();
    G4cout << "Source:\t"; TitlePrint();
    G4cout << "\t"; LinePrint();

    G4cout << "\t"; FragmentPrint( theNucleus );
    G4cout << "\t"; LinePrint();

    FragmentWrite( outFile, theNucleus );

    G4cout << "Start Loop...\n";

    const NOUT = 50;
    G4int nStep = N / NOUT;

    if( N > NOUT ) {
      G4cout << "    10%  20%  30%  40%  50%  60%  70%  80%  90%  100%\n";
      G4cout << "|----|----|----|----|----|----|----|----|----|----|\n";
    }

    G4int i;
    for( i = 0; i < N; i++ ) {
      G4FragmentVector* pFragments = theEvaporation.BreakItUp( theNucleus );
      if( !pFragments )
	exit( EXIT_FAILURE );
      G4int j;
      G4int nA = 0;
      G4int nZ = 0;
      G4LorentzVector lv;

      G4int nF = pFragments->entries();
      if( N > NOUT ) {
	if( !(i % nStep) )
	  G4cout << '|' << std::flush;
      } else {
	G4cout << "\t"; LinePrint();
	G4cout << "E:" << i << '\t'; TitlePrint();
  	G4cout << "\t"; LinePrint();
      }

      outFile << i << " " << nF << G4endl;

      for( j = 0; j < nF; j++ ){
	G4Fragment* pFragment = pFragments->at(j);
	if( !pFragment )
	  exit( j );
	if( N <= NOUT ) {
	  G4cout << "F:" << j << '\t'; FragmentPrint( *pFragment );
	}
	nA += pFragment->GetA();
	nZ += pFragment->GetZ();
	lv += pFragment->GetMomentum();

	FragmentWrite( outFile, *pFragment );
      }
      if( N <= NOUT ) {
	G4cout << "\t"; LinePrint();
	G4cout << "Sum:\t"; FragmentPrint( G4Fragment( A - nA, Z - nZ, theLV - lv ) );
      }
      pFragments->clearAndDestroy();
    }
    if( N > NOUT )
      G4cout << '|' << G4endl;
    G4cout << "...End Loop\n";
	
    return EXIT_SUCCESS; 
}
 /* The End */




