//
// G4InteractiveSolid.cc
//
// Implementation of a messenger for constructing solids interactively
//

#include "G4InteractiveSolid.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithPargs.hh"
#include "G4UIcmdPargDouble.hh"
#include "G4UIcmdPargInteger.hh"
#include "G4UIcmdPargListDouble.hh"

#include "G4Box.hh"
#include "G4Para.hh"
#include "G4Trap.hh"
#include "G4Trd.hh"
#include "G4Sphere.hh"
#include "G4Torus.hh"
#include "G4Tubs.hh"
#include "G4Cons.hh"
#include "G4Hype.hh"
#include "G4Polycone.hh"
#include "G4Polyhedra.hh"
#include "G4EllipticalTube.hh"
#include "G4SubtractionSolid.hh"

//
// Constructor
//
G4InteractiveSolid::G4InteractiveSolid( const G4String &prefix )
{
	//
	// Hmmm... well, how about a reasonable default?
	//
	solid = new G4Box( "InteractiveBox", 1.0*m, 1.0*m, 1.0*m );
	
	//
	// Declare messenger directory
	//
	volumeDirectory = new G4UIdirectory( prefix );
	volumeDirectory->SetGuidance( "Solid construction using parameters from the command line" );


	//
	// Declare G4Box
	//
	boxArgs[0] = new G4UIcmdPargDouble( "dx", 1.0, m );
	boxArgs[1] = new G4UIcmdPargDouble( "dy", 1.0, m );
	boxArgs[2] = new G4UIcmdPargDouble( "dz", 1.0, m );
	G4String boxPath = prefix+"G4Box";
	boxCmd = new G4UIcmdWithPargs( boxPath, this, boxArgs, 3 );
	boxCmd->SetGuidance( "Declare a G4Box solid" );

	//
	// Declare G4Para
	//
	paraArgs[0] = new G4UIcmdPargDouble( "dx", 1.0, m );
	paraArgs[1] = new G4UIcmdPargDouble( "dy", 1.0, m );
	paraArgs[2] = new G4UIcmdPargDouble( "dz", 1.0, m );
	paraArgs[3] = new G4UIcmdPargDouble( "alpha", 90, deg );
	paraArgs[4] = new G4UIcmdPargDouble( "theta", 90, deg );
	paraArgs[5] = new G4UIcmdPargDouble( "phi", 90, deg );
	G4String paraPath = prefix+"G4Para";
	paraCmd = new G4UIcmdWithPargs( paraPath, this, paraArgs, 6 );
	paraCmd->SetGuidance( "Declare a G4Para solid" );

	//
	// Declare G4Trap
	//
	trapArgs[ 0] = new G4UIcmdPargDouble( "dz",    1.0, m );
	trapArgs[ 1] = new G4UIcmdPargDouble( "theta",  90, deg );
	trapArgs[ 2] = new G4UIcmdPargDouble( "phi",   1.0, m );
	trapArgs[ 3] = new G4UIcmdPargDouble( "dy1",   1.0, m );
	trapArgs[ 4] = new G4UIcmdPargDouble( "dx1",   1.0, m );
	trapArgs[ 5] = new G4UIcmdPargDouble( "dx2",   1.0, m );
	trapArgs[ 6] = new G4UIcmdPargDouble( "alpha1", 90, deg );
	trapArgs[ 7] = new G4UIcmdPargDouble( "dy2",   1.0, m );
	trapArgs[ 8] = new G4UIcmdPargDouble( "dx3",   1.0, m );
	trapArgs[ 9] = new G4UIcmdPargDouble( "dx4",   1.0, m );
	trapArgs[10] = new G4UIcmdPargDouble( "alpha2", 90, deg );
	G4String trapPath = prefix+"G4Trap";
	trapCmd = new G4UIcmdWithPargs( trapPath, this, trapArgs, 11 );
	trapCmd->SetGuidance( "Declare a G4Trap solid" );

	//
	// Declare G4Trd
	//
	trdArgs[0] = new G4UIcmdPargDouble( "dx1", 1.0, m );
	trdArgs[1] = new G4UIcmdPargDouble( "dx2", 1.0, m );
	trdArgs[2] = new G4UIcmdPargDouble( "dy1", 1.0, m );
	trdArgs[3] = new G4UIcmdPargDouble( "dy2", 1.0, m );
	trdArgs[4] = new G4UIcmdPargDouble( "dz",  1.0, m );
	G4String trdPath = prefix+"G4Trd";
	trdCmd = new G4UIcmdWithPargs( trdPath, this, trdArgs, 5 );
	trdCmd->SetGuidance( "Declare a G4Trd solid" );

	//
	// Declare G4Sphere
	//
	sphereArgs[0] = new G4UIcmdPargDouble( "rmin", 1.0, m );
	sphereArgs[1] = new G4UIcmdPargDouble( "rmax", 1.0, m );
	sphereArgs[2] = new G4UIcmdPargDouble( "startPhi", 1.0, deg );
	sphereArgs[3] = new G4UIcmdPargDouble( "deltaPhi", 1.0, deg );
	sphereArgs[4] = new G4UIcmdPargDouble( "startTheta",  1.0, deg );
	sphereArgs[5] = new G4UIcmdPargDouble( "deltaTheta",  1.0, deg );
	G4String spherePath = prefix+"G4Sphere";
	sphereCmd = new G4UIcmdWithPargs( spherePath, this, sphereArgs, 6 );
	sphereCmd->SetGuidance( "Declare a G4Sphere solid" );

	//
	// Declare G4Torus
	//
	torusArgs[0] = new G4UIcmdPargDouble( "rmin", 1.0, m );
	torusArgs[1] = new G4UIcmdPargDouble( "rmax", 1.0, m );
	torusArgs[2] = new G4UIcmdPargDouble( "rtorus", 1.0, m );
	torusArgs[3] = new G4UIcmdPargDouble( "startPhi",  1.0, deg );
	torusArgs[4] = new G4UIcmdPargDouble( "deltaPhi",  1.0, deg );
	G4String torusPath = prefix+"G4Torus";
	torusCmd = new G4UIcmdWithPargs( torusPath, this, torusArgs, 5 );
	torusCmd->SetGuidance( "Declare a G4Torus solid" );

	//
	// Declare G4Tubs
	//
	tubsArgs[0] = new G4UIcmdPargDouble( "rmin", 1.0, m );
	tubsArgs[1] = new G4UIcmdPargDouble( "rmax", 1.0, m );
	tubsArgs[2] = new G4UIcmdPargDouble( "dz",   1.0, m );
	tubsArgs[3] = new G4UIcmdPargDouble( "startPhi",  1.0, deg );
	tubsArgs[4] = new G4UIcmdPargDouble( "deltaPhi",  1.0, deg );
	G4String tubsPath = prefix+"G4Tubs";
	tubsCmd = new G4UIcmdWithPargs( tubsPath, this, tubsArgs, 5 );
	tubsCmd->SetGuidance( "Declare a G4Tubs solid" );

	//
	// Declare G4Cons
	//
	consArgs[0] = new G4UIcmdPargDouble( "rmin1", 1.0, m );
	consArgs[1] = new G4UIcmdPargDouble( "rmax1", 1.0, m );
	consArgs[2] = new G4UIcmdPargDouble( "rmin2", 1.0, m );
	consArgs[3] = new G4UIcmdPargDouble( "rmax2", 1.0, m );
	consArgs[4] = new G4UIcmdPargDouble( "dz",    1.0, m );
	consArgs[5] = new G4UIcmdPargDouble( "startPhi",  1.0, deg );
	consArgs[6] = new G4UIcmdPargDouble( "deltaPhi",  1.0, deg );
	G4String consPath = prefix+"G4Cons";
	consCmd = new G4UIcmdWithPargs( consPath, this, consArgs, 7 );
	consCmd->SetGuidance( "Declare a G4Cons solid" );

	//
	// Declare G4Hype
	//
	hypeArgs[0] = new G4UIcmdPargDouble( "innerRadius", 1.0, m );
	hypeArgs[1] = new G4UIcmdPargDouble( "outerRadius", 1.0, m );
	hypeArgs[2] = new G4UIcmdPargDouble( "innerStereo", 1.0, deg );
	hypeArgs[3] = new G4UIcmdPargDouble( "outerStereo", 1.0, deg );
	hypeArgs[4] = new G4UIcmdPargDouble( "dz",  1.0, m );
	G4String hypePath = prefix+"G4Hype";
	hypeCmd = new G4UIcmdWithPargs( hypePath, this, hypeArgs, 5 );
	hypeCmd->SetGuidance( "Declare a G4Hype solid" );
	
	//
	// Declare G4Polycone
	//
	polyconeArgs[0] = new G4UIcmdPargDouble( "phiStart", 0, deg );
	polyconeArgs[1] = new G4UIcmdPargDouble( "phiTotal", 0, deg );
	polyconeArgs[2] = new G4UIcmdPargInteger( "numRZ", 1 );
	polyconeArgs[3] = new G4UIcmdPargListDouble( "r", 100, m );
	polyconeArgs[4] = new G4UIcmdPargListDouble( "z", 100, m );
	G4String polyconePath = prefix+"G4Polycone";
	polyconeCmd = new G4UIcmdWithPargs( polyconePath, this, polyconeArgs, 5 );
	polyconeCmd->SetGuidance( "Declare a G4Polycone solid" );
	
	//
	// Declare G4Polycone2
	//
	polycone2Args[0] = new G4UIcmdPargDouble( "phiStart", 0, deg );
	polycone2Args[1] = new G4UIcmdPargDouble( "phiTotal", 0, deg );
	polycone2Args[2] = new G4UIcmdPargInteger( "numRZ", 1 );
	polycone2Args[3] = new G4UIcmdPargListDouble( "z", 100, m );
	polycone2Args[4] = new G4UIcmdPargListDouble( "rin", 100, m );
	polycone2Args[5] = new G4UIcmdPargListDouble( "rout", 100, m );
	G4String polycone2Path = prefix+"G4Polycone2";
	polycone2Cmd = new G4UIcmdWithPargs( polycone2Path, this, polycone2Args, 6 );
	polycone2Cmd->SetGuidance( "Declare a G4Polycone solid (PCON style)" );
	
	//
	// Declare G4Polyhedra
	//
	polyhedraArgs[0] = new G4UIcmdPargDouble( "phiStart", 0, deg );
	polyhedraArgs[1] = new G4UIcmdPargDouble( "phiTotal", 0, deg );
	polyhedraArgs[2] = new G4UIcmdPargInteger( "numSides", 8 );
	polyhedraArgs[3] = new G4UIcmdPargInteger( "numRZ", 1 );
	polyhedraArgs[4] = new G4UIcmdPargListDouble( "r", 100, m );
	polyhedraArgs[5] = new G4UIcmdPargListDouble( "z", 100, m );
	G4String polyhedraPath = prefix+"G4Polyhedra";
	polyhedraCmd = new G4UIcmdWithPargs( polyhedraPath, this, polyhedraArgs, 6 );
	polyhedraCmd->SetGuidance( "Declare a G4Polyhedra solid" );
	
	//
	// Declare G4Polyhedra2
	//
	polyhedra2Args[0] = new G4UIcmdPargDouble( "phiStart", 0, deg );
	polyhedra2Args[1] = new G4UIcmdPargDouble( "phiTotal", 0, deg );
	polyhedra2Args[2] = new G4UIcmdPargInteger( "numSides", 8 );
	polyhedra2Args[3] = new G4UIcmdPargInteger( "numRZ", 1 );
	polyhedra2Args[4] = new G4UIcmdPargListDouble( "z", 100, m );
	polyhedra2Args[5] = new G4UIcmdPargListDouble( "rin", 100, m );
	polyhedra2Args[6] = new G4UIcmdPargListDouble( "rout", 100, m );
	G4String polyhedra2Path = prefix+"G4Polyhedra2";
	polyhedra2Cmd = new G4UIcmdWithPargs( polyhedra2Path, this, polyhedra2Args, 7 );
	polyhedra2Cmd->SetGuidance( "Declare a G4Polyhedra solid (PGON style)" );
	
	//
	// Declare G4EllipticalTube
	//
	ellipticalTubeArgs[0] = new G4UIcmdPargDouble( "dx", 1.0, m );
	ellipticalTubeArgs[1] = new G4UIcmdPargDouble( "dy", 1.0, m );
	ellipticalTubeArgs[2] = new G4UIcmdPargDouble( "dz", 1.0, m );
	G4String ellipticalTubePath = prefix+"G4EllipticalTube";
	ellipticalTubeCmd = new G4UIcmdWithPargs( ellipticalTubePath, this, ellipticalTubeArgs, 3 );
	ellipticalTubeCmd->SetGuidance( "Declare a G4EllipticalTube solid" );
	//
	// Declare DircTest
	//
	G4String dircTestPath = prefix+"DircTest";
	dircTestCmd = new G4UIcmdWithPargs( dircTestPath, this, 0, 0 );
	dircTestCmd->SetGuidance( "Declare a DircTest solid" );
}


//
// Destructor
//
G4InteractiveSolid::~G4InteractiveSolid()
{
	if (solid) delete solid;
	
	delete boxCmd;
	DeleteArgArray( boxArgs,       sizeof(      boxArgs)/sizeof(G4UIcmdParg**) );

	delete paraCmd;
	DeleteArgArray( paraArgs,      sizeof(     paraArgs)/sizeof(G4UIcmdParg**) );

	delete trapCmd;
	DeleteArgArray( trapArgs,      sizeof(     trapArgs)/sizeof(G4UIcmdParg**) );

	delete trdCmd;
	DeleteArgArray( trdArgs,       sizeof(      trdArgs)/sizeof(G4UIcmdParg**) );

	delete sphereCmd;
	DeleteArgArray( sphereArgs,    sizeof(   sphereArgs)/sizeof(G4UIcmdParg**) );

	delete torusCmd;
	DeleteArgArray( torusArgs,     sizeof(    torusArgs)/sizeof(G4UIcmdParg**) );

	delete tubsCmd;
	DeleteArgArray( tubsArgs,      sizeof(     tubsArgs)/sizeof(G4UIcmdParg**) );

	delete consCmd;
	DeleteArgArray( consArgs,      sizeof(     consArgs)/sizeof(G4UIcmdParg**) );

	delete hypeCmd;
	DeleteArgArray( hypeArgs,      sizeof(     hypeArgs)/sizeof(G4UIcmdParg**) );

	delete polyconeCmd;
	DeleteArgArray( polyconeArgs,  sizeof( polyconeArgs)/sizeof(G4UIcmdParg**) );

	delete polyhedraCmd;
	DeleteArgArray( polyhedraArgs, sizeof(polyhedraArgs)/sizeof(G4UIcmdParg**) );

	delete ellipticalTubeCmd;
	DeleteArgArray( ellipticalTubeArgs, sizeof(ellipticalTubeArgs)/sizeof(G4UIcmdParg**) );
}


//
// DeleteArgArray
//
void G4InteractiveSolid::DeleteArgArray( G4UIcmdParg **array, const G4int nItem )
{
	G4UIcmdParg **togo = array;
	while( togo < array+nItem ) {
		delete *togo;
		togo++;
	}
}


//
// MakeMeABox
//
void G4InteractiveSolid::MakeMeABox( G4String values )
{
	if (boxCmd->GetArguments( values )) {
		delete solid;
		
		G4UIcmdPargDouble *dxArg = (G4UIcmdPargDouble *)boxArgs[0],
				  *dyArg = (G4UIcmdPargDouble *)boxArgs[1],
				  *dzArg = (G4UIcmdPargDouble *)boxArgs[2];
		
		solid = new G4Box( "interactiveBox", dxArg->GetValue(),
					   	     dyArg->GetValue(),
					  	     dzArg->GetValue() );
	}
	else
		G4cerr << "G4Box not created" << G4endl;
}


//
// MakeMeAPara
//
void G4InteractiveSolid::MakeMeAPara( G4String values )
{
	if (paraCmd->GetArguments( values )) {
		delete solid;
		
		G4UIcmdPargDouble *dxArg = (G4UIcmdPargDouble *)paraArgs[0],
				  *dyArg = (G4UIcmdPargDouble *)paraArgs[1],
				  *dzArg = (G4UIcmdPargDouble *)paraArgs[2],
				  *alphaArg = (G4UIcmdPargDouble *)paraArgs[3],
				  *thetaArg = (G4UIcmdPargDouble *)paraArgs[4],
				  *phiArg   = (G4UIcmdPargDouble *)paraArgs[5];
		
		solid = new G4Para( "interactivePara", dxArg->GetValue(),
					    	       dyArg->GetValue(),
					  	       dzArg->GetValue(),
					  	       alphaArg->GetValue(),
					  	       thetaArg->GetValue(),
					  	       phiArg->GetValue()    );
	}
	else
		G4cerr << "G4Para not created" << G4endl;
}


//
// MakeMeATrap
//
void G4InteractiveSolid::MakeMeATrap( G4String values )
{
	if (trapCmd->GetArguments( values )) {
		delete solid;
		
		G4UIcmdPargDouble **dArg = (G4UIcmdPargDouble **)trapArgs;
		
		solid = new G4Trap( "interactiveTrap", dArg[0]->GetValue(),
						       dArg[1]->GetValue(),
						       dArg[2]->GetValue(),
						       dArg[3]->GetValue(),
						       dArg[4]->GetValue(),
						       dArg[5]->GetValue(),
						       dArg[6]->GetValue(),
						       dArg[7]->GetValue(),
						       dArg[8]->GetValue(),
						       dArg[9]->GetValue(),
						       dArg[10]->GetValue() );
	}
	else
		G4cerr << "G4Trap not created" << G4endl;
}


//
// MakeMeATrd
//
void G4InteractiveSolid::MakeMeATrd( G4String values )
{
	if (trdCmd->GetArguments( values )) {
		delete solid;
		
		G4UIcmdPargDouble *dx1Arg = (G4UIcmdPargDouble *)trdArgs[0],
				  *dx2Arg = (G4UIcmdPargDouble *)trdArgs[1],
				  *dy1Arg = (G4UIcmdPargDouble *)trdArgs[2],
				  *dy2Arg = (G4UIcmdPargDouble *)trdArgs[3],
				  *dzArg  = (G4UIcmdPargDouble *)trdArgs[4];
		
		solid = new G4Trd( "interactiveTrd", dx1Arg->GetValue(),
					   	     dx2Arg->GetValue(),
					  	     dy1Arg->GetValue(),
					  	     dy2Arg->GetValue(),
					  	      dzArg->GetValue() );
	}
	else
		G4cerr << "G4Trd not created" << G4endl;
}


//
// MakeMeASphere
//
void G4InteractiveSolid::MakeMeASphere( G4String values )
{
	if (sphereCmd->GetArguments( values )) {
		delete solid;
		
		G4UIcmdPargDouble **dArg = (G4UIcmdPargDouble **)sphereArgs;
		
		solid = new G4Sphere( "interactiveSphere", dArg[0]->GetValue(),
						           dArg[1]->GetValue(),
						           dArg[2]->GetValue(),
						           dArg[3]->GetValue(),
						           dArg[4]->GetValue(),
						           dArg[5]->GetValue() );
	}
	else
		G4cerr << "G4Sphere not created" << G4endl;
}


//
// MakeMeATorus
//
void G4InteractiveSolid::MakeMeATorus( G4String values )
{
	if (torusCmd->GetArguments( values )) {
		delete solid;
		
		G4UIcmdPargDouble **dArg = (G4UIcmdPargDouble **)torusArgs;
		
		solid = new G4Torus( "interactiveTorus", dArg[0]->GetValue(),
						         dArg[1]->GetValue(),
						         dArg[2]->GetValue(),
						         dArg[3]->GetValue(),
						         dArg[4]->GetValue() );
	}
	else
		G4cerr << "G4Torus not created" << G4endl;
}


//
// MakeMeATubs
//
void G4InteractiveSolid::MakeMeATubs( G4String values )
{
	if (tubsCmd->GetArguments( values )) {
		delete solid;
		
		G4UIcmdPargDouble **dArg = (G4UIcmdPargDouble **)tubsArgs;
		
		solid = new G4Tubs( "interactiveTubs", dArg[0]->GetValue(),
						       dArg[1]->GetValue(),
						       dArg[2]->GetValue(),
						       dArg[3]->GetValue(),
						       dArg[4]->GetValue() );
	}
	else
		G4cerr << "G4Tubs not created" << G4endl;
}


//
// MakeMeACons
//
void G4InteractiveSolid::MakeMeACons( G4String values )
{
	if (consCmd->GetArguments( values )) {
		delete solid;
		
		G4UIcmdPargDouble **dArg = (G4UIcmdPargDouble **)consArgs;
		
		solid = new G4Cons( "interactiveCons", dArg[0]->GetValue(),
						       dArg[1]->GetValue(),
						       dArg[2]->GetValue(),
						       dArg[3]->GetValue(),
						       dArg[4]->GetValue(),
						       dArg[5]->GetValue(),
						       dArg[6]->GetValue() );
	}
	else
		G4cerr << "G4Cons not created" << G4endl;
}


//
// MakeMeAHype
//
void G4InteractiveSolid::MakeMeAHype( G4String values )
{
	if (hypeCmd->GetArguments( values )) {
		delete solid;
		
		G4UIcmdPargDouble **dArg = (G4UIcmdPargDouble **)hypeArgs;
		
		solid = new G4Hype( "interactiveHype", dArg[0]->GetValue(),
						       dArg[1]->GetValue(),
						       dArg[2]->GetValue(),
						       dArg[3]->GetValue(),
						       dArg[4]->GetValue()  );
	}
	else
		G4cerr << "G4Hype not created" << G4endl;
}


//
// MakeMeAPolycone
//
void G4InteractiveSolid::MakeMeAPolycone( G4String values )
{
	if (polyconeCmd->GetArguments( values )) {
		G4UIcmdPargDouble	*phiStartArg = (G4UIcmdPargDouble	*)polyconeArgs[0];
		G4UIcmdPargDouble	*phiTotalArg = (G4UIcmdPargDouble	*)polyconeArgs[1];
		G4UIcmdPargInteger	*numRZArg    = (G4UIcmdPargInteger	*)polyconeArgs[2];
		G4UIcmdPargListDouble	*rArg        = (G4UIcmdPargListDouble	*)polyconeArgs[3];
		G4UIcmdPargListDouble	*zArg        = (G4UIcmdPargListDouble	*)polyconeArgs[4];
				  
		//
		// Check consistency
		//
		G4int numRZ = numRZArg->GetValue();
		if (numRZ != rArg->GetNItem() ||
		    numRZ != zArg->GetNItem()    ) {
		    	G4cerr << "numRZ inconsistent among polycone arguments" << G4endl;
			G4cerr << "G4Polycone not created" << G4endl;
			return;
		}
		
		delete solid;
		solid = new G4Polycone( "interactivePolycone", 
					phiStartArg->GetValue(),
					phiTotalArg->GetValue(),
					numRZ,
					rArg->GetValues(),
					zArg->GetValues() );
	}
	else
		G4cerr << "G4Polycone not created" << G4endl;
}


//
// MakeMeAPolycone2
//
void G4InteractiveSolid::MakeMeAPolycone2( G4String values )
{
	if (polycone2Cmd->GetArguments( values )) {
		G4UIcmdPargDouble	*phiStartArg = (G4UIcmdPargDouble	*)polycone2Args[0];
		G4UIcmdPargDouble	*phiTotalArg = (G4UIcmdPargDouble	*)polycone2Args[1];
		G4UIcmdPargInteger	*numRZArg    = (G4UIcmdPargInteger	*)polycone2Args[2];
		G4UIcmdPargListDouble	*zArg        = (G4UIcmdPargListDouble	*)polycone2Args[3];
		G4UIcmdPargListDouble	*rInArg      = (G4UIcmdPargListDouble	*)polycone2Args[4];
		G4UIcmdPargListDouble	*rOutArg     = (G4UIcmdPargListDouble	*)polycone2Args[5];
				  
		//
		// Check consistency
		//
		G4int numRZ = numRZArg->GetValue();
		if (numRZ != zArg->GetNItem() ||
		    numRZ != rInArg->GetNItem() ||
		    numRZ != rOutArg->GetNItem()    ) {
		    	G4cerr << "numRZ inconsistent among polycone arguments" << G4endl;
			G4cerr << "G4Polycone not created" << G4endl;
			return;
		}
		
		delete solid;
		solid = new G4Polycone( "interactivePolycone", 
					phiStartArg->GetValue(),
					phiTotalArg->GetValue(),
					numRZ,
					zArg->GetValues(),
					rInArg->GetValues(),
					rOutArg->GetValues() );
	}
	else
		G4cerr << "G4Polycone not created" << G4endl;
}


//
// MakeMeAPolyhedra
//
void G4InteractiveSolid::MakeMeAPolyhedra( G4String values )
{
	if (polyhedraCmd->GetArguments( values )) {
		G4UIcmdPargDouble	*phiStartArg = (G4UIcmdPargDouble	*)polyhedraArgs[0];
		G4UIcmdPargDouble	*phiTotalArg = (G4UIcmdPargDouble	*)polyhedraArgs[1];
		G4UIcmdPargInteger	*numSidesArg = (G4UIcmdPargInteger	*)polyhedraArgs[2];
		G4UIcmdPargInteger	*numRZArg    = (G4UIcmdPargInteger	*)polyhedraArgs[3];
		G4UIcmdPargListDouble	*rArg        = (G4UIcmdPargListDouble	*)polyhedraArgs[4];
		G4UIcmdPargListDouble	*zArg        = (G4UIcmdPargListDouble	*)polyhedraArgs[5];
				  
		//
		// Check consistency
		//
		G4int numRZ = numRZArg->GetValue();
		if (numRZ != rArg->GetNItem() ||
		    numRZ != zArg->GetNItem()    ) {
		    	G4cerr << "numRZ inconsistent among polyhedra arguments" << G4endl;
			G4cerr << "G4Polyhedra not created" << G4endl;
			return;
		}
		
		delete solid;
		solid = new G4Polyhedra( "interactivePolyhedra", 
					phiStartArg->GetValue(),
					phiTotalArg->GetValue(),
					numSidesArg->GetValue(),
					numRZ,
					rArg->GetValues(),
					zArg->GetValues() );
	}
	else
		G4cerr << "G4Polyhedra not created" << G4endl;
}


//
// MakeMeAPolyhedra2
//
void G4InteractiveSolid::MakeMeAPolyhedra2( G4String values )
{
	if (polyhedra2Cmd->GetArguments( values )) {
		G4UIcmdPargDouble	*phiStartArg = (G4UIcmdPargDouble	*)polyhedra2Args[0];
		G4UIcmdPargDouble	*phiTotalArg = (G4UIcmdPargDouble	*)polyhedra2Args[1];
		G4UIcmdPargInteger	*numSidesArg = (G4UIcmdPargInteger	*)polyhedra2Args[2];
		G4UIcmdPargInteger	*numRZArg    = (G4UIcmdPargInteger	*)polyhedra2Args[3];
		G4UIcmdPargListDouble	*zArg        = (G4UIcmdPargListDouble	*)polyhedra2Args[4];
		G4UIcmdPargListDouble	*rinArg      = (G4UIcmdPargListDouble	*)polyhedra2Args[5];
		G4UIcmdPargListDouble	*routArg     = (G4UIcmdPargListDouble	*)polyhedra2Args[6];
				  
		//
		// Check consistency
		//
		G4int numRZ = numRZArg->GetValue();
		if (numRZ != zArg->GetNItem() ||
		    numRZ != rinArg->GetNItem()  ||
		    numRZ != routArg->GetNItem()    ) {
		    	G4cerr << "numRZ inconsistent among polyhedra arguments" << G4endl;
			G4cerr << "G4Polyhedra not created" << G4endl;
			return;
		}
		
		delete solid;
		solid = new G4Polyhedra( "interactivePolyhedra", 
					phiStartArg->GetValue(),
					phiTotalArg->GetValue(),
					numSidesArg->GetValue(),
					numRZ,
					zArg->GetValues(),
					rinArg->GetValues(),
					routArg->GetValues() );
	}
	else
		G4cerr << "G4Polyhedra not created" << G4endl;
}


//
// MakeMeAnEllipticalTube
//
void G4InteractiveSolid::MakeMeAnEllipticalTube( G4String values )
{
	if (ellipticalTubeCmd->GetArguments( values )) {
		delete solid;
		
		G4UIcmdPargDouble *dxArg = (G4UIcmdPargDouble *)ellipticalTubeArgs[0],
				  *dyArg = (G4UIcmdPargDouble *)ellipticalTubeArgs[1],
				  *dzArg = (G4UIcmdPargDouble *)ellipticalTubeArgs[2];
		
		solid = new G4EllipticalTube( "interactiveBox", dxArg->GetValue(),
					   	                dyArg->GetValue(),
					  	                dzArg->GetValue() );
	}
	else
		G4cerr << "G4EllipticalTube not created" << G4endl;
}


//
// MakeMeDircTest
//
void G4InteractiveSolid::MakeMeDircTest()
{
	delete solid;

	G4Tubs *outside = new G4Tubs( "OuterFrame",	// name (arbitrary) 
				      1.0*m, 		// inner radius
				      1.1*m, 		// outer radius
				      0.01*m, 		// half-thickness in z
				      -15*deg, 		// start angle
				      30*deg );		// total angle
				      
	G4Box *cutout = new G4Box( "Cutout", 	// name (arbitrary)
				   0.02*m,	// half-width (x)
				   0.25*m,	// half-height (y)
				   0.01001*m );	// half-thickness (z)
				   
	G4Transform3D tran = G4Translate3D( 1.03*m, 0.0, 0.0 );
	
	solid = new G4SubtractionSolid( "drcExample", outside, cutout, tran );
}



//
// SetNewValue
//
// Invoked by the UI when the user enters a command
//
void G4InteractiveSolid::SetNewValue( G4UIcommand *command, G4String newValues )
{
	if (command == boxCmd) 
		MakeMeABox( newValues );
	else if (command == boxCmd) 
		MakeMeABox( newValues );
	else if (command == paraCmd) 
		MakeMeAPara( newValues );
	else if (command == trapCmd) 
		MakeMeATrap( newValues );
	else if (command == trdCmd) 
		MakeMeATrd( newValues );
	else if (command == sphereCmd) 
		MakeMeASphere( newValues );
	else if (command == torusCmd) 
		MakeMeATorus( newValues );
	else if (command == tubsCmd) 
		MakeMeATubs( newValues );
	else if (command == consCmd) 
		MakeMeACons( newValues );
	else if (command == hypeCmd) 
		MakeMeAHype( newValues );
	else if (command == polyconeCmd) 
		MakeMeAPolycone( newValues );
	else if (command == polycone2Cmd) 
		MakeMeAPolycone2( newValues );
	else if (command == polyhedraCmd) 
		MakeMeAPolyhedra( newValues );
	else if (command == polyhedra2Cmd) 
		MakeMeAPolyhedra2( newValues );
	else if (command == ellipticalTubeCmd) 
		MakeMeAnEllipticalTube( newValues );
	else if (command == dircTestCmd) 
		MakeMeDircTest();
	else
		G4Exception( "Unrecognized command" );
}


//
// GetCurrentValue
//
G4String G4InteractiveSolid::GetCurrentValue( G4UIcommand *command )
{
	if (command == boxCmd) 
		return ConvertArgsToString( 	   boxArgs, sizeof(	 boxArgs)/sizeof(G4UIcmdParg**) );
	else if (command == paraCmd)
		return ConvertArgsToString( 	  paraArgs, sizeof(	paraArgs)/sizeof(G4UIcmdParg**) );
	else if (command == trapCmd)
		return ConvertArgsToString( 	  trapArgs, sizeof(	trapArgs)/sizeof(G4UIcmdParg**) );
	else if (command == trdCmd)
		return ConvertArgsToString( 	   trdArgs, sizeof(	 trdArgs)/sizeof(G4UIcmdParg**) );
	else if (command == sphereCmd)
		return ConvertArgsToString( 	sphereArgs, sizeof(   sphereArgs)/sizeof(G4UIcmdParg**) );
	else if (command == torusCmd)
		return ConvertArgsToString( 	 torusArgs, sizeof(    torusArgs)/sizeof(G4UIcmdParg**) );
	else if (command == tubsCmd)
		return ConvertArgsToString( 	  tubsArgs, sizeof(	tubsArgs)/sizeof(G4UIcmdParg**) );
	else if (command == consCmd)
		return ConvertArgsToString( 	  consArgs, sizeof(	consArgs)/sizeof(G4UIcmdParg**) );
	else if (command == hypeCmd)
		return ConvertArgsToString( 	  hypeArgs, sizeof(	hypeArgs)/sizeof(G4UIcmdParg**) );
	else if (command == polyconeCmd)
		return ConvertArgsToString(   polyconeArgs, sizeof( polyconeArgs)/sizeof(G4UIcmdParg**) );
	else if (command == polycone2Cmd)
		return ConvertArgsToString(  polycone2Args, sizeof(polycone2Args)/sizeof(G4UIcmdParg**) );
	else if (command == polyhedraCmd)
		return ConvertArgsToString(  polyhedraArgs, sizeof(polyhedraArgs)/sizeof(G4UIcmdParg**) );
	else if (command == polyhedra2Cmd)
		return ConvertArgsToString(  polyhedra2Args, sizeof(polyhedra2Args)/sizeof(G4UIcmdParg**) );
	else if (command == ellipticalTubeCmd)
		return ConvertArgsToString(  ellipticalTubeArgs, sizeof(ellipticalTubeArgs)/sizeof(G4UIcmdParg**) );
	else if (command == dircTestCmd)
		return "";
	
	G4Exception( "Unrecognized command" );
	return "foo!";
}


//
// ConvertArgsToString
//
G4String G4InteractiveSolid::ConvertArgsToString( G4UIcmdParg **array, const G4int nItem )
{
	G4String answer = "(";
	
	G4UIcmdParg **togo = array;
	while( togo < array+nItem ) {
		if (togo > array) answer += ",";
		answer += (*togo)->ConvertToString();
		togo++;
	}
	
	return answer + ")";
}


