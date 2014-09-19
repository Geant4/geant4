//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// G4InteractiveSolid.cc
//
// Implementation of a messenger for constructing solids interactively
//

#include "SBTrun.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include "G4InteractiveSolid.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithPargs.hh"
#include "G4UIcmdPargDouble.hh"
#include "G4UIcmdPargInteger.hh" 
#include "G4UIcmdPargListDouble.hh"

#include "G4UIcmdWithAnInteger.hh" 

#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4Para.hh"
#include "G4Sphere.hh"
#include "G4Torus.hh"
#include "G4Trap.hh"
#include "G4Trd.hh"
#include "G4GenericTrap.hh"
#include "G4Paraboloid.hh"
#include "G4Tubs.hh"
#include "G4CutTubs.hh"
#include "G4Ellipsoid.hh"
#include "G4EllipticalCone.hh"
#include "G4EllipticalTube.hh"
#include "G4ExtrudedSolid.hh"
#include "G4Hype.hh"
#include "G4Polycone.hh"
#include "G4Polyhedra.hh"
#include "G4TessellatedSolid.hh"
#include "G4TriangularFacet.hh"
#include "G4QuadrangularFacet.hh"
#include "G4Tet.hh"
#include "G4TwistedBox.hh"
#include "G4TwistedTrap.hh"
#include "G4TwistedTrd.hh"
#include "G4TwistedTubs.hh"

#include "G4IntersectionSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"

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
	// Declare G4Orb
	//
	orbArgs[0] = new G4UIcmdPargDouble( "r", 1.0, m );
	G4String orbPath = prefix+"G4Orb";
	orbCmd = new G4UIcmdWithPargs( orbPath, this, orbArgs, 1 );
	orbCmd->SetGuidance( "Declare a G4Orb solid" );

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
	// Declare G4GenericTrap
	//
	gentrapArgs[0] = new G4UIcmdPargDouble( "dz", 1.0, m );
	gentrapArgs[1] = new G4UIcmdPargListDouble("pgonX", 100, m );
        gentrapArgs[2] = new G4UIcmdPargListDouble("pgonY", 100, m );
	G4String gentrapPath = prefix+"G4GenericTrap";
	gentrapCmd = new G4UIcmdWithPargs( gentrapPath, this, gentrapArgs, 3 );
	gentrapCmd->SetGuidance( "Declare a G4GenericTrap solid" );

        //
	// Declare G4Paraboloid
	//
	parabolArgs[0] = new G4UIcmdPargDouble( "dz" , 1.0, m );
      	parabolArgs[1] = new G4UIcmdPargDouble( "dr1", 1.0, m );
	parabolArgs[2] = new G4UIcmdPargDouble( "dr2", 1.0, m );
	G4String parabolPath = prefix+"G4Paraboloid";
	parabolCmd = new G4UIcmdWithPargs( parabolPath, this, parabolArgs, 3 );
	parabolCmd->SetGuidance( "Declare a G4Paraboloid solid" );

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
	// Declare G4CutTubs
	//
	cuttubsArgs[0] = new G4UIcmdPargDouble( "rmin", 1.0, m );
	cuttubsArgs[1] = new G4UIcmdPargDouble( "rmax", 1.0, m );
	cuttubsArgs[2] = new G4UIcmdPargDouble( "dz",   1.0, m );
	cuttubsArgs[3] = new G4UIcmdPargDouble( "startPhi",  1.0, deg );
	cuttubsArgs[4] = new G4UIcmdPargDouble( "deltaPhi",  1.0, deg );
        cuttubsArgs[5] = new G4UIcmdPargListDouble( "lowNorm", 3, mm );
	cuttubsArgs[6] = new G4UIcmdPargListDouble( "highNorm", 3, mm );
	G4String cuttubsPath = prefix+"G4CutTubs";
	cuttubsCmd = new G4UIcmdWithPargs( cuttubsPath, this, cuttubsArgs, 7 );
	cuttubsCmd->SetGuidance( "Declare a G4CutTubs solid" );

	//
	// Declare G4Ellipsoid
	//
	ellipsoidArgs[0] = new G4UIcmdPargDouble( "dx", 1.0, m );
	ellipsoidArgs[1] = new G4UIcmdPargDouble( "dy", 1.0, m );
	ellipsoidArgs[2] = new G4UIcmdPargDouble( "dz", 1.0, m );
	ellipsoidArgs[3] = new G4UIcmdPargDouble( "zBottomCut", 0, m );
	ellipsoidArgs[4] = new G4UIcmdPargDouble( "zTopCut", 0, m );
	G4String ellipsoidPath = prefix+"G4Ellipsoid";
	ellipsoidCmd = new G4UIcmdWithPargs( ellipsoidPath, this, ellipsoidArgs, 5 );
	ellipsoidCmd->SetGuidance( "Declare a G4Ellipsoid solid" );

	//
	// Declare G4EllipticalCone
	//
	elConeArgs[0] = new G4UIcmdPargDouble( "dx", 1.0, 1 );
	elConeArgs[1] = new G4UIcmdPargDouble( "dy", 1.0, 1 );
	elConeArgs[2] = new G4UIcmdPargDouble( "dz", 1.0, m );
	elConeArgs[3] = new G4UIcmdPargDouble( "zTopCut", 1.0, m );
	G4String elConePath = prefix+"G4EllipticalCone";
	elConeCmd = new G4UIcmdWithPargs( elConePath, this, elConeArgs, 4 );
	elConeCmd->SetGuidance( "Declare a G4EllipticalCone solid" );

	//
	// Declare G4EllipticalTube
	//
	elTubeArgs[0] = new G4UIcmdPargDouble( "dx", 1.0, m );
	elTubeArgs[1] = new G4UIcmdPargDouble( "dy", 1.0, m );
	elTubeArgs[2] = new G4UIcmdPargDouble( "dz", 1.0, m );
	G4String elTubePath = prefix+"G4EllipticalTube";
	elTubeCmd = new G4UIcmdWithPargs( elTubePath, this, elTubeArgs, 3 );
	elTubeCmd->SetGuidance( "Declare a G4EllipticalTube solid" );

	//
	// Declare G4ExtrudedSolid
	//
	extrudedArgs[0] = new G4UIcmdPargInteger( "numPoints", 8 );
	extrudedArgs[1] = new G4UIcmdPargListDouble( "pgonx", 100, m );
	extrudedArgs[2] = new G4UIcmdPargListDouble( "pgony", 100, m );
	extrudedArgs[3] = new G4UIcmdPargInteger( "numSides", 8 );
	extrudedArgs[4] = new G4UIcmdPargListDouble(     "z", 100, m );
	extrudedArgs[5] = new G4UIcmdPargListDouble(  "offx", 100, m );
	extrudedArgs[6] = new G4UIcmdPargListDouble(  "offy", 100, m );
	extrudedArgs[7] = new G4UIcmdPargListDouble( "scale", 100, 1 );
	G4String extrudedPath = prefix+"G4ExtrudedSolid";
	extrudedCmd = new G4UIcmdWithPargs( extrudedPath, this, extrudedArgs, 8 );
	extrudedCmd->SetGuidance( "Declare a G4ExtrudedSolid solid" );
	
	//
	// Declare G4Hype
	//
	hypeArgs[0] = new G4UIcmdPargDouble( "innerRadius", 1.0, m );
	hypeArgs[1] = new G4UIcmdPargDouble( "outerRadius", 1.0, m );
	hypeArgs[2] = new G4UIcmdPargDouble( "innerStereo", 1.0, rad );
	hypeArgs[3] = new G4UIcmdPargDouble( "outerStereo", 1.0, rad );
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
	// Declare G4TessellatedSolid
	//
	tesselArgs[0] = new G4UIcmdPargInteger("num3", 8 );
	tesselArgs[1] = new G4UIcmdPargListDouble("p1in3", 100, m );
	tesselArgs[2] = new G4UIcmdPargListDouble("p2in3", 100, m );
	tesselArgs[3] = new G4UIcmdPargListDouble("p3in3", 100, m );
	tesselArgs[4] = new G4UIcmdPargInteger("num4", 8 );

	tesselArgs[5] = new G4UIcmdPargListDouble("p1in4", 100, m );
	tesselArgs[6] = new G4UIcmdPargListDouble("p2in4", 100, m );
	tesselArgs[7] = new G4UIcmdPargListDouble("p3in4", 100, m );
	tesselArgs[8] = new G4UIcmdPargListDouble("p4in4", 100, m );
	G4String tesselPath = prefix+"G4TessellatedSolid";
	tesselCmd = new G4UIcmdWithPargs( tesselPath, this, tesselArgs, 9 );
	tesselCmd->SetGuidance( "Declare a G4TessellatedSolid solid" );
	
	//
	// Declare G4TessellatedSolid2
	//
	tessel2Args[0] = new G4UIcmdPargInteger( "numPoints", 8 );
	tessel2Args[1] = new G4UIcmdPargListDouble( "pgonx", 100, m );
	tessel2Args[2] = new G4UIcmdPargListDouble( "pgony", 100, m );
	tessel2Args[3] = new G4UIcmdPargInteger( "numSides", 8 );
	tessel2Args[4] = new G4UIcmdPargListDouble(     "z", 100, m );
	tessel2Args[5] = new G4UIcmdPargListDouble(  "offx", 100, m );
	tessel2Args[6] = new G4UIcmdPargListDouble(  "offy", 100, m );
	tessel2Args[7] = new G4UIcmdPargListDouble( "scale", 100, 1 );
	G4String tessel2Path = prefix+"G4TessellatedSolid2";
	tessel2Cmd = new G4UIcmdWithPargs( tessel2Path, this, tessel2Args, 8 );
	tessel2Cmd->SetGuidance( "Declare a G4TessellatedSolid solid as an extruded solid" );
	
	//
	// Declare G4Tet
	//
	tetArgs[0] = new G4UIcmdPargListDouble( "p1", 3, m );
	tetArgs[1] = new G4UIcmdPargListDouble( "p2", 3, m );
	tetArgs[2] = new G4UIcmdPargListDouble( "p3", 3, m );
	tetArgs[3] = new G4UIcmdPargListDouble( "p4", 3, m );
	G4String tetPath = prefix+"G4Tet";
	tetCmd = new G4UIcmdWithPargs( tetPath, this, tetArgs, 4 );
	tetCmd->SetGuidance( "Declare a G4Tet solid" );
	
	//
	// Declare G4TwistedBox
	//
	twistedBoxArgs[0] = new G4UIcmdPargDouble( "phi", 0.0, deg );
	twistedBoxArgs[1] = new G4UIcmdPargDouble( "dx",  1.0, m );
	twistedBoxArgs[2] = new G4UIcmdPargDouble( "dy",  1.0, m );
	twistedBoxArgs[3] = new G4UIcmdPargDouble( "dz",  1.0, m );
	G4String twistedBoxPath = prefix+"G4TwistedBox";
	twistedBoxCmd = new G4UIcmdWithPargs( twistedBoxPath, this, twistedBoxArgs, 4 );
	twistedBoxCmd->SetGuidance( "Declare a G4TwistedBox solid" );

	//
	// Declare regular G4TwistedTrap
        // 
	//
	twistedTrapArgs[0] = new G4UIcmdPargDouble( "phi", 0.0, deg );
	twistedTrapArgs[1] = new G4UIcmdPargDouble( "dx1", 1.0, m );
	twistedTrapArgs[2] = new G4UIcmdPargDouble( "dx2", 1.0, m );
	twistedTrapArgs[3] = new G4UIcmdPargDouble( "dy",  1.0, m );
	twistedTrapArgs[4] = new G4UIcmdPargDouble( "dz",  1.0, m );
	G4String twistedTrapPath = prefix+"G4TwistedTrap";
	twistedTrapCmd = new G4UIcmdWithPargs( twistedTrapPath, this, twistedTrapArgs, 5 );
	twistedTrapCmd->SetGuidance( "Declare a regular G4TwistedTrap solid" );
        

	//
	// Declare general G4TwistedTrap
	//
	twistedTrap2Args[ 0] = new G4UIcmdPargDouble( "phi",   0.0, deg );
	twistedTrap2Args[ 1] = new G4UIcmdPargDouble( "dz",    1.0, m );
	twistedTrap2Args[ 2] = new G4UIcmdPargDouble( "theta", 0.0, deg );
	twistedTrap2Args[ 3] = new G4UIcmdPargDouble( "phi",   1.0, deg );
	twistedTrap2Args[ 4] = new G4UIcmdPargDouble( "dy1",   1.0, m );
	twistedTrap2Args[ 5] = new G4UIcmdPargDouble( "dx1",   1.0, m );
	twistedTrap2Args[ 6] = new G4UIcmdPargDouble( "dx2",   1.0, m );
	twistedTrap2Args[ 7] = new G4UIcmdPargDouble( "dy2",   1.0, m );
	twistedTrap2Args[ 8] = new G4UIcmdPargDouble( "dx3",   1.0, m );
	twistedTrap2Args[ 9] = new G4UIcmdPargDouble( "dx4",   1.0, m );
	twistedTrap2Args[10] = new G4UIcmdPargDouble( "alpha", 0.0, deg );
	G4String twistedTrap2Path = prefix+"G4TwistedTrap2";
	twistedTrap2Cmd = new G4UIcmdWithPargs( twistedTrap2Path, this, twistedTrap2Args, 11 );
	twistedTrap2Cmd->SetGuidance( "Declare a general G4TwistedTrap solid" );

	//
	// Declare G4TwistedTrd
        // 
	//
	twistedTrdArgs[0] = new G4UIcmdPargDouble( "dx1", 1.0, m );
	twistedTrdArgs[1] = new G4UIcmdPargDouble( "dx2", 1.0, m );
	twistedTrdArgs[2] = new G4UIcmdPargDouble( "dy1",  1.0, m );
	twistedTrdArgs[3] = new G4UIcmdPargDouble( "dy2",  1.0, m );
	twistedTrdArgs[4] = new G4UIcmdPargDouble( "dz",  1.0, m );
	twistedTrdArgs[5] = new G4UIcmdPargDouble( "phi", 0.0, deg );
	G4String twistedTrdPath = prefix+"G4TwistedTrd";
	twistedTrdCmd = new G4UIcmdWithPargs( twistedTrdPath, this, twistedTrdArgs, 6 );
	twistedTrdCmd->SetGuidance( "Declare a regular G4TwistedTrd solid" );
        

	//
	// Declare G4TwistedTubs
	//
	twistedTubsArgs[0] = new G4UIcmdPargDouble( "phi", 0.0, deg );
	twistedTubsArgs[1] = new G4UIcmdPargDouble( "rmin", 1.0, m );
	twistedTubsArgs[2] = new G4UIcmdPargDouble( "rmax", 1.0, m );
	twistedTubsArgs[3] = new G4UIcmdPargDouble( "zneg", 1.0, m );
	twistedTubsArgs[4] = new G4UIcmdPargDouble( "zpos", 1.0, m );
	twistedTubsArgs[5] = new G4UIcmdPargInteger( "nseg", 1 );
	twistedTubsArgs[6] = new G4UIcmdPargDouble( "totphi", 360.0, deg );
	G4String twistedTubsPath = prefix+"G4TwistedTubs";
	twistedTubsCmd = new G4UIcmdWithPargs( twistedTubsPath, this, twistedTubsArgs, 7 );
	twistedTubsCmd->SetGuidance( "Declare a G4TwistedTubs solid" );

	//
	// Declare DircTest
	//
	G4String dircTestPath = prefix+"DircTest";
	dircTestCmd = new G4UIcmdWithPargs( dircTestPath, this, 0, 0 );
	dircTestCmd->SetGuidance( "Declare a DircTest solid" );

	//
	// Declare SimpleBooleanSolid
	//
	G4String SimpleBooleanSolidPath = prefix+"SimpleBooleanSolid";
	SimpleBooleanSolidCmd = new G4UIcmdWithAnInteger( SimpleBooleanSolidPath, this); 
	boxCmd->SetGuidance( "Declare a simple Boolean solid - with a type" );
	SimpleBooleanSolidCmd->SetParameterName("Type of Boolean", false, true); // 1) Cannot omit, 2) current is default
	SimpleBooleanSolidCmd->SetDefaultValue(1);   // Used only if omitable and default is false (ie not now.)
        //     "Type of Boolean (0: int, 1: sub, 2: union): ", 1 );
	// SimpleBooleanType = new G4UIcmdPargInteger( "Type of Boolean (0: int, 1: sub, 2: union): ", 1 );

	// BooleanSolid1Cmd = new G4UIcmdWithPargs( BooleanSolid1Path, this, 0, 0 );
	// BooleanSolid1Cmd->SetGuidance( "Declare a Boolean solid #1" );
}


//
// Destructor
//
G4InteractiveSolid::~G4InteractiveSolid()
{
	if (solid) delete solid;
	
	delete boxCmd;
	DeleteArgArray( boxArgs,       sizeof(      boxArgs)/sizeof(G4UIcmdParg**) );

	delete consCmd;
	DeleteArgArray( consArgs,      sizeof(     consArgs)/sizeof(G4UIcmdParg**) );

	delete orbCmd;
	DeleteArgArray( orbArgs,       sizeof(      orbArgs)/sizeof(G4UIcmdParg**) );

	delete paraCmd;
	DeleteArgArray( paraArgs,      sizeof(     paraArgs)/sizeof(G4UIcmdParg**) );

	delete sphereCmd;
	DeleteArgArray( sphereArgs,    sizeof(   sphereArgs)/sizeof(G4UIcmdParg**) );

	delete torusCmd;
	DeleteArgArray( torusArgs,     sizeof(    torusArgs)/sizeof(G4UIcmdParg**) );

	delete trapCmd;
	DeleteArgArray( trapArgs,      sizeof(     trapArgs)/sizeof(G4UIcmdParg**) );

	delete trdCmd;
	DeleteArgArray( trdArgs,       sizeof(      trdArgs)/sizeof(G4UIcmdParg**) );

        delete gentrapCmd;
	DeleteArgArray( gentrapArgs,   sizeof(  gentrapArgs)/sizeof(G4UIcmdParg**) );

        delete parabolCmd;
	DeleteArgArray( parabolArgs,   sizeof(  parabolArgs)/sizeof(G4UIcmdParg**) );

	delete tubsCmd;
	DeleteArgArray( tubsArgs,      sizeof(     tubsArgs)/sizeof(G4UIcmdParg**) );

	delete cuttubsCmd;
	DeleteArgArray( cuttubsArgs,   sizeof(  cuttubsArgs)/sizeof(G4UIcmdParg**) );

	delete ellipsoidCmd;
	DeleteArgArray( ellipsoidArgs, sizeof(ellipsoidArgs)/sizeof(G4UIcmdParg**) );

	delete elConeCmd;
	DeleteArgArray( elConeArgs,    sizeof(   elConeArgs)/sizeof(G4UIcmdParg**) );

	delete elTubeCmd;
	DeleteArgArray( elTubeArgs,    sizeof(   elTubeArgs)/sizeof(G4UIcmdParg**) );

	delete extrudedCmd;
	DeleteArgArray( extrudedArgs,  sizeof( extrudedArgs)/sizeof(G4UIcmdParg**) );

	delete hypeCmd;
	DeleteArgArray( hypeArgs,      sizeof(     hypeArgs)/sizeof(G4UIcmdParg**) );

	delete polyconeCmd;
	DeleteArgArray( polyconeArgs,  sizeof( polyconeArgs)/sizeof(G4UIcmdParg**) );

	delete polyhedraCmd;
	DeleteArgArray( polyhedraArgs, sizeof(polyhedraArgs)/sizeof(G4UIcmdParg**) );

	delete tesselCmd;
	DeleteArgArray( tesselArgs,    sizeof(   tesselArgs)/sizeof(G4UIcmdParg**) );

	delete tessel2Cmd;
	DeleteArgArray( tessel2Args,   sizeof(  tessel2Args)/sizeof(G4UIcmdParg**) );

	delete tetCmd;
	DeleteArgArray( tetArgs,       sizeof(      tetArgs)/sizeof(G4UIcmdParg**) );

	delete twistedBoxCmd;
	DeleteArgArray( twistedBoxArgs,   sizeof(   twistedBoxArgs)/sizeof(G4UIcmdParg**) );

	delete twistedTrapCmd;
	DeleteArgArray( twistedTrapArgs,  sizeof(  twistedTrapArgs)/sizeof(G4UIcmdParg**) );

	delete twistedTrap2Cmd;
	DeleteArgArray( twistedTrap2Args, sizeof( twistedTrap2Args)/sizeof(G4UIcmdParg**) );

	delete twistedTrdCmd;
	DeleteArgArray( twistedTrdArgs,  sizeof(  twistedTrdArgs)/sizeof(G4UIcmdParg**) );

	delete twistedTubsCmd;
	DeleteArgArray( twistedTubsArgs, sizeof( twistedTubsArgs)/sizeof(G4UIcmdParg**) );
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
// MakeMeAnOrb
//
void G4InteractiveSolid::MakeMeAnOrb( G4String values )
{
	if (orbCmd->GetArguments( values )) {
		delete solid;
		
		G4UIcmdPargDouble **dArg = (G4UIcmdPargDouble **)orbArgs;
		
		solid = new G4Orb( "interactiveOrb", dArg[0]->GetValue());
	}
	else
		G4cerr << "G4Orb not created" << G4endl;
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
// MakeMeAGenericTrap
//
void G4InteractiveSolid::MakeMeAGenericTrap( G4String values )
{
	if (gentrapCmd->GetArguments( values )) {
		delete solid;
                
		G4UIcmdPargDouble *dzArg = (G4UIcmdPargDouble*)gentrapArgs[0];
                G4UIcmdPargListDouble  *pgonxArg = (G4UIcmdPargListDouble *)gentrapArgs[1],
				       *pgonyArg = (G4UIcmdPargListDouble *)gentrapArgs[2];
               		      
                              
                std::vector<G4TwoVector> polygon;
                for ( G4int i=0; i<8; ++i ) {
                  polygon.push_back(G4TwoVector(pgonxArg->GetValues()[i], pgonyArg->GetValues()[i]));
		  //G4cout<<pgonxArg->GetValues()[i]<<G4endl;
                }   
               
		                              
                solid = new G4GenericTrap("interactiveGenericTrap",dzArg->GetValue(), polygon);                                        
	}
	else
		G4cerr << "G4GenericTrap not created" << G4endl;
}
//
// MakeMeAParaboloid
//
void G4InteractiveSolid::MakeMeAParaboloid( G4String values )
{
	if (parabolCmd->GetArguments( values )) {
		delete solid;
		
		G4UIcmdPargDouble *dzArg = (G4UIcmdPargDouble *)parabolArgs[0],
				  *dr1Arg = (G4UIcmdPargDouble *)parabolArgs[1],
		                  *dr2Arg = (G4UIcmdPargDouble *)parabolArgs[2];
				  
		
		solid = new G4Paraboloid( "interactiveParaboloid", 
                                                     dzArg->GetValue(),
					   	     dr1Arg->GetValue(),
					  	     dr2Arg->GetValue() );
	}
	else
		G4cerr << "G4Paraboloid not created" << G4endl;
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
// MakeMeACutTubs
//
void G4InteractiveSolid::MakeMeACutTubs( G4String values )
{
	if (cuttubsCmd->GetArguments( values )) {
		delete solid;
		
                G4UIcmdPargDouble *drMinArg     = (G4UIcmdPargDouble *)cuttubsArgs[0],
				  *drMaxArg     = (G4UIcmdPargDouble *)cuttubsArgs[1],
				  *dzArg        = (G4UIcmdPargDouble *)cuttubsArgs[2],
                                  *dphiStartArg = (G4UIcmdPargDouble *)cuttubsArgs[3],
                                  *dphiDeltaArg = (G4UIcmdPargDouble *)cuttubsArgs[4];
        	G4UIcmdPargListDouble    *lArg  = (G4UIcmdPargListDouble*)cuttubsArgs[5];
		G4UIcmdPargListDouble    *hArg  = (G4UIcmdPargListDouble*)cuttubsArgs[6];
				  

                G4ThreeVector lowNorm  = G4ThreeVector(lArg->GetValues()[0],
                                               lArg->GetValues()[1],lArg->GetValues()[2]);
                G4ThreeVector highNorm = G4ThreeVector(hArg->GetValues()[0],
                                               hArg->GetValues()[1],hArg->GetValues()[2]);

		
		solid = new G4CutTubs( "interactiveCutTubs", drMinArg->GetValue(),
						       drMaxArg->GetValue(),
						       dzArg->GetValue(),
						       dphiStartArg->GetValue(),
				                       dphiDeltaArg->GetValue(),
                                                       lowNorm,highNorm );
	}
	else
		G4cerr << "G4CutTubs not created" << G4endl;
}


//
// MakeMeAnEllipsoid
//
void G4InteractiveSolid::MakeMeAnEllipsoid( G4String values )
{
	if (ellipsoidCmd->GetArguments( values )) {
		delete solid;
		
		G4UIcmdPargDouble *dxArg = (G4UIcmdPargDouble *)ellipsoidArgs[0],
				  *dyArg = (G4UIcmdPargDouble *)ellipsoidArgs[1],
				  *dzArg = (G4UIcmdPargDouble *)ellipsoidArgs[2],
                                  *pzBottomCutArg = (G4UIcmdPargDouble *)ellipsoidArgs[3],
                                  *pzTopCutArg = (G4UIcmdPargDouble *)ellipsoidArgs[4];
		
		solid = new G4Ellipsoid( "interactiveEllipsoid", 
                                         dxArg->GetValue(),
					 dyArg->GetValue(),
					 dzArg->GetValue(),
                                         pzBottomCutArg->GetValue(),
                                         pzTopCutArg->GetValue() );
	}
	else
		G4cerr << "G4Ellipsoid not created" << G4endl;
}


//
// MakeMeAnEllipticalTube
//
void G4InteractiveSolid::MakeMeAnEllipticalCone( G4String values )
{
	if (elConeCmd->GetArguments( values )) {
		delete solid;
		
		G4UIcmdPargDouble *dxArg = (G4UIcmdPargDouble *)elConeArgs[0],
				  *dyArg = (G4UIcmdPargDouble *)elConeArgs[1],
				  *dzArg = (G4UIcmdPargDouble *)elConeArgs[2],
                                  *pzTopCutArg = (G4UIcmdPargDouble *)elConeArgs[3];
                                  
                G4cout << "Making G4EllipticalCone: " 
                       <<  dxArg->GetValue() << " "
                       <<  dyArg->GetValue() << " "
                       <<  dzArg->GetValue() << " "
                       <<  pzTopCutArg->GetValue() << G4endl;                            
		
		solid = new G4EllipticalCone( "interactiveEllipticalCone", 
                                              dxArg->GetValue(),
					      dyArg->GetValue(),
					      dzArg->GetValue(),
                                              pzTopCutArg->GetValue() );
	}
	else
		G4cerr << "G4EllipticalCone not created" << G4endl;
}

//
// MakeMeAnEllipticalTube
//
void G4InteractiveSolid::MakeMeAnEllipticalTube( G4String values )
{
	if (elTubeCmd->GetArguments( values )) {
		delete solid;
		
		G4UIcmdPargDouble *dxArg = (G4UIcmdPargDouble *)elTubeArgs[0],
				  *dyArg = (G4UIcmdPargDouble *)elTubeArgs[1],
				  *dzArg = (G4UIcmdPargDouble *)elTubeArgs[2];
		
		solid = new G4EllipticalTube( "interactiveEllipticalTube", 
                                              dxArg->GetValue(),
					      dyArg->GetValue(),
					      dzArg->GetValue() );
	}
	else
		G4cerr << "G4EllipticalTube not created" << G4endl;
}

//
// MakeMeAnExtrudedSolid
//
void G4InteractiveSolid::MakeMeAnExtrudedSolid( G4String values )
{
	if (extrudedCmd->GetArguments( values )) {
		delete solid;
                
		G4UIcmdPargInteger *numPointsArg = (G4UIcmdPargInteger*)extrudedArgs[0];
                G4UIcmdPargListDouble  *pgonxArg = (G4UIcmdPargListDouble *)extrudedArgs[1],
				       *pgonyArg = (G4UIcmdPargListDouble *)extrudedArgs[2];
                                  
                                  
		G4UIcmdPargInteger  *numSidesArg = (G4UIcmdPargInteger*)extrudedArgs[3];
                G4UIcmdPargListDouble      *zArg = (G4UIcmdPargListDouble *)extrudedArgs[4],
				        *offxArg = (G4UIcmdPargListDouble *)extrudedArgs[5],
				        *offyArg = (G4UIcmdPargListDouble *)extrudedArgs[6],
				       *scaleArg = (G4UIcmdPargListDouble *)extrudedArgs[7];
                
                std::vector<G4TwoVector> polygon;
                for ( G4int i=0; i<numPointsArg->GetValue(); ++i ) {
                  polygon.push_back(G4TwoVector(pgonxArg->GetValues()[i], pgonyArg->GetValues()[i]));
                  G4cout<<"extr="<<polygon[i]<<G4endl;
                }   
		
                std::vector<G4ExtrudedSolid::ZSection> zsections;
                for ( G4int i=0; i<numSidesArg->GetValue(); ++i ) {
                  zsections.push_back(G4ExtrudedSolid::ZSection(
                                        zArg->GetValues()[i],
                                        G4TwoVector(offxArg->GetValues()[i], offyArg->GetValues()[i]),
                                        scaleArg->GetValues()[i]));
                }
                
                solid = new G4ExtrudedSolid("interactiveExtrudedSolid", polygon, zsections);                                        
	}
	else
		G4cerr << "G4ExtrudedSolid not created" << G4endl;
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
                G4cout << *solid << G4endl;                                                       
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
// MakeMeATessellatedSolid
//
void G4InteractiveSolid::MakeMeATessellatedSolid( G4String values )
{
	if (tesselCmd->GetArguments( values )) {

		G4UIcmdPargInteger *num3Arg = (G4UIcmdPargInteger*)tesselArgs[0];
		G4UIcmdPargListDouble *p1in3Arg = (G4UIcmdPargListDouble   *)tesselArgs[1],
		                      *p2in3Arg = (G4UIcmdPargListDouble   *)tesselArgs[2],
	 	                      *p3in3Arg = (G4UIcmdPargListDouble   *)tesselArgs[3];
                
		G4UIcmdPargInteger *num4Arg = (G4UIcmdPargInteger*)tesselArgs[4];
		G4UIcmdPargListDouble *p1in4Arg = (G4UIcmdPargListDouble   *)tesselArgs[5],
		                      *p2in4Arg = (G4UIcmdPargListDouble   *)tesselArgs[6],
		                      *p3in4Arg = (G4UIcmdPargListDouble   *)tesselArgs[7],
		                      *p4in4Arg = (G4UIcmdPargListDouble   *)tesselArgs[8];
				  
		//
		// Check consistency
		//
		G4int num3 = num3Arg->GetValue();
		G4int nump1in3 = p1in3Arg->GetNItem();
		G4int nump2in3 = p2in3Arg->GetNItem();
		G4int nump3in3 = p3in3Arg->GetNItem();
		if ( nump1in3 != 3*num3 || nump2in3 != 3*num3 || nump3in3 != 3*num3 ) {
		    	G4cerr << "Wrong number of points coordinates among triangular tessel arguments" << G4endl;
			G4cerr << "G4TessellatedSolid not created" << G4endl;
			return;
		}
		G4int num4 = num4Arg->GetValue();
		G4int nump1in4 = p1in4Arg->GetNItem();
		G4int nump2in4 = p2in4Arg->GetNItem();
		G4int nump3in4 = p3in4Arg->GetNItem();
		G4int nump4in4 = p4in4Arg->GetNItem();
		if ( nump1in4 != 3*num4 || nump2in4 != 3*num4 || nump3in4 != 3*num4 || nump4in4 != 3*num4) {
		    	G4cerr << "Wrong number of points coordinates among quadrangular tessel arguments" << G4endl;
			G4cerr << "G4TessellatedSolid not created" << G4endl;
			return;
		}
		
		delete solid;
                                   
		G4TessellatedSolid* tessel 
                  = new G4TessellatedSolid( "interactiveTessellatedSolid");
                  
                for ( G4int i=0; i<num3; ++i) {
                  G4ThreeVector p1(p1in3Arg->GetValues()[3*i+0], p1in3Arg->GetValues()[3*i+1], p1in3Arg->GetValues()[3*i+2]);
                  G4ThreeVector p2(p2in3Arg->GetValues()[3*i+0], p2in3Arg->GetValues()[3*i+1], p2in3Arg->GetValues()[3*i+2]);
                  G4ThreeVector p3(p3in3Arg->GetValues()[3*i+0], p3in3Arg->GetValues()[3*i+1], p3in3Arg->GetValues()[3*i+2]);
                  tessel->AddFacet(new G4TriangularFacet(p1, p2, p3, ABSOLUTE));
                }  
                
                for ( G4int i=0; i<num4; ++i) {
                  G4ThreeVector p1(p1in4Arg->GetValues()[3*i+0], p1in4Arg->GetValues()[3*i+1], p1in4Arg->GetValues()[3*i+2]);
                  G4ThreeVector p2(p2in4Arg->GetValues()[3*i+0], p2in4Arg->GetValues()[3*i+1], p2in4Arg->GetValues()[3*i+2]);
                  G4ThreeVector p3(p3in4Arg->GetValues()[3*i+0], p3in4Arg->GetValues()[3*i+1], p3in4Arg->GetValues()[3*i+2]);
                  G4ThreeVector p4(p4in4Arg->GetValues()[3*i+0], p4in4Arg->GetValues()[3*i+1], p4in4Arg->GetValues()[3*i+2]);
                  tessel->AddFacet(new G4QuadrangularFacet(p1, p2, p3, p4, ABSOLUTE));
                }
                tessel->SetSolidClosed(true);
                G4cout << *tessel << G4endl;
                solid = tessel;  
               
	}
	else
		G4cerr << "G4TessellatedSolid not created" << G4endl;
}

//
// MakeMeAnTessellatedSolid2
//
void G4InteractiveSolid::MakeMeATessellatedSolid2( G4String values )
{
	if (tessel2Cmd->GetArguments( values )) {
		delete solid;
                
		G4UIcmdPargInteger *numPointsArg = (G4UIcmdPargInteger*)tessel2Args[0];
                G4UIcmdPargListDouble  *pgonxArg = (G4UIcmdPargListDouble *)tessel2Args[1],
				       *pgonyArg = (G4UIcmdPargListDouble *)tessel2Args[2];
                                  
                                  
		G4UIcmdPargInteger  *numSidesArg = (G4UIcmdPargInteger*)tessel2Args[3];
                G4UIcmdPargListDouble      *zArg = (G4UIcmdPargListDouble *)tessel2Args[4],
				        *offxArg = (G4UIcmdPargListDouble *)tessel2Args[5],
				        *offyArg = (G4UIcmdPargListDouble *)tessel2Args[6],
				       *scaleArg = (G4UIcmdPargListDouble *)tessel2Args[7];
                
                std::vector<G4TwoVector> polygon;
                for ( G4int i=0; i<numPointsArg->GetValue(); ++i ) {
                  polygon.push_back(G4TwoVector(pgonxArg->GetValues()[i], pgonyArg->GetValues()[i]));
                }   
		
                std::vector<G4ExtrudedSolid::ZSection> zsections;
                for ( G4int i=0; i<numSidesArg->GetValue(); ++i ) {
                  zsections.push_back(G4ExtrudedSolid::ZSection(
                                        zArg->GetValues()[i],
                                        G4TwoVector(offxArg->GetValues()[i], offyArg->GetValues()[i]),
                                        scaleArg->GetValues()[i]));
                }
                
                G4ExtrudedSolid* xtru
                  = new G4ExtrudedSolid("interactiveTessellatedSolid", polygon, zsections); 
//                solid = new G4TessellatedSolid(*xtru);
//                delete xtru;
                solid = xtru;
                                                       
	}
	else
		G4cerr << "G4TessellatedSolid not created" << G4endl;
}


//
// MakeMeATet
//
void G4InteractiveSolid::MakeMeATet( G4String values )
{
	if (tetCmd->GetArguments( values )) {
		G4UIcmdPargListDouble	*p1Arg = (G4UIcmdPargListDouble   *)tetArgs[0];
		G4UIcmdPargListDouble	*p2Arg = (G4UIcmdPargListDouble   *)tetArgs[1];
		G4UIcmdPargListDouble	*p3Arg = (G4UIcmdPargListDouble   *)tetArgs[2];
		G4UIcmdPargListDouble	*p4Arg = (G4UIcmdPargListDouble   *)tetArgs[3];
				  
		//
		// Check consistency
		//
		G4int numCoor1 = p1Arg->GetNItem();
		G4int numCoor2 = p2Arg->GetNItem();
		G4int numCoor3 = p3Arg->GetNItem();
		G4int numCoor4 = p4Arg->GetNItem();
		if (numCoor1 != 3 || numCoor2 != 3 || numCoor3 != 3 || numCoor4 != 3 ) {
		    	G4cerr << "Wrong number of points coordinates among tet arguments" << G4endl;
			G4cerr << "G4Tet not created" << G4endl;
			return;
		}
		
		delete solid;
                G4ThreeVector p1(p1Arg->GetValues()[0], p1Arg->GetValues()[1], p1Arg->GetValues()[2]);
                G4ThreeVector p2(p2Arg->GetValues()[0], p2Arg->GetValues()[1], p2Arg->GetValues()[2]);
                G4ThreeVector p3(p3Arg->GetValues()[0], p3Arg->GetValues()[1], p3Arg->GetValues()[2]);
                G4ThreeVector p4(p4Arg->GetValues()[0], p4Arg->GetValues()[1], p4Arg->GetValues()[2]);
                                   
		solid = new G4Tet( "interactiveTet", p1, p2, p3, p4 );
	}
	else
		G4cerr << "G4Tet not created" << G4endl;
}


//
// MakeMeABox
//
void G4InteractiveSolid::MakeMeATwistedBox( G4String values )
{
	if (twistedBoxCmd->GetArguments( values )) {
		delete solid;
		
		G4UIcmdPargDouble *phiArg = (G4UIcmdPargDouble *)twistedBoxArgs[0],
                                  *dxArg = (G4UIcmdPargDouble *)twistedBoxArgs[1],
				  *dyArg = (G4UIcmdPargDouble *)twistedBoxArgs[2],
				  *dzArg = (G4UIcmdPargDouble *)twistedBoxArgs[3];
		
		solid = new G4TwistedBox( "interactiveTwistedBox", 
                                          phiArg->GetValue(),
                                          dxArg->GetValue(),
					  dyArg->GetValue(),
					  dzArg->GetValue() );
                                          
                G4cout << *solid << G4endl;                                          
	}
	else
		G4cerr << "G4TwistedBox not created" << G4endl;
}

//
// MakeMeATwistedTrap
//
void G4InteractiveSolid::MakeMeATwistedTrap( G4String values )
{
	if (twistedTrapCmd->GetArguments( values )) {
		delete solid;
		
		G4UIcmdPargDouble *phiArg = (G4UIcmdPargDouble *)twistedTrapArgs[0],
				  *dx1Arg = (G4UIcmdPargDouble *)twistedTrapArgs[1],
				  *dx2Arg = (G4UIcmdPargDouble *)twistedTrapArgs[2],
				  *dyArg  = (G4UIcmdPargDouble *)twistedTrapArgs[3],
				  *dzArg  = (G4UIcmdPargDouble *)twistedTrapArgs[4];
		
		solid = new G4TwistedTrap( "interactiveTwistedTrap", 
                                           phiArg->GetValue(),
                                           dx1Arg->GetValue(),
					   dx2Arg->GetValue(),
					   dyArg->GetValue(),
					   dzArg->GetValue() );
	}
	else
		G4cerr << "G4TwistedTrap not created" << G4endl;
}


//
// MakeMeATwistedTrap
//
void G4InteractiveSolid::MakeMeATwistedTrap2( G4String values )
{
	if (twistedTrap2Cmd->GetArguments( values )) {
		delete solid;
		
		G4UIcmdPargDouble **dArg = (G4UIcmdPargDouble **)twistedTrap2Args;
		
		solid = new G4TwistedTrap( "interactiveTwistedTrap2", 
                                           dArg[0]->GetValue(),
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
		G4cerr << "G4TwistedTrap2 not created" << G4endl;
}

//
// MakeMeATwistedTrd
//
void G4InteractiveSolid::MakeMeATwistedTrd( G4String values )
{
	if (twistedTrdCmd->GetArguments( values )) {
		delete solid;
		
		G4UIcmdPargDouble *dx1Arg = (G4UIcmdPargDouble *)twistedTrdArgs[0],
				  *dx2Arg = (G4UIcmdPargDouble *)twistedTrdArgs[1],
				  *dy1Arg = (G4UIcmdPargDouble *)twistedTrdArgs[2],
				  *dy2Arg = (G4UIcmdPargDouble *)twistedTrdArgs[3],
				  *dzArg  = (G4UIcmdPargDouble *)twistedTrdArgs[4],
		                  *phiArg = (G4UIcmdPargDouble *)twistedTrdArgs[5];
				  
		solid = new G4TwistedTrd( "interactiveTwistedTrd", 
                                           dx1Arg->GetValue(),
					   dx2Arg->GetValue(),
					   dy1Arg->GetValue(),
					   dy2Arg->GetValue(),
					   dzArg->GetValue(),
                                           phiArg->GetValue() );
	} 
	else
		G4cerr << "G4TwistedTrd not created" << G4endl;
}


//
// MakeMeATwistedTubs
//
void G4InteractiveSolid::MakeMeATwistedTubs( G4String values )
{
	if (twistedTubsCmd->GetArguments( values )) {
		delete solid;
		
		G4UIcmdPargDouble  *phiArg    = (G4UIcmdPargDouble  *)twistedTubsArgs[0];
		G4UIcmdPargDouble  *rminArg   = (G4UIcmdPargDouble  *)twistedTubsArgs[1];
		G4UIcmdPargDouble  *rmaxArg   = (G4UIcmdPargDouble  *)twistedTubsArgs[2];
		G4UIcmdPargDouble  *znegArg   = (G4UIcmdPargDouble  *)twistedTubsArgs[3];
		G4UIcmdPargDouble  *zposArg   = (G4UIcmdPargDouble  *)twistedTubsArgs[4];
		G4UIcmdPargInteger *nsegArg   = (G4UIcmdPargInteger *)twistedTubsArgs[5];
		G4UIcmdPargDouble  *totphiArg = (G4UIcmdPargDouble  *)twistedTubsArgs[6];
		
		solid = new G4TwistedTubs( "interactiveTwistedTubs", 
                                           phiArg->GetValue(),
					   rminArg->GetValue(),
					   rmaxArg->GetValue(),
					   znegArg->GetValue(),
					   zposArg->GetValue(),
                                           nsegArg->GetValue(),
                                           totphiArg->GetValue() );
                                           
                G4cout << *solid << G4endl;                                           
	}
	else
		G4cerr << "G4TwistedTubs not created" << G4endl;
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
// BooleanSolid1Test
//
void G4InteractiveSolid::MakeMeASimpleBooleanSolid(G4String values)
{
  /*
    G4IntersectionSolid.hh  G4SubtractionSolid.hh  G4UnionSolid.hh
    all CSGs : Box Tubs Sphere Cons Torus
    So: Boolean type operation and 2 CSG Objects with parameters for each (..)
    plus a transformation to apply to the second solid
   */
        delete solid;

	BooleanOp OperationType;
	// OperationType = INTERSECTION;
	OperationType = SUBTRACTION;

	// if ( SimpleBooleanSolidCmd->GetArguments( values )) {
        G4int intValue= SimpleBooleanSolidCmd->GetNewIntValue(values); 

	// if (     intValue <= std::max( SUBTRACTION, std::max( UNION, INTERSECTION ) ) 
        //  &&  intValue >= std::min( SUBTRACTION, std::min( UNION, INTERSECTION ) ) )

        if( (intValue == SUBTRACTION) || (intValue == INTERSECTION) || (intValue == UNION) )
        {
           BooleanOp booleanOpType;

           booleanOpType = (BooleanOp) intValue; 

           //  if( intValue  == SUBTRACTION ) { booleanOpType = SUBTRACTION; }
           // else if ( intValue == INTERSECTION ) { booleanOpType = INTERSECTION; }
           // else if ( intValue == INTERSECTION ) { booleanOpType = UNION; }

           if( ( booleanOpType == SUBTRACTION ) || ( booleanOpType == INTERSECTION ) || ( booleanOpType == UNION ) )
           {
              OperationType= booleanOpType;
           }else{
              G4cerr << " ERROR> Error in converting value for boolean operation type.  The value " 
                     << booleanOpType << " is not valid. " << G4endl;
           }
        }
        else
        {
           G4cerr << " ERROR> Requested value for boolean operation type is out of range. " << G4endl;
           G4cerr << "  Value requested = " << intValue << G4endl;
           G4cerr << "  Allowed values are: " << G4endl;
           G4cerr << "    Intersection: " << INTERSECTION << G4endl
                  << "    Subtraction:  " << SUBTRACTION << G4endl
                  << "    Union:        " << UNION <<  G4endl;
           G4cerr << "  Using default value: " << OperationType << G4endl;
        }

	/*
	G4Tubs *outside = new G4Tubs( "OuterFrame",	// name (arbitrary) 
				      1.0*m, 		// inner radius
				      1.1*m, 		// outer radius
				      0.50*m, 		// half-thickness in z
				      0*deg, 		// start angle
				      180*deg );		// total angle
	*/
	/*
	G4Torus *outside = new G4Torus( "interactiveTorus",
					0.2*m,
				        0.8*m,
				        1.4*m,
				        0*deg,
				        360*deg );
	*/
	
	G4Cons *outside = new G4Cons( "OuterFrame",
				      0.6*m, // pRmin1
				      1.0*m, // pRmax1
				      0.2*m, // pRmin2
				      0.8*m, // pRmax2
				      0.2*m,
				      0*deg,
				      180*deg );
		
	/* Dirctest Box cutout
	G4Box *cutout = new G4Box( "Cutout", 	// name (arbitrary)
				   0.02*m,	// half-width (x)
				   0.25*m,	// half-height (y)
				   0.01001*m );	// half-thickness (z)
	*/
	
	/*
	G4Tubs *cutout = new G4Tubs("AnotherTubs",
				    1.0*m,
				    1.1*m,
				    0.50*m,
				    0*deg,
				    180*deg
				    );
	*/

	
	 G4Cons *cutout = new G4Cons( "OuterFrame",
				      0.6*m, // pRmin1
				      1.0*m, // pRmax1
				      0.2*m, // pRmin2
				      0.8*m, // pRmax2
				      0.2*m,
				      0*deg,
				      180*deg );
	
	/*
	G4Torus *cutout = new G4Torus( "interactiveTorus",
					0.2*m,
				        0.8*m,
				        1.4*m,
				        0*deg,
				        360*deg );

	*/
	
	
	G4RotationMatrix rm;
	rm.rotateY(pi/4.0);

	G4Transform3D tran = G4Transform3D(rm,G4ThreeVector(0.0,0.0,0.0));
	
	/* G4Transform3D tran = G4Translate3D( 1.03*m, 0.0, 0.0 ); */

	switch (OperationType) {
	case INTERSECTION:
	  solid = new G4IntersectionSolid( "drcExample", outside, cutout, tran );
	  break;	
	case SUBTRACTION:
	  solid = new G4SubtractionSolid( "drcExample", outside, cutout, tran );
	  break;	
	case UNION:
	  solid = new G4UnionSolid( "drcExample", outside, cutout, tran );
	  break;	
	}
	
}


// G4VPhysicalVolume* ExN03DetectorConstruction::MakeMeDircTest()
  
//
// SetNewValue
//
// Invoked by the UI when the user enters a command
//
void G4InteractiveSolid::SetNewValue( G4UIcommand *command, G4String newValues )
{
  /*
    We want to retrieve the current solid
    So we keep the current solid command
   */


  G4String CurrentSolid = command->GetCommandPath() + " " + newValues ;

  SBTrun::SetCurrentSolid (CurrentSolid);

	if (command == boxCmd) 
		MakeMeABox( newValues );
	else if (command == consCmd) 
		MakeMeACons( newValues );
	else if (command == orbCmd) 
		MakeMeAnOrb( newValues );
	else if (command == paraCmd) 
		MakeMeAPara( newValues );
	else if (command == sphereCmd) 
		MakeMeASphere( newValues );
	else if (command == torusCmd) 
		MakeMeATorus( newValues );
	else if (command == trapCmd) 
		MakeMeATrap( newValues );
	else if (command == trdCmd) 
		MakeMeATrd( newValues );
        else if (command == gentrapCmd) 
		MakeMeAGenericTrap( newValues );
	else if (command == parabolCmd) 
		MakeMeAParaboloid( newValues );
	else if (command == tubsCmd) 
		MakeMeATubs( newValues );
	else if (command == cuttubsCmd) 
		MakeMeACutTubs( newValues );
	else if (command == ellipsoidCmd) 
		MakeMeAnEllipsoid( newValues );
	else if (command == elConeCmd) 
		MakeMeAnEllipticalCone( newValues );
	else if (command == elTubeCmd) 
		MakeMeAnEllipticalTube( newValues );
	else if (command == extrudedCmd) 
		MakeMeAnExtrudedSolid( newValues );
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
	else if (command == tesselCmd) 
		MakeMeATessellatedSolid( newValues );
	else if (command == tessel2Cmd) 
		MakeMeATessellatedSolid2( newValues );
	else if (command == tetCmd) 
		MakeMeATet( newValues );
	else if (command == twistedBoxCmd) 
		MakeMeATwistedBox( newValues );
	else if (command == twistedTrapCmd) 
		MakeMeATwistedTrap( newValues );
	else if (command == twistedTrap2Cmd) 
		MakeMeATwistedTrap2( newValues );
	else if (command == twistedTrdCmd) 
		MakeMeATwistedTrd( newValues );
	else if (command == dircTestCmd) 
		MakeMeDircTest();
	else if (command == twistedTubsCmd) 
		MakeMeATwistedTubs( newValues );
	else if (command == SimpleBooleanSolidCmd) 
		MakeMeASimpleBooleanSolid(newValues);
	// else if (command == TestBooleanSolidCmd) a
	//	MakeMeTestBooleanSolid(newValues);
	// else if (command == BooleanSolid1Cmd) 
	//	StoreConstituentSolid(newValues);

	//else if (command == BooleanSolid1Cmd) 
	//	MakeMeGeneralBooleanSolid(newValues);

	/* Here to add new boolean solids */

	else
	  G4Exception("G4InteractiveSolid","SBT007", FatalErrorInArgument,
		      "Unrecognized command");
}


//
// GetCurrentValue
//
G4String G4InteractiveSolid::GetCurrentValue( G4UIcommand *command )
{
	if (command == boxCmd) 
		return ConvertArgsToString( 	   boxArgs, sizeof(	 boxArgs)/sizeof(G4UIcmdParg**) );
	else if (command == consCmd)
		return ConvertArgsToString( 	  consArgs, sizeof(	consArgs)/sizeof(G4UIcmdParg**) );
	else if (command == orbCmd)
		return ConvertArgsToString( 	   orbArgs, sizeof(      orbArgs)/sizeof(G4UIcmdParg**) );
	else if (command == paraCmd)
		return ConvertArgsToString( 	  paraArgs, sizeof(	paraArgs)/sizeof(G4UIcmdParg**) );
	else if (command == sphereCmd)
		return ConvertArgsToString( 	sphereArgs, sizeof(   sphereArgs)/sizeof(G4UIcmdParg**) );
	else if (command == torusCmd)
		return ConvertArgsToString( 	 torusArgs, sizeof(    torusArgs)/sizeof(G4UIcmdParg**) );
	else if (command == trapCmd)
		return ConvertArgsToString( 	  trapArgs, sizeof(	trapArgs)/sizeof(G4UIcmdParg**) );
	else if (command == trdCmd)
		return ConvertArgsToString( 	   trdArgs, sizeof(	 trdArgs)/sizeof(G4UIcmdParg**) );
        else if (command == gentrapCmd)
		return ConvertArgsToString( 	gentrapArgs, sizeof(  gentrapArgs)/sizeof(G4UIcmdParg**) );
        else if (command == parabolCmd)
		return ConvertArgsToString(    parabolArgs, sizeof(  parabolArgs)/sizeof(G4UIcmdParg**) );
	else if (command == tubsCmd)
		return ConvertArgsToString( 	  tubsArgs, sizeof(	tubsArgs)/sizeof(G4UIcmdParg**) );
        else if (command == cuttubsCmd)
		return ConvertArgsToString(    cuttubsArgs, sizeof(   cuttubsArgs)/sizeof(G4UIcmdParg**) );
	else if (command == ellipsoidCmd)
		return ConvertArgsToString(  ellipsoidArgs, sizeof( ellipsoidArgs)/sizeof(G4UIcmdParg**) );
	else if (command == elConeCmd)
		return ConvertArgsToString(     elConeArgs, sizeof(    elConeArgs)/sizeof(G4UIcmdParg**) );
	else if (command == elTubeCmd)
		return ConvertArgsToString(     elTubeArgs, sizeof(    elTubeArgs)/sizeof(G4UIcmdParg**) );
	else if (command == extrudedCmd)
		return ConvertArgsToString(   extrudedArgs, sizeof(  extrudedArgs)/sizeof(G4UIcmdParg**) );
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
	else if (command == tesselCmd)
		return ConvertArgsToString(     tesselArgs, sizeof(    tesselArgs)/sizeof(G4UIcmdParg**) );
	else if (command == tessel2Cmd)
		return ConvertArgsToString(    tessel2Args, sizeof(   tessel2Args)/sizeof(G4UIcmdParg**) );
	else if (command == tetCmd)
		return ConvertArgsToString(        tetArgs, sizeof(       tetArgs)/sizeof(G4UIcmdParg**) );
	else if (command == twistedBoxCmd) 
		return ConvertArgsToString( twistedBoxArgs, sizeof(  twistedBoxArgs)/sizeof(G4UIcmdParg**) );
	else if (command == twistedTrapCmd)
		return ConvertArgsToString(twistedTrapArgs, sizeof( twistedTrapArgs)/sizeof(G4UIcmdParg**) );
	else if (command == twistedTrap2Cmd)
		return ConvertArgsToString(twistedTrap2Args,sizeof(twistedTrap2Args)/sizeof(G4UIcmdParg**) );
	else if (command == twistedTrdCmd)
		return ConvertArgsToString( twistedTrdArgs, sizeof(  twistedTrdArgs)/sizeof(G4UIcmdParg**) );
	else if (command == twistedTubsCmd)
		return ConvertArgsToString(twistedTubsArgs, sizeof( twistedTubsArgs)/sizeof(G4UIcmdParg**) );
	else if (command == dircTestCmd)
		return "";
	
	G4Exception("G4InteractiveSolid","SBT008", FatalErrorInArgument,
		    "Unrecognized command");
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


