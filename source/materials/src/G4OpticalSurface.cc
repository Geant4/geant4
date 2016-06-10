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
// $Id: G4OpticalSurface.cc 80747 2014-05-09 12:53:43Z gcosmo $
//
// 
////////////////////////////////////////////////////////////////////////
// Optical Surface Class Implementation
////////////////////////////////////////////////////////////////////////
//
// File:        G4OpticalSurface.cc
// Description: An optical surface class for use in G4OpBoundaryProcess
// Version:     2.0
// Created:     1997-06-26
// Author:      Peter Gumplinger
// mail:        gum@triumf.ca
//
////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <fstream>

//#include "G4ios.hh"
#include "globals.hh"
#include "G4OpticalSurface.hh"

/////////////////////////
// Class Implementation
/////////////////////////

        //////////////
        // Operators
        //////////////

G4OpticalSurface& G4OpticalSurface::operator=(const G4OpticalSurface& right)
{
  if (this != &right)
    {
      theName                    = right.theName;
      theType                    = right.theType;
      theModel                   = right.theModel;
      theFinish                  = right.theFinish;
      sigma_alpha                = right.sigma_alpha;
      polish                     = right.polish;
      theMaterialPropertiesTable = right.theMaterialPropertiesTable;
      if (AngularDistribution) delete [] AngularDistribution;
      AngularDistribution =
                       new G4float[incidentIndexMax*thetaIndexMax*phiIndexMax];
      *(AngularDistribution)     = *(right.AngularDistribution);
      if (DichroicVector) delete DichroicVector;
      DichroicVector = new G4Physics2DVector();
      *DichroicVector            = *(right.DichroicVector);
     } 
  return *this;
}

        /////////////////
        // Constructors
        /////////////////

G4OpticalSurface::G4OpticalSurface(const G4String& name,
				   G4OpticalSurfaceModel model,
				   G4OpticalSurfaceFinish finish,
				   G4SurfaceType type,
				   G4double value)
                                   : G4SurfaceProperty(name,type),
				     theModel(model),
				     theFinish(finish),
				     theMaterialPropertiesTable(0)
{
	if (model == glisur ){
		polish = value;
		sigma_alpha = 0.0;
	}
	else if ( model == unified ) {
		sigma_alpha = value;
		polish = 0.0;
	}
        else if ( model == LUT ) {
                sigma_alpha = value;
                polish = 0.0;
        }
        else if ( model == dichroic ) {
                sigma_alpha = value;
                polish = 0.0;
        }
	else {
                G4Exception("G4OpticalSurface::G4OpticalSurface()", "mat309",
                            FatalException,
			    "Constructor called with INVALID model.");
	}

        AngularDistribution = NULL;
        DichroicVector      = NULL;

        if (type == dielectric_LUT) {
           AngularDistribution =
                       new G4float[incidentIndexMax*thetaIndexMax*phiIndexMax];
           ReadLUTFile();
        }

        if (type == dielectric_dichroic) {
           DichroicVector = new G4Physics2DVector();
           ReadDichroicFile();
        }
}

G4OpticalSurface::~G4OpticalSurface()
{
        if (AngularDistribution) delete [] AngularDistribution;
        if (DichroicVector) delete DichroicVector;
}

G4OpticalSurface::G4OpticalSurface(const G4OpticalSurface &right)
  : G4SurfaceProperty(right.theName,right.theType)
{
       *this = right;
       this->theName = right.theName;
       this->theType = right.theType;
       this->theModel = right.theModel;
       this->theFinish = right.theFinish;
       this->sigma_alpha = right.sigma_alpha;
       this->polish = right.polish;
       this->theMaterialPropertiesTable = right.theMaterialPropertiesTable;
       if (this->AngularDistribution) delete [] AngularDistribution;
       this->AngularDistribution =
                       new G4float[incidentIndexMax*thetaIndexMax*phiIndexMax];
       *(this->AngularDistribution) = *(right.AngularDistribution);
       if (this->DichroicVector) delete DichroicVector;
       this->DichroicVector = new G4Physics2DVector();
       *(this->DichroicVector) = *(right.DichroicVector);
}

G4int G4OpticalSurface::operator==(const G4OpticalSurface &right) const
{
	return (this == (G4OpticalSurface *) &right);
}

G4int G4OpticalSurface::operator!=(const G4OpticalSurface &right) const
{
	return (this != (G4OpticalSurface *) &right);
}
        ////////////
        // Methods
        ////////////

void G4OpticalSurface::DumpInfo() const 
{

	// Dump info for surface

	G4cout << 
        "  Surface type   = " << G4int(theType)   << G4endl <<
	"  Surface finish = " << G4int(theFinish) << G4endl <<
	"  Surface model  = " << G4int(theModel)  << G4endl;

	G4cout << G4endl;

	G4cout << "  Surface parameter " << G4endl;
	G4cout << "  ----------------- " << G4endl;
	if (theModel == glisur ){
		G4cout << polish      << G4endl;
	}
        else if (theModel == LUT ){
                G4cout << sigma_alpha << G4endl;
        }
	else {
		G4cout << sigma_alpha << G4endl;
	}
	G4cout << G4endl;
}

void G4OpticalSurface::SetType(const G4SurfaceType& type)
{
  theType = type;
  if (type == dielectric_LUT) {
     if (!AngularDistribution) AngularDistribution =
                       new G4float[incidentIndexMax*thetaIndexMax*phiIndexMax];
     ReadLUTFile();
  }
  if (type == dielectric_dichroic) {
     if (!DichroicVector) DichroicVector = new G4Physics2DVector();
     ReadDichroicFile();
  }
}

void G4OpticalSurface::SetFinish(const G4OpticalSurfaceFinish finish)
{
  theFinish = finish;
  if (theType == dielectric_LUT) {
     if (!AngularDistribution) AngularDistribution =
                       new G4float[incidentIndexMax*thetaIndexMax*phiIndexMax];
     ReadLUTFile();
  }
  if (theType == dielectric_dichroic) {
     if (!DichroicVector) DichroicVector = new G4Physics2DVector();
     ReadDichroicFile();
  }
}

void G4OpticalSurface::ReadLUTFile()
{
  G4String readLUTFileName = " ";

  if (theFinish == polishedlumirrorglue) {
     readLUTFileName = "PolishedLumirrorGlue.dat";
  }
  else if (theFinish == polishedlumirrorair) {
     readLUTFileName = "PolishedLumirror.dat";
  }
  else if (theFinish == polishedteflonair) {
     readLUTFileName = "PolishedTeflon.dat";
  }
  else if (theFinish == polishedtioair) {
     readLUTFileName = "PolishedTiO.dat";
  }
  else if (theFinish == polishedtyvekair) {
     readLUTFileName = "PolishedTyvek.dat";
  }
  else if (theFinish == polishedvm2000glue) {
     readLUTFileName = "PolishedVM2000Glue.dat";
  }
  else if (theFinish == polishedvm2000air) {
     readLUTFileName = "PolishedVM2000.dat";
  }
  else if (theFinish == etchedlumirrorglue) {
     readLUTFileName = "EtchedLumirrorGlue.dat";
  }
  else if (theFinish == etchedlumirrorair) {
     readLUTFileName = "EtchedLumirror.dat";
  }
  else if (theFinish == etchedteflonair) {
     readLUTFileName = "EtchedTeflon.dat";
  }
  else if (theFinish == etchedtioair) {
     readLUTFileName = "EtchedTiO.dat";
  }
  else if (theFinish == etchedtyvekair) {
     readLUTFileName = "EtchedTyvek.dat";
  }
  else if (theFinish == etchedvm2000glue) {
     readLUTFileName = "EtchedVM2000Glue.dat";
  }
  else if (theFinish == etchedvm2000air) {
     readLUTFileName = "EtchedVM2000.dat";
  }
  else if (theFinish == groundlumirrorglue) {
     readLUTFileName = "GroundLumirrorGlue.dat";
  }
  else if (theFinish == groundlumirrorair) {
     readLUTFileName = "GroundLumirror.dat";
  }
  else if (theFinish == groundteflonair) {
     readLUTFileName = "GroundTeflon.dat";
  }
  else if (theFinish == groundtioair) {
     readLUTFileName = "GroundTiO.dat";
  }
  else if (theFinish == groundtyvekair) {
     readLUTFileName = "GroundTyvek.dat";
  }
  else if (theFinish == groundvm2000glue) {
     readLUTFileName = "GroundVM2000Glue.dat";
  }
  else if (theFinish == groundvm2000air) {
     readLUTFileName = "GroundVM2000.dat";
  }

  if (readLUTFileName == " ") return;

  char* path = getenv("G4REALSURFACEDATA");
  if (!path) {
     G4String excep =
        "G4OpBoundaryProcess - G4REALSURFACEDATA environment variable not set";
     G4Exception("G4OpticalSurface::ReadLUTFile()", "mat310",
                 FatalException, excep);
     return;
  }
  G4String pathString(path);

  readLUTFileName = pathString + "/" + readLUTFileName;

  std::ifstream readLUTFileHandle(readLUTFileName, std::ios::in);

  if (readLUTFileHandle) {
     G4int idxmax = incidentIndexMax*thetaIndexMax*phiIndexMax;
     for (G4int i = 0; i<idxmax; i++) {
       if (readLUTFileHandle.eof()) break;
       readLUTFileHandle >> AngularDistribution[i];
     }
     if (!readLUTFileHandle.bad()) {
        G4cout <<"LUT - data file: "<< readLUTFileName <<" read in! "<< G4endl;
     }
     else {
        G4String excep="LUT - data file: "+readLUTFileName+" not read propery";
        G4Exception("G4OpticalSurface::ReadLUTFile()", "mat312",
                    FatalException, excep);
        return;
     }
  }
  else {
     G4String excep ="LUT - data file: "+readLUTFileName+" not found";
     G4Exception("G4OpticalSurface::ReadLUTFile()", "mat311",
                 FatalException, excep);
     return;
  }
  readLUTFileHandle.close();
}

void G4OpticalSurface::ReadDichroicFile()
{
  const char* datadir = getenv("G4DICHROICDATA");

  if(!datadir) {
    G4Exception("G4OpticalSurface::ReadDichroicFile()","mat313",
    FatalException,"Environment variable G4DICHROICDATA not defined");
    return;
  }

  std::ostringstream ost;
  ost << datadir;
  std::ifstream fin(ost.str().c_str());
  if( !fin.is_open()) {
    G4ExceptionDescription ed;
    ed << "Dichroic surface data file <" << ost.str().c_str()
       << "> is not opened!" << G4endl;
    G4Exception("G4OpticalSurface::ReadDichroicFile()","mat314",
    FatalException,ed," ");
    return;
  }

  if( !(DichroicVector->Retrieve(fin)) ) {
    G4ExceptionDescription ed;
    ed << "Dichroic surface data file <" << ost.str().c_str()
       << "> is not opened!" << G4endl;
    G4Exception("G4OpticalSurface::ReadDichroicFile()","mat315",
    FatalException,ed," ");
    return;
  }

//  DichroicVector->SetBicubicInterpolation(true);

  G4cout << " *** Dichroic surface data file *** " << G4endl;

  G4int numberOfXNodes = DichroicVector->GetLengthX();
  G4int numberOfYNodes = DichroicVector->GetLengthY();

  G4cout << "numberOfXNodes: " << numberOfXNodes << G4endl;
  G4cout << "numberOfYNodes: " << numberOfYNodes << G4endl;

  if (0 > numberOfXNodes || numberOfXNodes >= INT_MAX) numberOfXNodes = 0;
  if (0 > numberOfYNodes || numberOfYNodes >= INT_MAX) numberOfYNodes = 0;

  G4PV2DDataVector  xVector;
  G4PV2DDataVector  yVector;

  xVector.resize(numberOfXNodes,0.);
  yVector.resize(numberOfYNodes,0.);

  for(G4int i = 0; i<numberOfXNodes; ++i) {
     G4cout << "i: " << DichroicVector->GetX(i) << G4endl;
     xVector[i] = DichroicVector->GetX(i);
  }
  for(G4int j = 0; j<numberOfYNodes; ++j) {
     G4cout << "j: " << DichroicVector->GetY(j) << G4endl;
     yVector[j] = DichroicVector->GetY(j);
  }

  for(G4int j = 0; j<numberOfYNodes; ++j) {
     for(G4int i = 0; i<numberOfXNodes; ++i) {
        G4cout << " i: " << i << " j: " << j << " "
               << DichroicVector->GetValue(i,j) << G4endl;
     }
  }

//  G4int idx, idy;

//  for(G4int j = 0; j<numberOfYNodes-1; ++j) {
//     G4double y = (yVector[j] + yVector[j+1])/2.;
//     for(G4int i = 0; i<numberOfXNodes-1; ++i) {
//        G4double x = (xVector[i] + xVector[i+1])/2.;
//        G4cout << " x: " << x << " y: " << y << " "
//               << DichroicVector->Value(x,y,idx,idy) << G4endl;
//     }
//  }

}
