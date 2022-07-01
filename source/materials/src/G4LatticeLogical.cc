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
/// \file materials/src/G4LatticeLogical.cc
/// \brief Implementation of the G4LatticeLogical class
//
//
#include "G4LatticeLogical.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include <cmath>
#include <fstream>


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4LatticeLogical::G4LatticeLogical()
  : verboseLevel(0), fVresTheta(0), fVresPhi(0), fDresTheta(0), fDresPhi(0),
    fA(0), fB(0), fLDOS(0), fSTDOS(0), fFTDOS(0),
    fBeta(0), fGamma(0), fLambda(0), fMu(0) {
  for (G4int i=0; i<3; i++) {
    for (G4int j=0; j<MAXRES; j++) {
      for (G4int k=0; k<MAXRES; k++) {
	fMap[i][j][k] = 0.;
	fN_map[i][j][k].set(0.,0.,0.);
      }
    }
  }
}

G4LatticeLogical::~G4LatticeLogical() {;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


///////////////////////////////////////////
//Load map of group velocity scalars (m/s)
////////////////////////////////////////////
G4bool G4LatticeLogical::LoadMap(G4int tRes, G4int pRes,
				 G4int polarizationState, G4String map) {
  if (tRes>MAXRES || pRes>MAXRES) {
    G4cerr << "G4LatticeLogical::LoadMap exceeds maximum resolution of "
           << MAXRES << " by " << MAXRES << ". terminating." << G4endl;
    return false; 		//terminate if resolution out of bounds.
  }

  std::ifstream fMapFile(map.data());
  if(!fMapFile.is_open())
  {
    return false;
  }

  G4double vgrp = 0.;
  for (G4int theta = 0; theta<tRes; theta++) {
    for (G4int phi = 0; phi<pRes; phi++) {
      fMapFile >> vgrp;
      fMap[polarizationState][theta][phi] = vgrp * (m/s);
    }
  }

  if(verboseLevel != 0)
  {
    G4cout << "\nG4LatticeLogical::LoadMap(" << map << ") successful"
	   << " (Vg scalars " << tRes << " x " << pRes << " for polarization "
	   << polarizationState << ")." << G4endl;
  }

  fVresTheta=tRes; //store map dimensions
  fVresPhi=pRes;
  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


////////////////////////////////////
//Load map of group velocity unit vectors
///////////////////////////////////
G4bool G4LatticeLogical::Load_NMap(G4int tRes, G4int pRes,
				   G4int polarizationState, G4String map) {
  if (tRes>MAXRES || pRes>MAXRES) {
    G4cerr << "G4LatticeLogical::LoadMap exceeds maximum resolution of "
           << MAXRES << " by " << MAXRES << ". terminating." << G4endl;
    return false; 		//terminate if resolution out of bounds.
  }

  std::ifstream fMapFile(map.data());
  if(!fMapFile.is_open())
  {
    return false;
  }

  G4double x,y,z;	// Buffers to read coordinates from file
  G4ThreeVector dir;
  for (G4int theta = 0; theta<tRes; theta++) {
    for (G4int phi = 0; phi<pRes; phi++) {
      fMapFile >> x >> y >> z;
      dir.set(x,y,z);
      fN_map[polarizationState][theta][phi] = dir.unit();	// Enforce unity
    }
  }

  if(verboseLevel != 0)
  {
    G4cout << "\nG4LatticeLogical::Load_NMap(" << map << ") successful"
	   << " (Vdir " << tRes << " x " << pRes << " for polarization "
	   << polarizationState << ")." << G4endl;
  }

  fDresTheta=tRes; //store map dimensions
  fDresPhi=pRes;
  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

//Given the phonon wave vector k, phonon physical volume Vol 
//and polarizationState(0=LON, 1=FT, 2=ST), 
//returns phonon velocity in m/s

G4double G4LatticeLogical::MapKtoV(G4int polarizationState,
				   const G4ThreeVector& k) const {
  G4double theta, phi, tRes, pRes;

  tRes=pi/fVresTheta;
  pRes=twopi/fVresPhi;
  
  theta=k.getTheta();
  phi=k.getPhi();

  if(phi < 0)
  {
    phi = phi + twopi;
  }
  if(theta > pi)
  {
    theta = theta - pi;
  }

  G4double Vg = fMap[polarizationState][int(theta/tRes)][int(phi/pRes)];

  if(Vg == 0){
      G4cout<<"\nFound v=0 for polarization "<<polarizationState
            <<" theta "<<theta<<" phi "<<phi<< " translating to map coords "
            <<"theta "<< int(theta/tRes) << " phi " << int(phi/pRes)<<G4endl;
  }

  if (verboseLevel>1) {
    G4cout << "G4LatticeLogical::MapKtoV theta,phi=" << theta << " " << phi
	   << " : ith,iph " << int(theta/tRes) << " " << int(phi/pRes)
	   << " : V " << Vg << G4endl;
  }

  return Vg;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

//Given the phonon wave vector k, phonon physical volume Vol 
//and polarizationState(0=LON, 1=FT, 2=ST), 
//returns phonon propagation direction as dimensionless unit vector

G4ThreeVector G4LatticeLogical::MapKtoVDir(G4int polarizationState,
					   const G4ThreeVector& k) const {  
  G4double theta, phi, tRes, pRes;

  tRes=pi/(fDresTheta-1);//The summant "-1" is required:index=[0:array length-1]
  pRes=2*pi/(fDresPhi-1);

  theta=k.getTheta();
  phi=k.getPhi();

  if(theta > pi)
  {
    theta = theta - pi;
  }
  //phi=[0 to 2 pi] in accordance with DMC //if(phi>pi/2) phi=phi-pi/2;
  if(phi < 0)
  {
    phi = phi + 2 * pi;
  }

  G4int iTheta = int(theta/tRes+0.5);
  G4int iPhi = int(phi/pRes+0.5);

  if (verboseLevel>1) {
    G4cout << "G4LatticeLogical::MapKtoVDir theta,phi=" << theta << " " << phi
	   << " : ith,iph " << iTheta << " " << iPhi
	   << " : dir " << fN_map[polarizationState][iTheta][iPhi] << G4endl;
  }

  return fN_map[polarizationState][iTheta][iPhi];
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

// Dump structure in format compatible with reading back

void G4LatticeLogical::Dump(std::ostream& os) const {
  os << "dyn " << fBeta << " " << fGamma << " " << fLambda << " " << fMu
     << "\nscat " << fB << " decay " << fA
     << "\nLDOS " << fLDOS << " STDOS " << fSTDOS << " FTDOS " << fFTDOS
     << std::endl;

  Dump_NMap(os, 0, "LVec.ssv");
  Dump_NMap(os, 1, "FTVec.ssv");
  Dump_NMap(os, 2, "STVec.ssv");

  DumpMap(os, 0, "L.ssv");
  DumpMap(os, 1, "FT.ssv");
  DumpMap(os, 2, "ST.ssv");
}

void G4LatticeLogical::DumpMap(std::ostream& os, G4int pol,
			       const G4String& name) const {
  os << "VG " << name << " " << (pol==0?"L":pol==1?"FT":pol==2?"ST":"??")
     << " " << fVresTheta << " " << fVresPhi << std::endl;

  for (G4int iTheta=0; iTheta<fVresTheta; iTheta++) {
    for (G4int iPhi=0; iPhi<fVresPhi; iPhi++) {
      os << fMap[pol][iTheta][iPhi] << std::endl;
    }
  }
}

void G4LatticeLogical::Dump_NMap(std::ostream& os, G4int pol,
				 const G4String& name) const {
  os << "VDir " << name << " " << (pol==0?"L":pol==1?"FT":pol==2?"ST":"??")
     << " " << fDresTheta << " " << fDresPhi << std::endl;

  for (G4int iTheta=0; iTheta<fDresTheta; iTheta++) {
    for (G4int iPhi=0; iPhi<fDresPhi; iPhi++) {
      os << fN_map[pol][iTheta][iPhi].x()
	 << " " << fN_map[pol][iTheta][iPhi].y()
	 << " " << fN_map[pol][iTheta][iPhi].z()
	 << std::endl;
    }
  }
}

