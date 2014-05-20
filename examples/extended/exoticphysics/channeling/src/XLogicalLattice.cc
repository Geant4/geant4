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

#include "XLogicalLattice.hh"
#include "G4PhysicalConstants.hh"
#include "Randomize.hh"
#include <cmath>

XLogicalLattice::XLogicalLattice(){
  fVresTheta=0;
  fVresPhi=0;
  fDresTheta=0;
  fDresPhi=0;
  fA=0;
  fB=0;
  fDosL=0;
  fDosST=0;
  fDosFT=0;
  fBeta=0;
  fLambda=0;
  fGamma=0;
  fMu=0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


XLogicalLattice::~XLogicalLattice(){;}


void XLogicalLattice::SetDynamicalConstants(double Beta,
                                            double Gamma,
                                            double Lambda,
                                            double Mu)
{
  fBeta=Beta;
  fGamma=Gamma;
  fLambda=Lambda;
  fMu=Mu;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


void XLogicalLattice::SetScatteringConstant(G4double b)
{
  //Constant governing the rate of isotope scattering, for use with
  //XPhononScatteringProcess
  fB=b;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


void XLogicalLattice::SetAnhDecConstant(G4double a)
{
  //Constant governing rate of anharmonic down conversion of L-phonons,
  //for use with XPhononDownconversionProcess
  fA=a;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XLogicalLattice::SetLDOS(double LDOS)
{
  //Longitudinal phonon density of states
  fDosL=LDOS;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


void XLogicalLattice::SetSTDOS(double STDOS)
{
  //Slow transverse phonon density of states
  fDosST=STDOS;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


void XLogicalLattice::SetFTDOS(double FTDOS)
{
  //Fast transverse phonon density of states
  fDosFT=FTDOS;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


double XLogicalLattice::GetBeta()
{
  return fBeta;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


double XLogicalLattice::GetGamma()
{
  return fGamma;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


double XLogicalLattice::GetLambda()
{
  return fLambda;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


double XLogicalLattice::GetMu()
{
  return fMu;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


G4double XLogicalLattice::GetScatteringConstant()
{
  return fB;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


G4double XLogicalLattice::GetAnhDecConstant()
{
  return fA;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


double XLogicalLattice::GetLDOS()
{
  return fDosL;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


double XLogicalLattice::GetSTDOS()
{
  return fDosST;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


double XLogicalLattice::GetFTDOS()
{
  return fDosFT;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


///////////////////////////////////////////
//Load map of group velocity scalars (m/s)
////////////////////////////////////////////
bool XLogicalLattice::LoadMap(int tRes,
                              int pRes,
                              int polarizationState,
                              string m)
{
  if((tRes>MAXRES)||(pRes>MAXRES)){
    G4cout<<"\nk-V fMap exceeds maximum resolution of "<<MAXRES<<
      " by "<<MAXRES<<". terminating."<<endl;
    return false; //terminate if resolution out of bounds.
  }


  fMapFile.clear();

  fThetaRes=tRes;
  fPhiRes=pRes;
  fMapFile.open(m.data());
  if(!fMapFile.is_open()) return false;
  
  for(int theta = 0; theta<fThetaRes; theta++){
    for(int phi = 0; phi<fPhiRes; phi++){
        fMapFile>>fMap[polarizationState][theta][phi];        
    }
  }
  fMapFile.close();
  G4cout<<"\nXLogicalLattice::LoadMap() sucessful (Vg scalars)."<<endl;
  fVresTheta=tRes; //store map dimensions
  fVresPhi=pRes;
  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


////////////////////////////////////
//Load map of group velocity unit vectors
///////////////////////////////////
bool XLogicalLattice::Load_NMap(int tRes,
                                int pRes,
                                int polarizationState,
                                string m)
{
  if((tRes>MAXRES)||(pRes>MAXRES)){
    G4cout<<"\nk-V map exceeds maximum resolution of "<<MAXRES<<
      " by "<<MAXRES<<" terminating."<<endl;
    return false; //terminate if resolution out of bounds.
  }


  fMapFile.clear();

  fThetaRes=tRes;
  fPhiRes=pRes;
  fMapFile.open(m.data());
  if(!fMapFile.is_open()) return false;
  
  for(int theta = 0; theta<fThetaRes; theta++){
    for(int phi = 0; phi<fPhiRes; phi++){
      for(int coord = 0; coord<3; coord++){
        fMapFile>>fN_map[polarizationState][theta][phi][coord];
      }
    }
  }
  G4cout<<"\n";
  fMapFile.close();
  G4cout<<"\nXLogicalLattice::Load_NMap() sucessful"<<endl;
  fDresTheta=tRes; //store map dimesnions
  fDresPhi=pRes;
  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....



double XLogicalLattice::MapKtoV(int polarizationState,G4ThreeVector k)
{
  //Given the phonon wave vector k, phonon physical volume Vol 
  //and polarizationState(0=LON, 1=FT, 2=ST), 
  //returns phonon velocity in m/s

  double theta, phi, tRes, pRes;

  tRes=pi/(fVresTheta);//The summant "-1" is required:index=[0:array length-1]
  pRes=2*pi/(fVresPhi);

  
  theta=k.getTheta();
  phi=k.getPhi();

  if(phi<0) phi = phi + 2*pi;
  if(theta>pi) theta=theta-pi;
  //phi=[0 to 2 pi] in accordance with DMC //if(phi>pi/2) phi=phi-pi/2;
  if(fMap[polarizationState][int(theta/tRes)][int(phi/pRes)]==0){
      G4cout<<"\nFound v=0 for polarization "<<polarizationState
      <<" theta "<<theta<<" phi "<<phi<< " translating to map coords "
      << "theta "<< int(theta/tRes) << " phi " << int(phi/pRes)<<endl;
  }

  return fMap[polarizationState][int(theta/tRes)][int(phi/pRes)];
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


G4ThreeVector XLogicalLattice::MapKtoVDir(int polarizationState,G4ThreeVector k)
{  
  //Given the phonon wave vector k, phonon physical volume Vol 
  //and polarizationState(0=LON, 1=FT, 2=ST), 
  //returns phonon propagation direction as dimensionless unit vector

  double theta, phi, tRes, pRes;

  tRes=pi/(fDresTheta-1);//The summant "-1" is required:index=[0:array length-1]
  pRes=2*pi/(fDresPhi-1);

  theta=k.getTheta();
  phi=k.getPhi(); 

  if(theta>pi) theta=theta-pi;
  //phi=[0 to 2 pi] in accordance with DMC //if(phi>pi/2) phi=phi-pi/2;
  if(phi<0) phi = phi + 2*pi;

  G4int iTheta = int(theta/tRes+0.5);
  G4int iPhi = int(phi/pRes+0.5);

  
  G4ThreeVector v(fN_map[polarizationState][iTheta][iPhi][0],
                  fN_map[polarizationState][iTheta][iPhi][1],
                  fN_map[polarizationState][iTheta][iPhi][2]);


  //////debugging purposes only//////////////
  //v.Set(sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta));
  //v.Set(k.getX(), k.getY(), k.getZ());

  return v.unit();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

