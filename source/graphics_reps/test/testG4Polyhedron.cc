// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: testG4Polyhedron.cc,v 1.2 1999-12-15 14:50:37 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

// Create polyhedron for different shapes

#include "G4ios.hh"
#include "G4Polyhedron.hh"

int main() {

  G4Polyhedron polyhedron;

//   B O X

  G4cout << "=== G4PolyhedronBox" << G4endl;
  polyhedron = G4PolyhedronBox(100., 200., 400.);

//   T R D 1

  G4cout << "=== G4PolyhedronTrd1" << G4endl;
  polyhedron = G4PolyhedronTrd1(100., 150., 200., 400.);

//   T R D 2

  G4cout << "=== G4PolyhedronTrd2" << G4endl;
  polyhedron = G4PolyhedronTrd2(2., 3., 4., 5., 6.);

//   P A R A

  G4cout << "=== G4PolyhedronPara" << G4endl;
  polyhedron = G4PolyhedronPara(100., 200., 400., 15.*deg, 30.*deg, 30.*deg);

//   T R A P

  G4cout << "=== G4PolyhedronTrap" << G4endl;
  polyhedron = G4PolyhedronTrap(390., 0.*deg, 0.*deg,
				60., 40., 90., 15.*deg,
				120., 80., 180., 15.*deg);

//   T U B E

  G4cout << "=== G4PolyhedronTube" << G4endl;
  polyhedron = G4PolyhedronTube(100., 200., 400.);

  G4cout << "=== G4PolyhedronTube(Rmin = 0)" << G4endl;
  polyhedron = G4PolyhedronTube(0., 200., 400.);

//   T U B S

  G4cout << "=== G4PolyhedronTubs" << G4endl;
  polyhedron = G4PolyhedronTubs(0., 200., 400., 200.*deg, 140.*deg);

  G4cout << "=== G4PolyhedronTubs(Rmin = 0, Dphi=180)" << G4endl;
  polyhedron = G4PolyhedronTubs(0., 200., 400., 200.*deg, 180.*deg);

//   C O N E

  G4cout << "=== G4PolyhedronCone" << G4endl;
  polyhedron = G4PolyhedronCone(50., 100., 150., 200., 400.);

  G4cout << "=== G4PolyhedronCone(Rmin1 = Rmin2 = 0)" << G4endl;
  polyhedron = G4PolyhedronCone(0., 100., 0., 200., 400.);

  G4cout << "=== G4PolyhedronCone(Rmin1 = Rmax1 = 0)" << G4endl;
  polyhedron = G4PolyhedronCone(0., 0., 150., 200., 400.);

  G4cout << "=== G4PolyhedronCone(Rmin1 = Rmax1)" << G4endl;
  polyhedron = G4PolyhedronCone(100., 100., 150., 200., 400.);

  G4cout << "=== G4PolyhedronCone(Rmin2 = Rmax2 = 0.)" << G4endl;
  polyhedron = G4PolyhedronCone(50., 100., 0., 0., 400.);

  G4cout << "=== G4PolyhedronCone(Rmin2 = Rmax2)" << G4endl;
  polyhedron = G4PolyhedronCone(50., 100., 200., 200., 400.);

//   C O N S

  G4cout << "=== G4PolyhedronCons" << G4endl;
  polyhedron = G4PolyhedronCons(50.,100.,150.,200.,400.,200.*deg,140.*deg);

  G4cout << "=== G4PolyhedronCons(Rmin1 = Rmin2 = 0)" << G4endl;
  polyhedron = G4PolyhedronCons(0.,100.,0.,200.,400.,200.*deg,140.*deg);

  G4cout << "=== G4PolyhedronCons(Rmin1 = Rmin2 = 0, Dphi=180)" << G4endl;
  polyhedron = G4PolyhedronCons(0.,100.,0.,200.,400.,200.*deg,180.*deg);

  G4cout << "=== G4PolyhedronCons(Rmin1 = Rmax1 = 0)" << G4endl;
  polyhedron = G4PolyhedronCons(0.,0.,150.,200.,400.,200.*deg,180.*deg);

  G4cout << "=== G4PolyhedronCons(Rmin1 = Rmax1)" << G4endl;
  polyhedron = G4PolyhedronCons(100.,100.,150.,200.,400.,200.*deg,180.*deg);

  G4cout << "=== G4PolyhedronCons(Rmin2 = Rmax2 = 0.)" << G4endl;
  polyhedron = G4PolyhedronCons(50.,100.,0.,0.,400.,200.*deg,180.*deg);

  G4cout << "=== G4PolyhedronCons(Rmin2 = Rmax2)" << G4endl;
  polyhedron = G4PolyhedronCons(50.,100.,200.,200.,400.,200.*deg,180.*deg);

//   S P H E R E

  G4cout << "=== G4PolyhedronSphere" << G4endl;
  polyhedron = G4PolyhedronSphere(100.,200.,0.*deg,360.*deg,0.*deg,180.*deg);

  G4cout << "=== G4PolyhedronSphere(Rmin=0)" << G4endl;
  polyhedron = G4PolyhedronSphere(0.,200.,0.*deg,360.*deg,0.*deg,180.*deg);

  G4cout << "=== G4PolyhedronSphere(Dphi=180)" << G4endl;
  polyhedron = G4PolyhedronSphere(100.,200.,0.*deg,180.*deg,0.*deg,180.*deg);

  G4cout << "=== G4PolyhedronSphere(Rmin=0, Dphi=180)" << G4endl;
  polyhedron = G4PolyhedronSphere(0.,200.,5.*deg,180.*deg,0.*deg,180.*deg);

  G4cout << "=== G4PolyhedronSphere(Dthe=0-90)" << G4endl;
  polyhedron = G4PolyhedronSphere(100.,200.,0.*deg,360.*deg,0.*deg,90.*deg);

  G4cout << "=== G4PolyhedronSphere(Rmin=0, Dthe=0-90)" << G4endl;
  polyhedron = G4PolyhedronSphere(0.,200.,5.*deg,360.*deg,0.*deg,90.*deg);

  G4cout << "=== G4PolyhedronSphere(Dthe=90-180)" << G4endl;
  polyhedron = G4PolyhedronSphere(100.,200.,0.*deg,360.*deg,90.*deg,90.*deg);

  G4cout << "=== G4PolyhedronSphere(Rmin=0, Dthe=90-180)" << G4endl;
  polyhedron = G4PolyhedronSphere(0.,200.,5.*deg,360.*deg,90.*deg,90.*deg);

  G4cout << "=== G4PolyhedronSphere(Dphi=180, Dthe=0-90)" << G4endl;
  polyhedron = G4PolyhedronSphere(100.,200.,5.*deg,180.*deg,0.*deg,90.*deg);

  G4cout << "=== G4PolyhedronSphere(Rmin=0, Dphi=180, Dthe=0-900)" << G4endl;
  polyhedron = G4PolyhedronSphere(0.,200.,5.*deg,180.*deg,0.*deg,90.*deg);

  G4cout << "=== G4PolyhedronSphere(Dphi=180, Dthe=90-180)" << G4endl;
  polyhedron = G4PolyhedronSphere(100.,200.,5.*deg,180.*deg,90.*deg,90.*deg);

  G4cout << "=== G4PolyhedronSphere(Rmin=0, Dphi=180, Dthe=90-180)" << G4endl;
  polyhedron = G4PolyhedronSphere(0.,200.,5.*deg,180.*deg,90.*deg,90.*deg);

  G4cout << "=== G4PolyhedronSphere(Dphi=45-135, Dthe=45-135)" << G4endl;
  polyhedron = G4PolyhedronSphere(100.,200.,45.*deg,90.*deg,45.*deg,90.*deg);

  G4cout << "=== G4PolyhedronSphere(Rmin=0, Dphi=30-120, Dthe=30-120)" << G4endl;
  polyhedron = G4PolyhedronSphere(0.,200.,30.*deg,90.*deg,30.*deg,90.*deg);

//   T O R U S

  G4cout << "=== G4PolyhedronTorus" << G4endl;
  polyhedron = G4PolyhedronTorus(100.,200.,400.,0.*deg,360.*deg);
  
  G4cout << "=== G4PolyhedronTorus(Rmin=0)" << G4endl;
  polyhedron = G4PolyhedronTorus(0.,200.,400.,0.*deg,360.*deg);
  
  G4cout << "=== G4PolyhedronTorus(Dphi=180)" << G4endl;
  polyhedron = G4PolyhedronTorus(100.,200.,400.,5.*deg,180.*deg);
  
  G4cout << "=== G4PolyhedronTorus(Rmin=0, Dphi=180)" << G4endl;
  polyhedron = G4PolyhedronTorus(0.,200.,400.,5.*deg,180.*deg);
  
//   P G O N

  G4cout << "=== G4PolyhedronPgon(Nz=4)" << G4endl;
  G4double rmax01[4] = {  150.,  200., 200., 150.};
  G4double rmin01[4] = {   50.,  100., 100.,  50.};
  G4double z01[4]    = { -200., -100., 100., 200.};
  polyhedron = G4PolyhedronPgon(5.*deg, 45.*deg, 2, 4, z01, rmin01, rmax01);

  G4cout << "=== G4PolyhedronPgon(Nz=4, N=1)" << G4endl;
  polyhedron = G4PolyhedronPgon(5.*deg, 45.*deg, 1, 4, z01, rmin01, rmax01);

  G4cout << "=== G4PolyhedronPgon(Nz=4, Rmin=Rmax)" << G4endl;
  G4double rmax02[4] = {  150.,  200., 200., 150.};
  G4double rmin02[4] = {  150.,  100., 100., 150.};
  G4double z02[4]    = { -200., -100., 100., 200.};
  polyhedron = G4PolyhedronPgon(5.*deg, 45.*deg, 2, 4, z02, rmin02, rmax02);

  G4cout << "=== G4PolyhedronPgon(Nz=4, N=1, Rmin=Rmax)" << G4endl;
  polyhedron = G4PolyhedronPgon(5.*deg, 45.*deg, 1, 4, z02, rmin02, rmax02);

  G4cout << "=== G4PolyhedronPgon(Nz=4, Rmin=Rmax=0)" << G4endl;
  G4double rmax03[4] = {  0.,  200., 200., 0.};
  G4double rmin03[4] = {  0.,  100., 100., 0.};
  G4double z03[4]    = { -200., -100., 100., 200.};
  polyhedron = G4PolyhedronPgon(5.*deg, 45.*deg, 2, 4, z03, rmin03, rmax03);

  G4cout << "=== G4PolyhedronPgon(Nz=4, N=1, Rmin=Rmax=0)" << G4endl;
  polyhedron = G4PolyhedronPgon(5.*deg, 45.*deg, 1, 4, z03, rmin03, rmax03);

  G4cout << "=== G4PolyhedronPgon(Nz=4, Rmin=100)" << G4endl;
  G4double rmax04[4] = {  150.,  200., 200., 150.};
  G4double rmin04[4] = {  100.,  100., 100., 100.};
  G4double z04[4]    = { -200., -100., 100., 200.};
  polyhedron = G4PolyhedronPgon(5.*deg, 120.*deg, 6, 4, z04, rmin04, rmax04);

  G4cout << "=== G4PolyhedronPgon(Nz=4, Rmin=0, Dphi=180)" << G4endl;
  G4double rmax05[4] = {  150.,  200., 200., 150.};
  G4double rmin05[4] = {    0.,    0.,   0.,   0.};
  G4double z05[4]    = { -200., -100., 100., 200.};
  polyhedron = G4PolyhedronPgon(5.*deg, 180.*deg, 6, 4, z05, rmin05, rmax05);

  G4cout << "=== G4PolyhedronPgon(Nz=3, Rmin=Rmax)" << G4endl;
  G4double rmax06[3] = {  100.,  200., 200.};
  G4double rmin06[3] = {  100.,  100., 100.};
  G4double z06[3]    = { -200., -100., 100.};
  polyhedron = G4PolyhedronPgon(5.*deg, 120.*deg, 6, 3, z06, rmin06, rmax06);

//   P C O N

  G4cout << "=== G4PolyhedronPcon(Nz=4)" << G4endl;
  polyhedron = G4PolyhedronPcon(5.*deg, 45.*deg, 4, z01, rmin01, rmax01);

  G4cout << "=== G4PolyhedronPcon(Nz=4, Rmin=Rmax)" << G4endl;
  polyhedron = G4PolyhedronPcon(5.*deg, 45.*deg, 4, z02, rmin02, rmax02);

  G4cout << "=== G4PolyhedronPcon(NZ=4, Rmin=Rmax=0)" << G4endl;
  polyhedron = G4PolyhedronPcon(5.*deg, 45.*deg, 4, z03, rmin03, rmax03);

  G4cout << "=== G4PolyhedronPcon(Nz=4, Rmin=100)" << G4endl;
  polyhedron = G4PolyhedronPcon(5.*deg, 120.*deg, 4, z04, rmin04, rmax04);

  G4cout << "=== G4PolyhedronPcon(Nz=4, Rmin=0, Dphi=180)" << G4endl;
  polyhedron = G4PolyhedronPcon(5.*deg, 180.*deg, 4, z05, rmin05, rmax05);

  G4cout << "=== G4PolyhedronPcon(Nz=3, Rmin=Rmax)" << G4endl;
  polyhedron = G4PolyhedronPcon(5.*deg, 120.*deg, 3, z06, rmin06, rmax06);

  return 0;
}
     
      
