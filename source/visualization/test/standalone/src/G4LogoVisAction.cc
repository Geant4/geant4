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
// $Id: G4LogoVisAction.cc,v 1.2 2005-03-16 17:25:15 allison Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

#include "G4LogoVisAction.hh"

#include "G4VVisManager.hh"
#include "G4VisAttributes.hh"
#include "G4Polyhedron.hh"
#include "G4Tubs.hh"
#include "G4Box.hh"
#include "G4UnionSolid.hh"
#include "G4SubtractionSolid.hh"

G4LogoVisAction::G4LogoVisAction() {

  fpVisAtts = new G4VisAttributes(G4Colour(0.,1.,0.));  // Green.
  fpVisAtts->SetForceSolid(true);                      // Always solid.
  // The vis attributes object has to have a life at least as long as
  // the objects to which it is attached; it may not be a
  // local/automatic object.

  const G4double height = 1.*m;
  const G4double& h =  height;
  const G4double h2  = 0.5 * h;   // Half height.
  const G4double ri  = 0.25 * h;  // Inner radius.
  const G4double ro  = 0.5 * h;   // Outer radius.
  const G4double ro2 = 0.5 * ro;  // Half outer radius.
  const G4double w   = ro - ri;   // Width.
  const G4double w2  = 0.5 * w;   // Half width.
  const G4double d2  = 0.2 * h;   // Half depth.
  const G4double f1  = 0.05 * h;  // left edge of stem of "4".
  const G4double f2  = -0.3 * h;  // bottom edge of cross of "4".
  const G4double e = 1.e-4 * h;   // epsilon.
  const G4double xt = f1, yt = h2;      // Top of slope.
  const G4double xb = -h2, yb = f2 + w; // Bottom of slope.
  const G4double dx = xt - xb, dy = yt - yb;
  const G4double angle = atan2(dy,dx);
  G4RotationMatrix rm;
  rm.rotateZ(angle*rad);
  const G4double d = sqrt(dx * dx + dy * dy);
  const G4double s = h;  // Half height of square subtractor
  const G4double y8 = s; // Choose y of subtractor for outer slope.
  const G4double x8 = ((-s * d - dx * (yt - y8)) / dy) + xt;
  G4double y9 = s; // Choose y of subtractor for inner slope.
  G4double x9 = ((-(s - w) * d - dx * (yt - y8)) / dy) + xt;
  // But to get inner, we make a triangle translated by...
  const G4double xtr = s - f1, ytr = -s - f2 -w;
  x9 += xtr; y9 += ytr;

  // G...
  G4Tubs tG("tG",ri,ro,d2,0.15*pi,1.85*pi);
  G4Box bG("bG",w2,ro2,d2);
  G4UnionSolid logoG("logoG",&tG,&bG,G4Translate3D(ri+w2,-ro2,0.));
  fpG = logoG.CreatePolyhedron();
  fpG->SetVisAttributes(fpVisAtts);
  fpG->Transform(G4Translate3D(-0.55*h,0.,0.));

  // 4...
  G4Box b1("b1",h2,h2,d2);
  G4Box bS("bS",s,s,d2+e);  // Subtractor.
  G4Box bS2("bS2",s,s,d2+2.*e);  // 2nd Subtractor.
  G4SubtractionSolid s1("s1",&b1,&bS,G4Translate3D(f1-s,f2-s,0.));
  G4SubtractionSolid s2("s2",&s1,&bS,G4Translate3D(f1+s+w,f2-s,0.));
  G4SubtractionSolid s3("s3",&s2,&bS,G4Translate3D(f1+s+w,f2+s+w,0.));
  G4SubtractionSolid s4
    ("s4",&s3,&bS,G4Transform3D(rm,G4ThreeVector(x8,y8,0.)));
  G4SubtractionSolid s5    // Triangular hole.
    ("s5",&bS,&bS2,G4Transform3D(rm,G4ThreeVector(x9,y9,0.)));
  G4SubtractionSolid logo4("logo4",&s4,&s5,G4Translate3D(-xtr,-ytr,0.));
  fp4 = logo4.CreatePolyhedron();
  /* Experiment with creating own polyhedron...
  int nNodes = 4;
  int nFaces = 4;
  double xyz[][3] = {{0,0,0},{1*m,0,0},{0,1*m,0},{0,0,1*m}};
  int faces[][4] = {{1,3,2,0},{1,2,4,0},{1,4,3,0},{2,3,4,0}};
  fp4 = new G4Polyhedron();
  fp4->createPolyhedron(nNodes,nFaces,xyz,faces);
  */
  fp4->SetVisAttributes(fpVisAtts);
  fp4->Transform(G4Translate3D(0.55*h,0.,0.));
}

G4LogoVisAction::~G4LogoVisAction() {
  delete fpG;
  delete fp4;
}

void G4LogoVisAction::Draw() {
  G4VVisManager* pVisManager = G4VVisManager::GetConcreteInstance();
  if (pVisManager) {
    pVisManager->Draw(*fpG);
    pVisManager->Draw(*fp4);
  }
}
