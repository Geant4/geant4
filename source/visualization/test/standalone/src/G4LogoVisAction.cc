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
// $Id: G4LogoVisAction.cc,v 1.1 2005-03-03 16:45:34 allison Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

#include "G4LogoVisAction.hh"

#include "G4VVisManager.hh"
#include "G4VisAttributes.hh"
#include "G4Polyhedron.hh"
#include "G4Tubs.hh"
#include "G4Box.hh"
#include "G4UnionSolid.hh"
#include "G4SubtractionSolid.hh"

void G4LogoVisAction::Draw() {
  G4VVisManager* pVisManager = G4VVisManager::GetConcreteInstance();
  if (pVisManager) {
    const G4double height = 1*m;
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

    G4VisAttributes visAtts;
    visAtts.SetColour(G4Colour(0.,1.,0.));  // Green.
    visAtts.SetForceSolid(true);         // Always solid.

    // G...
    G4Tubs logo1("logo1",ri,ro,d2,0.15*pi,1.85*pi);
    G4Box  logo2("logo2",w2,ro2,d2);
    G4UnionSolid logoG("logoG",&logo1,&logo2,G4Translate3D(ri+w2,-ro2,0.));
    pVisManager->Draw(logoG,visAtts,G4Translate3D(-1.1*ro,0.,0.));

    // 4...
    G4Box logo3("logo3",h2,h2,d2);
    G4Box logoS("logoS",s,s,d2+e);  // Subtractor.
    G4SubtractionSolid logo5("logo5",&logo3,&logoS,
			      G4Translate3D(f1-s,f2-s,0.));
    G4SubtractionSolid logo6("logo6",&logo5,&logoS,
			      G4Translate3D(f1+s+w,f2-s,0.));
    G4SubtractionSolid logo7("logo7",&logo6,&logoS,
			      G4Translate3D(f1+s+w,f2+s+w,0.));
    G4SubtractionSolid logo8("logo8",&logo7,&logoS,
			     G4Transform3D(rm,G4ThreeVector(x8,y8,0.)));
    G4SubtractionSolid logo9("logo9",&logoS,&logoS,  // Triangular hole.
			      G4Transform3D(rm,G4ThreeVector(x9,y9,0.)));
    G4SubtractionSolid logo4("logo4",&logo8,&logo9,
			      G4Translate3D(-xtr,-ytr,0.));
    pVisManager->Draw(logo4,visAtts,G4Translate3D(1.1*ro,0.,0.));
  }
}
