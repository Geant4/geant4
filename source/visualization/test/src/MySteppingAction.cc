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
// $Id: MySteppingAction.cc,v 1.16 2006-06-29 21:34:34 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// 16/Apr/1997  J. Allison:  For visualization/test/test19.

#include "MySteppingAction.hh"

#include "G4SteppingManager.hh"
#include "G4ParticleDefinition.hh"
#include "G4VVisManager.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4Polyline.hh"
#include "G4Circle.hh"
#include "G4Square.hh"

void MySteppingAction::UserSteppingAction(const G4Step* pStep) {
  // User Action Example - begin snippet.
  static int coutCount = 0;
  if (coutCount < 10) {
    coutCount++;
    G4cout << "MySteppingAction::UserSteppingAction called." << G4endl;
  }

  G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
  if(pVVisManager) {
    G4Polyline pl;
    G4Colour c;
    G4double chg =
      pStep -> GetTrack () -> GetDefinition () -> GetPDGCharge ();
    if (chg != 0.) c = G4Colour (1., 0., 0.);
    else           c = G4Colour (0., 1., 0.);
    G4VisAttributes va (c);
    pl.SetVisAttributes (&va);
    pl.push_back (pStep -> GetPreStepPoint () -> GetPosition ());
    pl.push_back (pStep -> GetPostStepPoint () -> GetPosition ());
    static int coutCount = 0;
    if (coutCount < 10) {
      coutCount++;
      G4cout << "MySteppingAction::UserSteppingAction: drawing "
	   <<  pl << G4endl;
    }
    pVVisManager -> Draw (pl);
    G4Circle intercept (pStep -> GetPostStepPoint () -> GetPosition ());
    G4VisAttributes iva (G4Colour(1.,0.,1));
    intercept.SetVisAttributes (&iva);
    pVVisManager -> Draw (intercept);

    const G4double pixels = 1.;

    G4Circle circle22(G4Point3D(100.*cm,100.*cm,0.));
    circle22.SetScreenDiameter(50.*pixels);
    circle22.SetFillStyle(G4Circle::noFill);
    G4Colour colour22;
    G4VisAttributes attribs22;
    colour22=G4Colour(0.,1.,0.);
    attribs22=G4VisAttributes(colour22);
    circle22.SetVisAttributes(attribs22);
    pVVisManager->Draw(circle22);
    
    G4Square square222(G4Point3D(200.*cm,100.*cm,0.));
    square222.SetWorldDiameter(100.*cm);
    //square222.SetFillStyle(G4Square::filled);
    G4Colour colour222;
    G4VisAttributes attribs222;
    colour222=G4Colour(0.,1.,0.);
    attribs222=G4VisAttributes(colour222);
    square222.SetVisAttributes(attribs222);
    pVVisManager->Draw(square222);
    
    G4Square square33(G4Point3D(-100.*cm,100.*cm,0.));
    square33.SetScreenDiameter(50.*pixels);
    square33.SetFillStyle(G4Square::hashed);
    G4Colour colour33;
    G4VisAttributes attribs33;
    colour33=G4Colour(1.,0.,0.);
    attribs33=G4VisAttributes(colour33);
    square33.SetVisAttributes(attribs33);
    pVVisManager->Draw(square33);

    G4Circle circle333(G4Point3D(-200.*cm,100.*cm,0.));
    circle333.SetWorldDiameter(100.*cm);
    circle333.SetFillStyle(G4Circle::filled);
    G4Colour colour333;
    G4VisAttributes attribs333;
    colour333=G4Colour(1.,0.,0.);
    attribs333=G4VisAttributes(colour333);
    circle333.SetVisAttributes(attribs333);
    pVVisManager->Draw(circle333);
  }
}
