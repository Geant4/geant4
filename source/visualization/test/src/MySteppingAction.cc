// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: MySteppingAction.cc,v 1.7 2000-10-18 13:52:20 allison Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// 16/Apr/1997  J. Allison:  For visualization/test/test19.

#include "MySteppingAction.hh"

#ifdef G4VIS_USE_OPACS
#include "MyRunAction.hh"
#endif

#include "G4SteppingManager.hh"
#include "G4ParticleDefinition.hh"
#include "G4VVisManager.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4Polyline.hh"
#include "G4Circle.hh"
#include "G4Square.hh"

void MySteppingAction::UserSteppingAction(const G4Step* pStep) {
#ifdef G4VIS_USE_OPACS
  OHistogram h1 = MyRunAction::get_1d();
  OHistogram h2 = MyRunAction::Get2d();

  HepRandom::setTheEngine (&theJamesEngine);
  double     james = G4RandGauss::shoot(0.3,0.1);
  OHistogramFillOneDimensional(h1,james,0.01);

  HepRandom::setTheEngine(&theDRand48Engine);
  double     d48 = G4RandGauss::shoot(0.7,0.1);
  OHistogramFillTwoDimensional(h2,james,d48,0.01);   
#endif

  // User Action Example - begin snippet.
  G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
  static int coutCount = 0;
  if (coutCount < 10) {
    coutCount++;
    G4cout << "MySteppingAction::UserSteppingAction called." << G4endl;
  }
  if(pVVisManager) {
    G4Polyline pl;
    G4Colour c;
    G4double chg =
      pStep -> GetTrack () -> GetDefinition () -> GetPDGCharge ();
    if (chg != 0.) c = G4Colour (1., 0., 0.);
    else           c = G4Colour (0., 1., 0.);
    G4VisAttributes va (c);
    pl.SetVisAttributes (&va);
    pl.append (pStep -> GetPreStepPoint () -> GetPosition ());
    pl.append (pStep -> GetPostStepPoint () -> GetPosition ());
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

    G4Circle circle22(G4Point3D(100.*cm,100.*cm,0.));
    circle22.SetScreenDiameter(50.*cm);
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
    square33.SetScreenDiameter(50.*cm);
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
  // User Action Example - end snippet.
}

