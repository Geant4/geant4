// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: MySteppingAction.cc,v 1.2 1999-04-28 14:09:19 johna Exp $
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

void MySteppingAction::UserSteppingAction(const G4Step* pStep) {
#ifdef G4VIS_USE_OPACS
  OHistogram h1 = MyRunAction::get_1d();
  OHistogram h2 = MyRunAction::Get2d();

  HepRandom::setTheEngine (&theJamesEngine);
  double     james = RandGauss::shoot(0.3,0.1);
  OHistogramFillOneDimensional(h1,james,0.01);

  HepRandom::setTheEngine(&theDRand48Engine);
  double     d48 = RandGauss::shoot(0.7,0.1);
  OHistogramFillTwoDimensional(h2,james,d48,0.01);   
#endif

  // User Action Example - begin snippet.
  G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
  static int coutCount = 0;
  if (coutCount < 10) {
    coutCount++;
    G4cout << "MySteppingAction::UserSteppingAction called." << endl;
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
	   <<  pl << endl;
    }
    pVVisManager -> Draw (pl);
  }
  // User Action Example - end snippet.
}




