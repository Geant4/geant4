// Rich advanced example for Geant4
// RichTbHit.cc for Rich of LHCb
// History:
// Created: Sajan Easo (Sajan.Easo@cern.ch)
// Revision and changes: Patricia Mendez (Patricia.Mendez@cern.ch)
/////////////////////////////////////////////////////////////////////////////
#include "RichTbHit.hh"
#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "G4Transform3D.hh"
#include "G4LogicalVolume.hh"
#include "RichTbVisManager.hh"

G4Allocator<RichTbHit> RichTbHitAllocator;

RichTbHit::RichTbHit()
{;}

RichTbHit::~RichTbHit()
{;}

RichTbHit::RichTbHit(const RichTbHit &right)
{
  edep = right.edep;
  posAtSilicon = right.posAtSilicon;
  posAtPhotoCathode=right.posAtPhotoCathode;
  CurHpdNum=right.CurHpdNum;
  CurSectNum=right.CurSectNum;
  CurPixelNum=right.CurPixelNum;
}

const RichTbHit& RichTbHit::operator=(const RichTbHit &right)
{
  edep = right.edep;
  posAtSilicon = right.posAtSilicon;
  posAtPhotoCathode=right.posAtPhotoCathode;
  CurHpdNum=right.CurHpdNum;
  CurSectNum=right.CurSectNum;
  CurPixelNum=right.CurPixelNum;
  
  return *this;
}

int RichTbHit::operator==(const RichTbHit &right) const
{
  return 0;
}

void RichTbHit::Draw()
{

  // The folowing does not work anymore .. SE 26-04-01
  G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
  if(pVVisManager)
   {

  G4Circle circle(posAtPhotoCathode);
  circle.SetScreenSize(0.04);
  circle.SetFillStyle(G4Circle::filled);
  G4Colour colour(1.0,0.0,0.0);
  G4VisAttributes attribs(colour);
  circle.SetVisAttributes(attribs);
  pVVisManager->Draw(circle);
  }
}
void RichTbHit::DrawWithVisM(RichTbVisManager* pVisManager)
{

  G4VisManager* pVVisManager = pVisManager;
  if(pVVisManager)
   {

  G4Circle circle(posAtPhotoCathode);
  circle.SetScreenSize(0.04);
  circle.SetFillStyle(G4Circle::filled);
  G4Colour colour(1.0,0.0,0.0);
  G4VisAttributes attribs(colour);
  circle.SetVisAttributes(attribs);
  pVVisManager->Draw(circle);
  }
}
void RichTbHit::Print()
{;}


// This is a forward declarations of an instantiated G4Allocator<Type> object.
// It has been added in order to make code portable for the GNU g++ 
// (release 2.7.2) compiler. 
// Whenever a new Type is instantiated via G4Allocator, it has to be forward
// declared to make object code (compiled with GNU g++) link successfully. 
// 
#ifdef GNU_GCC
  template class G4Allocator<RichTbHit>;
#endif


