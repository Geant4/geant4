#include "Hits.hh"
#include "G4UnitsTable.hh"
#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"

G4Allocator<Hits> HitAllocator;

//Funzioni che gestiscono le informazioni contenute in ogni hit


Hits::Hits() {}


Hits::~Hits() {}

//overload del costruttore (assegna all'oggetto uno preesistente)
Hits::Hits(const Hits &right) : 
  G4VHit()
{
  trackID   = right.trackID;
  edep      = right.edep;
  pos       = right.pos;
  parentID  = right.parentID;
}


//ridefinizione dell'operatore = (copia oggetti di tipo hit)
const Hits& Hits::operator=(const Hits& right)
{
  trackID   = right.trackID;
  edep      = right.edep;
  pos       = right.pos;
  parentID  = right.parentID;
  return *this;
}


int Hits::operator==(const Hits& right) const
{
  if (trackID != right.trackID || edep != right.edep || 
    pos != right.pos || parentID != right.parentID)
   return 0;
  return 1;
}


void Hits::Draw()
{
  G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
  if(pVVisManager)
  {
    G4Circle circle(pos);
    circle.SetScreenSize(0.04);
    circle.SetFillStyle(G4Circle::filled);
    G4Colour colour(1.,0.,0.);
    G4VisAttributes attribs(colour);
    circle.SetVisAttributes(attribs);
    pVVisManager->Draw(circle);
  }
}



//viene fatto step per step
void Hits::Print()
{
  G4cout << "  trackID: " << trackID << G4endl;
  G4cout << "  energy deposit: " << G4BestUnit(edep,"Energy")
	 << "  position: " << G4BestUnit(pos,"Length") << G4endl;
}


