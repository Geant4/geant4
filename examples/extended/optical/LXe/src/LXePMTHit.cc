#include "LXePMTHit.hh"
#include "G4ios.hh"
#include "G4VVisManager.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"

G4Allocator<LXePMTHit> LXePMTHitAllocator;

//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
LXePMTHit::LXePMTHit()
  :pmtNumber(-1),photons(0),physVol(0),drawit(false)
{}

//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
LXePMTHit::~LXePMTHit()
{}

//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
LXePMTHit::LXePMTHit(const LXePMTHit &right)
  : G4VHit()
{
  pmtNumber=right.pmtNumber;
  photons=right.photons;
  physVol=right.physVol;
  drawit=right.drawit;
}

//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
const LXePMTHit& LXePMTHit::operator=(const LXePMTHit &right){
  pmtNumber = right.pmtNumber;
  photons=right.photons;
  physVol=right.physVol;
  drawit=right.drawit;
  return *this;
}

//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
G4int LXePMTHit::operator==(const LXePMTHit &right) const{
  return (pmtNumber==right.pmtNumber);
}

//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
void LXePMTHit::Draw(){
  if(drawit&&physVol){ //ReDraw only the PMTs that have hit counts > 0
    //Also need a physical volume to be able to draw anything
    G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
    if(pVVisManager){//Make sure that the VisManager exists
      G4VisAttributes attribs(G4Colour(1.,0.,0.));
      attribs.SetForceSolid(true);
      G4RotationMatrix rot;
      if(physVol->GetRotation())//If a rotation is defined use it
	rot=*(physVol->GetRotation());
      G4Transform3D trans(rot,physVol->GetTranslation());//Create transform
      pVVisManager->Draw(*physVol,attribs,trans);//Draw it
    }
  }
}

//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
void LXePMTHit::Print(){
}









