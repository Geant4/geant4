#include "G4VParallelStepper.hh"
#include "G4ImportanceFinder.hh"
#include "G4VIStore.hh"
#include "g4std/strstream"
#include "G4PStepStream.hh"

G4double G4ImportanceFinder::
GetIPre_over_IPost(const G4PTouchableKey &prekey,
		   const G4PTouchableKey &postkey) const{

  
  G4double  ipre = fIStore.GetImportance(prekey);
  G4double ipost = fIStore.GetImportance(postkey);

  if (ipre <= 0 || ipost <=0 ) {
    ostrstream os;
    os << "ipre <= 0 || ipost <=0, preTouchableKey = " << prekey 
       << ", postTouchableKey = " << postkey << '\0';
    Error(os.str());
  }
  return ipre/ipost;
}
