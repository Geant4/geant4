///////////////////////////////////////////////////////////////////////////////
// File: G4Able.cc
// Date: 04/98 I. Gonzalez
// Modified: 11/98 I.G. => Added Sensitivity features
//           02/99 I.G. => Added Output for visibility features
//           03/00 S.B. => Directory reorganisation
///////////////////////////////////////////////////////////////////////////////
//Comment/Uncomment next line to unset/set debug information printing
//#define debug
//#define ddebug

#ifdef ddebug
  #include "G4Timer.hh"
#endif
#ifdef debug
  #include "utils.hh"
#endif

#include "GeometryConfiguration.hh"
#include "G4Able.hh"
#include "G4GeometryConfiguration.hh"

#include "G4Color.hh"
#include "G4VisAttributes.hh"

G4Able::G4Able(G4String name):
  detPhysicalVolume(0),
  visProperties(G4GeometryConfiguration::getInstance()->getFileName(name)+".vis"),
  sensitivity(false),
  g4ableName(name){
    //Initialize g4VisAtt pointers
    for (int i=0; i<Visualisable::TotalVisTypes; i++) {
      g4VisAtt[i]=0;
    }
    sensitivity = 
      G4GeometryConfiguration::getInstance()->getSensitiveFlag(name);
}

G4VPhysicalVolume* G4Able::PhysicalVolume(G4VPhysicalVolume* pv) {
  //If detPhysicalVolume is not (nil) the volume has already been built
  //so return it. In other case, construct it and its daughters, then
  //check for sensitivity and build it if set.
#ifdef ddebug
  G4Timer timer;
  timer.Start();
#endif
  if (GeometryConfiguration::getInstance()->getConstructFlag(G4Name())!=0) {
    if (!detPhysicalVolume) {
      detPhysicalVolume = constructIn(pv);
      for(int i = 0; i < theG4DetectorsInside.size(); i++) {
	theG4DetectorsInside[i]->PhysicalVolume(detPhysicalVolume);
      }
      if (sensitivity) {
#ifdef debug
	G4cout << "==> Making " << detPhysicalVolume->GetName() << " sensitive..." 
	       << endl;
#endif
	constructSensitive();
      } //if sensitivity
    } //if sensitive
  } //if construct
  else {
    G4cout << "NOTE: You decided to skip the construction of " 
	   << G4Name() << endl;
  }
#ifdef ddebug
  timer.Stop();
  G4cout << tab << "G4Able::PhysicalVolume(...) --> time spent: " 
	 << timer << endl;
#endif
  return detPhysicalVolume;
}

void G4Able::AddG4Able(G4Able* det) {
  theG4DetectorsInside.push_back(det);
}

void G4Able::setVisType(Visualisable::visType vt, G4LogicalVolume* log) {
  if (!g4VisAtt[vt]) {
#ifdef debug
    G4cout << "G4Able::setVisType: Constructing G4VisAttributes for " 
	   << log->GetName() << " as " << vt << endl;
#endif
    G4Color col(visProperties.colorRed(vt),
		visProperties.colorGreen(vt),
		visProperties.colorBlue(vt));
    G4bool wf      = visProperties.isWireFrame(vt);
    G4bool visible = visProperties.isVisible(vt);

#ifdef debug
    G4cout << "Color: " 
	   << visProperties.colorRed(vt)   << ", " 
	   << visProperties.colorGreen(vt) << ", "
	   << visProperties.colorBlue(vt)  << tab
	   << "Wireframe: " << wf << tab
	   << "Visible: " << visible << endl;
#endif
    g4VisAtt[vt] = new G4VisAttributes(col);
    g4VisAtt[vt]->SetForceWireframe(wf);
    g4VisAtt[vt]->SetVisibility(visible);
  }
  log->SetVisAttributes(g4VisAtt[vt]);
}



G4bool G4Able::operator==(const G4Able& right) const {
  return detPhysicalVolume==right.detPhysicalVolume;
}



//========================================================================
//Protected and private methods.

//========================================================================
//Global operators
ostream& operator<<(ostream& os, const G4Able& det) {
  if (det.detPhysicalVolume)
    os << "Physical volume already constructed." << endl;
  else
    os << "Physical volume still not constructed." << endl;

  if (det.isSensitive())
    os << "and it is Sensitive" << endl;
  else
    os << "and it is not Sensitive" << endl;
  
  return os;
}
