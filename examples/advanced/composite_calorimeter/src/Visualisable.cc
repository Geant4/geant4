///////////////////////////////////////////////////////////////////////////////
// File: Visualisable.cc
// Date: 03/98 I. Gonzalez
// Modification: 27/03/00 SB In OSCAR
///////////////////////////////////////////////////////////////////////////////
#include "Visualisable.hh"

#include <fstream.h>
#include "utils.hh"

//Comment/Uncomment next line to unset/set debug information printing
//#define debug

const char* visEnvName = "OSCARVISPATH";
const char* Visualisable::pathName=0;

Visualisable::Visualisable(G4String file):visFile(file) {
  if (!pathName)
    setPath();
  setDefault();
  readFile();
}

bool Visualisable::readFile(G4String file) { 
  visFile=file;
  return readFile();
}

void Visualisable::setDefault(){
  theParameters[Sensitive] = visParameters(true);
  theParameters[Electronics] = visParameters(false);
  theParameters[Support] = visParameters(false);
  theParameters[Cable] = visParameters(false);
  theParameters[OtherServices] = visParameters(false);
  theParameters[PseudoVolumes] = visParameters(false);
}

void Visualisable::setColor(visType v, double r, double g, double b){
  theParameters[v].rColor=r;
  theParameters[v].gColor=g;
  theParameters[v].bColor=b;
}

void Visualisable::setPath() {
  pathName = getenv(visEnvName);
  if (!pathName) {
    cerr << "ERROR: " << visEnvName << " environmental variable not set!" 
	 << endl;
    cerr << "       Set it and restart." << endl;
    exit(-333);
  }
}

bool Visualisable::readFile() {
  if (visFile=="") {
    cerr << "ERROR: No file was specified from which to read Visualisation parameters" 
	 << endl;
    return false;
  }

  ///////////////////////////////////////////////////////////////
  //Let's open the file
  cout << " ==> Opening file " << visFile 
       << " to read visualisation parameters..."
       << endl;

  G4String pathname(pathName);
  ifstream is;
#ifdef debug
  cout << "Viualisable : Path " << pathname << " FIle " << visFile << endl;
#endif
  bool ok = openGeomFile(is, pathname, visFile);
  if (!ok) {
    cout << "WARNING: Could not read " << visFile << endl;
    cout << "         Default visualization parameters will be used." << endl;
    return false;
  }
  
  while (is) {
    G4String viewvol;
    readName(is,viewvol);

    visType vt = Undefined;
    if (viewvol=="Sensitive")
      vt=Sensitive;
    else if (viewvol=="Electronics")
      vt=Electronics;
    else if (viewvol=="Support")
      vt=Support;
    else if (viewvol=="Absorber")
      vt=Absorber;
    else if (viewvol=="Cable")
      vt=Cable;
    else if (viewvol=="OtherServices")
      vt=OtherServices;
    else if (viewvol=="PseudoVolumes")
      vt=PseudoVolumes;
    else {
      vt=Undefined;
      cerr << "WARNING: Unknown type of visualisable object \"" << viewvol
	   << "\"." << endl;
    }


    int isvisible, wireframe;
    double r, g, b;
    
    is >> isvisible >> r >> g >> b >> wireframe >> jump;

    r=checkColorRange(r,'R');
    g=checkColorRange(g,'G');
    b=checkColorRange(b,'B');

    if (vt!=Undefined) {
#ifdef debug
      cout << tab << viewvol << tab << isvisible << tab 
	   << r << " " << g << " "<< b << tab
	   << wireframe << endl;
#endif
      theParameters[vt]=visParameters(isvisible, r, g, b, wireframe);
    }

  }

  ///////////////////////////////////////////////////////////////
  // Close the file
  cout << " ==> Closing file " << visFile << endl;
  is.close();

  return true;
}

double Visualisable::checkColorRange(double cvalue, char ctype) const {
  if (cvalue>1) {
    cerr << "ERROR: In " << visFile << ". Color " << ctype << "=" 
	 << cvalue << " > 1" << endl;
    cerr << "       It will be reset to 1." << endl;
    return 1.;
  }
  if (cvalue<0) {
    cerr << "ERROR: In " << visFile << ". Color " << ctype << "=" 
	 << cvalue << " < 0" << endl;
    cerr << "       It will be reset to 0." << endl;
    return 0.;
  }
  return cvalue;
}
