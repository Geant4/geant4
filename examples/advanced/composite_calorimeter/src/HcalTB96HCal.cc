///////////////////////////////////////////////////////////////////////////////
// File: HcalTB96HCal.cc
// Date: 08/00 S.Banerjee
// Modifications:
///////////////////////////////////////////////////////////////////////////////
#include "HcalTB96HCal.hh"

#include <fstream.h>
#include "utils.hh"

//#define debug

HcalTB96HCal::HcalTB96HCal(const G4String &name):
  CMSDetector(name),xposBox(0),typeLayerScnt(0),mothLayerScnt(0),
  xposLayerScnt(0),typeLayerAbs(0),mothLayerAbs(0),xposLayerAbs(0),
  dx_2Absorber(0),dy_2ScntLayer(0),dx_2ScntLayer(0),dx_2Wrapper(0),
  dx_2FrontPlastic(0),dx_2BackPlastic(0),dx_2Scintillator(0) {}

HcalTB96HCal::~HcalTB96HCal() {
  if (xposBox)
    delete[] xposBox;
  if (typeLayerScnt)
    delete[] typeLayerScnt;
  if (mothLayerScnt)
    delete[] mothLayerScnt;
  if (xposLayerScnt)
    delete[] xposLayerScnt;
  if (typeLayerAbs)
    delete[] typeLayerAbs;
  if (mothLayerAbs)
    delete[] mothLayerAbs;
  if (xposLayerAbs)
    delete[] xposLayerAbs;
  if (dx_2Absorber)
    delete[] dx_2Absorber;
  if (dy_2ScntLayer)
    delete[] dy_2ScntLayer;
  if (dx_2ScntLayer)
    delete[] dx_2ScntLayer;
  if (dx_2Wrapper)
    delete[] dx_2Wrapper;
  if (dx_2FrontPlastic)
    delete[] dx_2FrontPlastic;
  if (dx_2BackPlastic)
    delete[] dx_2BackPlastic;
  if (dx_2Scintillator)
    delete[] dx_2Scintillator;
}

int HcalTB96HCal::readFile() {
  ///////////////////////////////////////////////////////////////
  //Let's open the file
  cout << " ==> Opening file " << File() << " to read elements..."
       << endl;

  ifstream is;
  bool ok = openGeomFile(is, pathName, File());
  if (!ok)
    return 0;

  // Find *DO HCal
  findDO(is, G4String("HCal"));

  // Calorimeter boundaries
  readName(is,genMaterial);
  is >> dy_2Cal >> dx_2Cal >> xposCal >> jump;
#ifdef debug
  cout << tab << "General material: " << genMaterial  << " Size " << dy_2Cal
       << ", " << dx_2Cal << " Position " << xposCal << endl;
#endif

  // Boxes
  readName(is,boxMaterial);
  is >> nBox >> dy_2Box >> dx_2Box >> wallThickBox;
  int i = 0;
  xposBox = new double[nBox];
  for (i=0; i<nBox; i++)
    is >> xposBox[i];
#ifdef debug
  cout << tab << "Box material: " << boxMaterial  << " Size " << dy_2Box
       << ", " << dx_2Box << " Wall Thickness " << wallThickBox << " number "
       << nBox << " position ";
  for (i=0; i<nBox; i++)
    cout << i << " " << xposBox[i] << " ";
  cout << endl;
#endif

  // Layers of scintillators
  G4String rubbish;
  readName(is,rubbish);
  is >> nLayerScnt;
  typeLayerScnt = new int[nLayerScnt];
  mothLayerScnt = new int[nLayerScnt];
  xposLayerScnt = new double[nLayerScnt];
  for (i=0; i<nLayerScnt; i++)
    is >> typeLayerScnt[i] >> mothLayerScnt[i] >> xposLayerScnt[i];
#ifdef debug
  cout << tab << nLayerScnt << " Layers of scintillators of type/mother box/"
       << "position" << endl;
  for (i=0; i<nLayerScnt; i++)
    cout << tab << i << " " << typeLayerScnt[i] << " " << mothLayerScnt[i] 
	 << " " << xposLayerScnt[i] << endl;
#endif

  // Layers of absorbers
  readName(is,rubbish);
  is >> nLayerAbs;
  typeLayerAbs = new int[nLayerAbs];
  mothLayerAbs = new int[nLayerAbs];
  xposLayerAbs = new double[nLayerAbs];
  for (i=0; i<nLayerAbs; i++)
    is >> typeLayerAbs[i] >> mothLayerAbs[i] >> xposLayerAbs[i];
#ifdef debug
  cout << tab << nLayerAbs << " Layers of absorbers of type/mother box/"
       << "position" << endl;
  for (i=0; i<nLayerAbs; i++)
    cout << tab << i << " " << typeLayerAbs[i] << " " << mothLayerAbs[i] 
	 << " " << xposLayerAbs[i] << endl;
#endif

  // Absorber parameters
  readName(is,absMaterial);
  is >> nAbsorber >> dy_2Absorber;
  dx_2Absorber = new double[nAbsorber];
  for (i=0; i<nAbsorber; i++)
    is >> dx_2Absorber[i];
#ifdef debug
  cout << "\tAbsorber mad of " << absMaterial << " with " << nAbsorber
       << " types and size " << dy_2Absorber;
  for (i=0;  i<nAbsorber; i++)
    cout << "  " << i << " " << dx_2Absorber[i];
  cout << endl;
#endif

  // Scintillator parameters
  readName(is,scntMaterial);
  readName(is,wrapMaterial);
  readName(is,plasMaterial);
  is >> nScintillator;
  dy_2ScntLayer    = new double[nScintillator];
  dx_2ScntLayer    = new double[nScintillator];
  dx_2Wrapper      = new double[nScintillator];
  dx_2FrontPlastic = new double[nScintillator];
  dx_2BackPlastic  = new double[nScintillator];
  dx_2Scintillator = new double[nScintillator];
  for (i=0; i<nScintillator; i++)
    is >> dy_2ScntLayer[i] >> dx_2ScntLayer[i] >> dx_2Wrapper[i]
       >> dx_2FrontPlastic[i] >> dx_2BackPlastic[i] >> dx_2Scintillator[i];
#ifdef debug
  cout << tab << nScintillator << " Scintillator layers made of " 
       << scntMaterial << " " << wrapMaterial << " and " << plasMaterial 
       << " of sizes " << endl;
  for (i=0; i<nScintillator; i++)
    cout << tab << i << " " << dy_2ScntLayer[i] << " " << dx_2ScntLayer[i] 
	 << " " << dx_2Wrapper[i] << " " << dx_2FrontPlastic[i] << " " 
	 << dx_2BackPlastic[i] << " " << dx_2Scintillator[i] << endl;
#endif
  
  ///////////////////////////////////////////////////////////////
  // Close the file
  cout << " ==> Closing file " << File() << endl;
  is.close();

  return 1;

}

void HcalTB96HCal::constructDaughters() {
}
