///////////////////////////////////////////////////////////////////////////////
// File: HcalTB96.hh
// Date: 08/00 
// Modifications: 
// Description: Equipped to construct the geometry of the 96 Test Beam
///////////////////////////////////////////////////////////////////////////////
#ifndef HcalTB96_h
#define HcalTB96_h 1

#include "CMSDetector.hh"

class HcalTB96: public CMSDetector {
public:
  //Constructor and Destructor
  HcalTB96(const G4String &name);
  virtual ~HcalTB96();

  //Get Methods
  G4String getMaterial()                  const {return genMaterial;}
  double   getDy_2Hall()                  const {return dy_2Hall;}
  double   getDx_2Hall()                  const {return dx_2Hall;}

protected:
  virtual int readFile();
  virtual void constructDaughters();

private:
  G4String genMaterial;            //General material
  double   dy_2Hall;               //Half width     of the Experimental Hall
  double   dx_2Hall;               //Half thickness of the Experimental Hall
};

#endif
