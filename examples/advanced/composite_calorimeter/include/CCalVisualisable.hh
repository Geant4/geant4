//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
///////////////////////////////////////////////////////////////////////////////
// File: CCalVisualisable.hh
// Description: Sets visualisable attributes from the information in flat file
///////////////////////////////////////////////////////////////////////////////
#ifndef CCalVisualisable_hh
#define CCalVisualisable_hh 1

#include "globals.hh"

#if __SUNPRO_CC==0x420
   typedef G4bool bool;
#endif

class CCalVisualisable {
public:

  //Here we define the different type of volumes we consider.
  enum visType {Sensitive=0,
		Electronics=1,
		Support=2,
		Cable=3, 
		Absorber=4,
		OtherServices=5,
		PseudoVolumes=6,
		TotalVisTypes=7,
		Undefined=-1};

private:

  //This class groups the important visualization parameters.
  class visParameters {
  public:
    visParameters(bool v=false, double r=1,double g=1,double b=1, bool w=true):
      visibility(v),rColor(r),gColor(g),bColor(b),wireframe(w) {}
    bool visibility;
    double rColor;
    double gColor;
    double bColor;
    bool wireframe;
  };
  
public:
  //Constructs this object from this file
  CCalVisualisable(G4String file);
  
  virtual ~CCalVisualisable() {}

  //Reads this object from file
  bool readFile(G4String file);

  //Sets visibility to true for Sensitive and to false otherwise.
  void setDefault();

  //Get & Set methods.
  bool isVisible (visType v) const {return theParameters[v].visibility;}
  void setVisible(visType v, bool flag=true){theParameters[v].visibility=flag;}

  double colorRed  (visType v) const {return theParameters[v].rColor;}
  double colorGreen(visType v) const {return theParameters[v].gColor;}
  double colorBlue (visType v) const {return theParameters[v].bColor;}
  void setColorRed  (visType v, double r) {theParameters[v].rColor=r;}
  void setColorGreen(visType v, double g) {theParameters[v].gColor=g;}
  void setColorBlue (visType v, double b) {theParameters[v].bColor=b;}
  void setColor(visType v, double r, double g, double b);

  bool isWireFrame (visType v) const {return theParameters[v].wireframe;}
  void setWireFrame(visType v, bool wf=true){theParameters[v].wireframe=wf;}

  
  
protected:
  //Read this object from visFile
  static void setPath();
  bool readFile();

private:
  static const char* pathName;                //Path in which to look for files
  G4String visFile;                           //File with visualization info
  visParameters theParameters[TotalVisTypes]; //Visualisation parameters

  double checkColorRange(double color, char type) const;
};

#endif



