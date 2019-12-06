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

class CCalVisualisable
{
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
  class visParameters
  {
   public:
     visParameters(G4bool v=false, G4double r=1,G4double g=1,
                   G4double b=1, G4bool w=true)
       : visibility(v),rColor(r),gColor(g),bColor(b),wireframe(w) {}
     G4bool visibility;
     G4double rColor;
     G4double gColor;
     G4double bColor;
     G4bool wireframe;
  };
  
public:
  //Constructs this object from this file
  CCalVisualisable(G4String file);
  
  virtual ~CCalVisualisable() {}

  //Reads this object from file
  G4bool readFile(G4String file);

  //Sets visibility to true for Sensitive and to false otherwise.
  void setDefault();

  //Get & Set methods.
  G4bool isVisible (visType v) const
    {return theParameters[v].visibility;}
  void setVisible(visType v, G4bool flag=true)
    {theParameters[v].visibility=flag;}

  G4double colorRed  (visType v) const {return theParameters[v].rColor;}
  G4double colorGreen(visType v) const {return theParameters[v].gColor;}
  G4double colorBlue (visType v) const {return theParameters[v].bColor;}
  void setColorRed  (visType v, G4double r) {theParameters[v].rColor=r;}
  void setColorGreen(visType v, G4double g) {theParameters[v].gColor=g;}
  void setColorBlue (visType v, G4double b) {theParameters[v].bColor=b;}
  void setColor(visType v, G4double r, G4double g, G4double b);

  G4bool isWireFrame (visType v) const {return theParameters[v].wireframe;}
  void setWireFrame(visType v, G4bool wf=true){theParameters[v].wireframe=wf;}

protected:
  //Read this object from visFile
  static void setPath();
  G4bool readFile();

private:
  static const char* pathName;                //Path in which to look for files
  G4String visFile;                           //File with visualization info
  visParameters theParameters[TotalVisTypes]; //Visualisation parameters

  G4double checkColorRange(G4double color, char type) const;
};

#endif
