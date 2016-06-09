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
// This is the *BASIC* version of Hadrontherapy, a Geant4-based application
// See more at: http://g4advancedexamples.lngs.infn.it/Examples/hadrontherapy
//
// Visit the Hadrontherapy web site (http://www.lns.infn.it/link/Hadrontherapy) to request 
// the *COMPLETE* version of this program, together with its documentation;
// Hadrontherapy (both basic and full version) are supported by the Italian INFN
// Institute in the framework of the MC-INFN Group
//

#ifndef HadrontherapyModulator_H
#define HadrontherapyModulator_H 1

class G4Tubs;
class G4LogicalVolume;
class G4VPhysicalVolume;


class HadrontherapyModulator
{
public:

  HadrontherapyModulator();
  ~HadrontherapyModulator();

  void BuildModulator(G4VPhysicalVolume*);  
  void SetModulatorAngle(G4double);

private:

  G4RotationMatrix* rm;
  G4VPhysicalVolume* physiMotherMod;

  G4Tubs*            solidMod0;    
  G4LogicalVolume*   logicMod0;   
  G4VPhysicalVolume* physiMod0;
  
  
  G4Tubs*            solidMod1;   
  G4LogicalVolume*   logicMod1;   
  G4VPhysicalVolume* physiMod1;
  
  G4Tubs*            solidMod2;   
  G4LogicalVolume*   logicMod2;   
  G4VPhysicalVolume* physiMod2;
  
  G4Tubs*            solidMod3;   
  G4LogicalVolume*   logicMod3;   
  G4VPhysicalVolume* physiMod3;
  
  G4Tubs*            solidMod4;  
  G4LogicalVolume*   logicMod4;  
  G4VPhysicalVolume* physiMod4;

  G4Tubs*            solidMod5;  
  G4LogicalVolume*   logicMod5;   
  G4VPhysicalVolume* physiMod5;
  
  G4Tubs*            solidMod6;  
  G4LogicalVolume*   logicMod6;  
  G4VPhysicalVolume* physiMod6;
  
  G4Tubs*            solidMod7;   
  G4LogicalVolume*   logicMod7;   
  G4VPhysicalVolume* physiMod7;
  
  G4Tubs*            solidMod8;  
  G4LogicalVolume*   logicMod8;   
  G4VPhysicalVolume* physiMod8;
  
  G4Tubs*            solidMod9;   
  G4LogicalVolume*   logicMod9;   
  G4VPhysicalVolume* physiMod9;
  
  G4Tubs*            solidMod10;   
  G4LogicalVolume*   logicMod10;  
  G4VPhysicalVolume* physiMod10;
  
  G4Tubs*            solidMod11;   
  G4LogicalVolume*   logicMod11; 
  G4VPhysicalVolume* physiMod11;
  
  G4Tubs*            solidMod12;   
  G4LogicalVolume*   logicMod12;  
  G4VPhysicalVolume* physiMod12;
  
  G4Tubs*            solidMod13;  
  G4LogicalVolume*   logicMod13;   
  G4VPhysicalVolume* physiMod13;
  
  G4Tubs*            solidMod14;   
  G4LogicalVolume*   logicMod14;  
  G4VPhysicalVolume* physiMod14;
  
  G4Tubs*            solidMod15;  
  G4LogicalVolume*   logicMod15;  
  G4VPhysicalVolume* physiMod15;
  
  G4Tubs*            solidMod16;   
  G4LogicalVolume*   logicMod16;   
  G4VPhysicalVolume* physiMod16;
  
  G4Tubs*            solidMod17;  
  G4LogicalVolume*   logicMod17; 
  G4VPhysicalVolume* physiMod17;
  
  G4Tubs*            solidMod18;   
  G4LogicalVolume*   logicMod18;  
  G4VPhysicalVolume* physiMod18;
  
  G4Tubs*            solidMod20;  
  G4LogicalVolume*   logicMod20;  
  G4VPhysicalVolume* physiMod20;

  G4Tubs*            solidMod21;  
  G4LogicalVolume*   logicMod21;   
  G4VPhysicalVolume* physiMod21;

  G4Tubs*            solidMod22;    
  G4LogicalVolume*   logicMod22;  
  G4VPhysicalVolume* physiMod22;
  
  G4Tubs*            solidMod23;   
  G4LogicalVolume*   logicMod23; 
  G4VPhysicalVolume* physiMod23;
  
  G4Tubs*            solidMod24;  
  G4LogicalVolume*   logicMod24;  
  G4VPhysicalVolume* physiMod24;

  G4Tubs*            solidMod25;   
  G4LogicalVolume*   logicMod25;  
  G4VPhysicalVolume* physiMod25;

  G4Tubs*            solidMod26;    
  G4LogicalVolume*   logicMod26;  
  G4VPhysicalVolume* physiMod26;
  
  G4Tubs*            solidMod27;  
  G4LogicalVolume*   logicMod27;  
  G4VPhysicalVolume* physiMod27;
  
  G4Tubs*            solidMod28;  
  G4LogicalVolume*   logicMod28;  
  G4VPhysicalVolume* physiMod28;
  
  G4Tubs*            solidMod29;   
  G4LogicalVolume*   logicMod29;  
  G4VPhysicalVolume* physiMod29;
  
  G4Tubs*            solidMod30;  
  G4LogicalVolume*   logicMod30;  
  G4VPhysicalVolume* physiMod30;
  
  G4Tubs*            solidMod31;   
  G4LogicalVolume*   logicMod31;  
  G4VPhysicalVolume* physiMod31;
  
  G4Tubs*            solidMod32;  
  G4LogicalVolume*   logicMod32;   
  G4VPhysicalVolume* physiMod32;
  
  G4Tubs*            solidMod33;  
  G4LogicalVolume*   logicMod33;   
  G4VPhysicalVolume* physiMod33;
  
  G4Tubs*            solidMod34;   
  G4LogicalVolume*   logicMod34;  
  G4VPhysicalVolume* physiMod34;
  
  G4Tubs*            solidMod35;  
  G4LogicalVolume*   logicMod35;  
  G4VPhysicalVolume* physiMod35;

  G4Tubs*            solidMod36;  
  G4LogicalVolume*   logicMod36;  
  G4VPhysicalVolume* physiMod36;
  
  G4Tubs*            solidMod37;   
  G4LogicalVolume*   logicMod37;   
  G4VPhysicalVolume* physiMod37;

  G4Tubs*            solidMod38;  
  G4LogicalVolume*   logicMod38;   
  G4VPhysicalVolume* physiMod38;
  
  
  G4Tubs*            solidMod40;   // pointer to the 
  G4LogicalVolume*   logicMod40;  
  G4VPhysicalVolume* physiMod40;
  
  G4Tubs*            solidMod41;    
  G4LogicalVolume*   logicMod41;  
  G4VPhysicalVolume* physiMod41;
  
  G4Tubs*            solidMod42;  
  G4LogicalVolume*   logicMod42;  
  G4VPhysicalVolume* physiMod42;
  
  G4Tubs*            solidMod43;   
  G4LogicalVolume*   logicMod43;  
  G4VPhysicalVolume* physiMod43;
  
  G4Tubs*            solidMod44;  
  G4LogicalVolume*   logicMod44;   
  G4VPhysicalVolume* physiMod44;
  
  G4Tubs*            solidMod45;    
  G4LogicalVolume*   logicMod45;   
  G4VPhysicalVolume* physiMod45;

  G4Tubs*            solidMod46;   
  G4LogicalVolume*   logicMod46;   
  G4VPhysicalVolume* physiMod46;
  
  G4Tubs*            solidMod47;  
  G4LogicalVolume*   logicMod47; 
  G4VPhysicalVolume* physiMod47;
  
  G4Tubs*            solidMod48;   
  G4LogicalVolume*   logicMod48;  
  G4VPhysicalVolume* physiMod48;
  
  G4Tubs*            solidMod49;   
  G4LogicalVolume*   logicMod49;  
  G4VPhysicalVolume* physiMod49;
  
  G4Tubs*            solidMod50;   
  G4LogicalVolume*   logicMod50; 
  G4VPhysicalVolume* physiMod50;
  
  G4Tubs*            solidMod51; 
  G4LogicalVolume*   logicMod51;  
  G4VPhysicalVolume* physiMod51;

  G4Tubs*            solidMod52;  
  G4LogicalVolume*   logicMod52;   
  G4VPhysicalVolume* physiMod52;
  
  G4Tubs*            solidMod53;  
  G4LogicalVolume*   logicMod53;  
  G4VPhysicalVolume* physiMod53;
  
  G4Tubs*            solidMod54;  
  G4LogicalVolume*   logicMod54;   
  G4VPhysicalVolume* physiMod54;

  G4Tubs*            solidMod55;    
  G4LogicalVolume*   logicMod55;   
  G4VPhysicalVolume* physiMod55;
  
  G4Tubs*            solidMod56;   
  G4LogicalVolume*   logicMod56;  
  G4VPhysicalVolume* physiMod56;

  G4Tubs*            solidMod57;   
  G4LogicalVolume*   logicMod57;  
  G4VPhysicalVolume* physiMod57;
  
  G4Tubs*            solidMod58;    
  G4LogicalVolume*   logicMod58;   
  G4VPhysicalVolume* physiMod58;

  G4Tubs*            solidMod60;   
  G4LogicalVolume*   logicMod60;   
  G4VPhysicalVolume* physiMod60;
  
  G4Tubs*            solidMod61;   
  G4LogicalVolume*   logicMod61;   
  G4VPhysicalVolume* physiMod61;
  
  G4Tubs*            solidMod62;  
  G4LogicalVolume*   logicMod62;   
  G4VPhysicalVolume* physiMod62;
    
  G4Tubs*            solidMod63;  
  G4LogicalVolume*   logicMod63;   
  G4VPhysicalVolume* physiMod63;

  G4Tubs*            solidMod64;  
  G4LogicalVolume*   logicMod64;  
  G4VPhysicalVolume* physiMod64;
  
  G4Tubs*            solidMod65;  
  G4LogicalVolume*   logicMod65;  
  G4VPhysicalVolume* physiMod65;
  
  G4Tubs*            solidMod66;  
  G4LogicalVolume*   logicMod66; 
  G4VPhysicalVolume* physiMod66;
  
  G4Tubs*            solidMod67;  
  G4LogicalVolume*   logicMod67; 
  G4VPhysicalVolume* physiMod67;
  
  G4Tubs*            solidMod68;  
  G4LogicalVolume*   logicMod68;   
  G4VPhysicalVolume* physiMod68;
  
  G4Tubs*            solidMod69;   
  G4LogicalVolume*   logicMod69;  
  G4VPhysicalVolume* physiMod69;
  
  G4Tubs*            solidMod70;  
  G4LogicalVolume*   logicMod70;   
  G4VPhysicalVolume* physiMod70;
  
  G4Tubs*            solidMod71;  
  G4LogicalVolume*   logicMod71;   
  G4VPhysicalVolume* physiMod71;
  
  G4Tubs*            solidMod72;  
  G4LogicalVolume*   logicMod72;  
  G4VPhysicalVolume* physiMod72;
    
  G4Tubs*            solidMod73;   
  G4LogicalVolume*   logicMod73;  
  G4VPhysicalVolume* physiMod73;
  
  G4Tubs*            solidMod74;   
  G4LogicalVolume*   logicMod74;   
  G4VPhysicalVolume* physiMod74;
  
  G4Tubs*            solidMod75; 
  G4LogicalVolume*   logicMod75; 
  G4VPhysicalVolume* physiMod75;
  
  G4Tubs*            solidMod76;    
  G4LogicalVolume*   logicMod76;  
  G4VPhysicalVolume* physiMod76;

  G4Tubs*            solidMod77;  
  G4LogicalVolume*   logicMod77;  
  G4VPhysicalVolume* physiMod77;
  
  G4Tubs*            solidMod78;   
  G4LogicalVolume*   logicMod78;  
  G4VPhysicalVolume* physiMod78;
   
};
#endif
