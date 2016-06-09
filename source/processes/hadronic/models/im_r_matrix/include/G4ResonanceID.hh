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
#ifndef G4ResonanceID_h
#define G4ResonanceID_h

#include "G4Types.hh"

class G4ResonanceID
{
  public:
    G4bool IsDelta1232(G4int pdg) {return pdg==2224||pdg==2214||pdg==2114||pdg==1114;}
    G4bool IsDelta1600(G4int pdg) {return pdg==31114||pdg==32114||pdg==32214||pdg==32224;}
    G4bool IsDelta1620(G4int pdg) {return pdg==1112||pdg==1212||pdg==2122||pdg==2222;}
    G4bool IsDelta1700(G4int pdg) {return pdg==11114||pdg==12114||pdg==12214||pdg==12224;}
    G4bool IsDelta1900(G4int pdg) {return pdg==11112||pdg==11212||pdg==12122||pdg==12222;}
    G4bool IsDelta1905(G4int pdg) {return pdg==1116||pdg==1216||pdg==2126||pdg==2226;}
    G4bool IsDelta1910(G4int pdg) {return pdg==21112||pdg==21212||pdg==22122||pdg==22222;}
    G4bool IsDelta1920(G4int pdg) {return pdg==21114||pdg==22114||pdg==22214||pdg==22224;}
    G4bool IsDelta1930(G4int pdg) {return pdg==11116||pdg==11216||pdg==12126||pdg==12226;}
    G4bool IsDelta1950(G4int pdg) {return pdg==1118||pdg==2118||pdg==2218||pdg==2228;}
    
    G4bool IsN1440(G4int pdg) {return pdg==12112||pdg==12212;}
    G4bool IsN1520(G4int pdg) {return pdg==1214||pdg==2124;}
    G4bool IsN1535(G4int pdg) {return pdg==22112||pdg==22212;}
    G4bool IsN1650(G4int pdg) {return pdg==32112||pdg==32212;}
    G4bool IsN1675(G4int pdg) {return pdg==2116||pdg==2216;}
    G4bool IsN1680(G4int pdg) {return pdg==12116||pdg==12216;}
    G4bool IsN1700(G4int pdg) {return pdg==22124||pdg==21214;}
    G4bool IsN1710(G4int pdg) {return pdg==42212||pdg==42112;}
    G4bool IsN1720(G4int pdg) {return pdg==32124||pdg==31214;}
    G4bool IsN1900(G4int pdg) {return pdg==42124||pdg==41214;}
    G4bool IsN1990(G4int pdg) {return pdg==12218||pdg==12118;}
    G4bool IsN2090(G4int pdg) {return pdg==52214||pdg==52114;}
    G4bool IsN2190(G4int pdg) {return pdg==2128||pdg==1218;}
    G4bool IsN2220(G4int pdg) {return pdg==100002210||pdg==100002110;}
    G4bool IsN2250(G4int pdg) {return pdg==100012210||pdg==100012110;}
    
    G4bool IsL1405(G4int pdg) {return pdg==13122;}
    G4bool IsL1520(G4int pdg) {return pdg==3124;}
    G4bool IsL1600(G4int pdg) {return pdg==23122;}
    G4bool IsL1670(G4int pdg) {return pdg==33122;}
    G4bool IsL1690(G4int pdg) {return pdg==13124;}
    G4bool IsL1800(G4int pdg) {return pdg==43122;}
    G4bool IsL1810(G4int pdg) {return pdg==53122;}
    G4bool IsL1820(G4int pdg) {return pdg==3126;}
    G4bool IsL1830(G4int pdg) {return pdg==13126;}
    G4bool IsL1890(G4int pdg) {return pdg==23124;}
    G4bool IsL2100(G4int pdg) {return pdg==3128;}
    G4bool IsL2110(G4int pdg) {return pdg==23126;}
    
    G4bool IsS1192(G4int pdg) {return pdg==3222||pdg==3212||pdg==3112;}
    G4bool IsS1385(G4int pdg) {return pdg==3214||pdg==3114||pdg==3224;}
    G4bool IsS1660(G4int pdg) {return pdg==13212||pdg==13112||pdg==13222;}
    G4bool IsS1670(G4int pdg) {return pdg==13214||pdg==13114||pdg==13224;}
    G4bool IsS1750(G4int pdg) {return pdg==23212||pdg==23112||pdg==23222;}
    G4bool IsS1775(G4int pdg) {return pdg==3216||pdg==3116||pdg==3226;}
    G4bool IsS1915(G4int pdg) {return pdg==13216||pdg==13116||pdg==13226;}
    G4bool IsS1940(G4int pdg) {return pdg==23214||pdg==23114||pdg==23224;}
    G4bool IsS2030(G4int pdg) {return pdg==3218||pdg==3118||pdg==3228;}
    
    G4bool IsX1530(G4int pdg) {return pdg==3324||pdg==3314;}
    G4bool IsX1690(G4int pdg) {return pdg==23324||pdg==23314;}
    G4bool IsX1820(G4int pdg) {return pdg==13324||pdg==13314;}
    G4bool IsX1950(G4int pdg) {return pdg==33324||pdg==33314;}
    G4bool IsX2030(G4int pdg) {return pdg==13326||pdg==13316;}
};

#endif
