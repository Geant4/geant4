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
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4HtmlPPReporter.cc 67971 2013-03-13 10:13:24Z gcosmo $
//
// 
// ---------------------------------------------------------------
#include "G4HtmlPPReporter.hh"
#include "G4ios.hh"
#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleDefinition.hh"  
#include "G4ParticleTable.hh"  
#include "G4DecayTable.hh"  
#include "G4VDecayChannel.hh"  
#include "G4Tokenizer.hh"
#include <iomanip>

 G4HtmlPPReporter::G4HtmlPPReporter():G4VParticlePropertyReporter()
{
 
}

 G4HtmlPPReporter::~G4HtmlPPReporter()
{
}    

 void G4HtmlPPReporter::Print(const G4String& option)
{
  SparseOption( option );

  GenerateIndex();

  for (size_t i=0; i< pList.size(); i++){
    G4ParticleDefinition* particle  = G4ParticleTable::GetParticleTable()->FindParticle( pList[i]->GetParticleName() ); 
    GeneratePropertyTable(particle);
  }
}    


void G4HtmlPPReporter::SparseOption(const G4String& option)
{
  G4Tokenizer savedToken( option );
  
  // 1st option : base directory
  baseDir = savedToken();
  if (!baseDir.isNull()) {
    if(baseDir(baseDir.length()-1)!='/') {
      baseDir += "/";
    }
  }
  comment =  savedToken();
}

 void G4HtmlPPReporter::GenerateIndex()
{
  //--- open index file -----
  G4String fileName = baseDir + "index.html";
  std::ofstream outFile(fileName, std::ios::out );
  outFile.setf( std::ios:: scientific, std::ios::floatfield );
  
  // header
  PrintHeader(outFile);
  
  // comment
  outFile << "<! -- " << comment << " -- !> " << G4endl;
  outFile << G4endl;
  
  
  outFile << sTABLE << '"' << "80%" << '"' << " > " << G4endl;
  
  // Raw #1
  outFile << sTR;
  outFile << sTD << sLFONT << "Code" << eLFONT<< eTD; 
  outFile << sTD << sLFONT << "Name" << eLFONT<< eTD; 
  outFile << sTD << sLFONT << "Mass" << eLFONT << eTD;
  outFile << sTD << sLFONT << "Charge" << eLFONT << eTD;
  outFile << sTD << sLFONT << "Life Time" << eLFONT << eTD;
  outFile << sTD << sLFONT << "Anti-Particle" << eLFONT<< eTD; 
  outFile << eTR << G4endl;;
  
  // Raw #2
  outFile << sTR;
  outFile << sTD << " " << eTD; 
  outFile << sTD << " " << eTD; 
  outFile << sTD <<  " [GeV/c" << sSUP << "2" << eSUP << "]" << eTD;
  outFile << sTD << " " << eTD; 
  outFile << sTD <<  " [ns]" << eTD;
  outFile << sTD << " " << eTD; 
  outFile << eTR << G4endl;;

  for (size_t i=0; i< pList.size(); i++){
    if (pList[i]->GetPDGEncoding()<0) continue;

    outFile << sTR << G4endl;;
    // column 1  : endcoding
    outFile << sTD << pList[i]->GetPDGEncoding() << eTD << G4endl;; 
    // column 2  : name 
    G4String name = pList[i]->GetParticleName();
    
    G4String fname = name +".html";
    // exception
    if (name == "J/psi") fname = "jpsi.html";
    
    outFile << sTD;
    outFile << "<A HREF=" << '"' << fname << '"' << ">";
    outFile << name << "</A>" << eTD << G4endl;
   
    // column 3 mass
    outFile << sTD <<  pList[i]->GetPDGMass()/GeV << eTD << G4endl;

    // column 4 charge
    outFile << sTD <<  pList[i]->GetPDGCharge()/eplus << eTD << G4endl;

    // column 5 life time
    outFile << sTD <<  pList[i]->GetPDGLifeTime()/ns << eTD << G4endl;
    
    // column 6 AntiParticle
    if  ( (pList[i]->GetAntiPDGEncoding()!= 0) &&
          (pList[i]->GetAntiPDGEncoding() != pList[i]->GetPDGEncoding() ) ) {
      G4ParticleDefinition* anti_particle  = G4ParticleTable::GetParticleTable()->FindParticle( pList[i]->GetAntiPDGEncoding() ); 
      
      outFile << sTD <<  anti_particle->GetParticleName() << eTD << G4endl;;
    }

    // end raw
    outFile << eTR << G4endl;;
  }
  
  outFile << eTABLE << G4endl;
  
  // footer
  PrintFooter(outFile); 
  
}

 void  G4HtmlPPReporter::GeneratePropertyTable(const G4ParticleDefinition* particle)
{
  if (particle->GetPDGEncoding()<0) return;

  G4String name = particle->GetParticleName();
  //--- open index file -----
  G4String fileName = baseDir + name + ".html";
  // exception
  if (name == "J/psi") fileName = baseDir +"jpsi.html";
  std::ofstream outFile(fileName, std::ios::out );
  outFile.setf( std::ios:: scientific, std::ios::floatfield );
  outFile << std::setprecision(7) << G4endl;

  PrintHeader(outFile);
   
  // particle name
  outFile << "<H2>" << name << "</H2>" << G4endl;
  outFile << "<HR>" << G4endl;

  // encoding, type
  outFile << sTABLE << '"' << "40%" << '"' << " > " << G4endl;
  outFile << sTR << sTD << sB << "PDG encoding" << eB << eTD;
  outFile << sTD <<  particle->GetPDGEncoding() << eTD << eTR << G4endl;
  outFile << sTR << sTD << sB << "Type" << eB << eTD;
  outFile << sTD <<  particle->GetParticleType() << eTD << eTR << G4endl;
  outFile << eTABLE << G4endl;
  outFile << "<HR>" << G4endl;

  // Properties  
  outFile << sTABLE << '"' << "60%" << '"' << " > " << G4endl;
  // mass
  outFile << sTR << sTD << sB << "Mass" << eB << eTD;
  outFile << sTD <<  particle->GetPDGMass()/GeV;
  outFile << " [GeV/c" << sSUP << "2" << eSUP << "]" << eTD << eTR << G4endl;
  // width
  outFile << sTR << sTD << sB << "Width" << eB << eTD;
  outFile << sTD <<  particle->GetPDGWidth()/GeV;
  outFile << " [GeV/c" << sSUP << "2" << eSUP << "]" << eTD << eTR << G4endl;
  // IJPC
  outFile << sTR << sTD << sB << "I J" << sSUP << "PC"<< eSUP << eB << eTD;
  if ( particle->GetPDGiIsospin() <0 ) {
    outFile << sTD << "    ";
  } else if ( particle->GetPDGiIsospin() == 1) {
    outFile << sTD << "1/2 ";
  } else if ( particle->GetPDGiIsospin() == 3) {
    outFile << sTD << "3/2 "; 
  } else {
    outFile << sTD << particle->GetPDGiIsospin()/2 << " ";
  }
  if ( particle->GetPDGiSpin() == 1) {
    outFile << "1/2";
  } else if ( particle->GetPDGiSpin() == 3) {
    outFile << "3/2"; 
  } else if ( particle->GetPDGiSpin() == 5) {
    outFile << "5/2"; 
  } else if ( particle->GetPDGiSpin() == 7) {
    outFile << "7/2"; 
  } else if ( particle->GetPDGiSpin() == 9) {
    outFile << "9/2"; 
  } else if ( particle->GetPDGiSpin() == 11) {
    outFile << "11/2"; 
  } else if ( particle->GetPDGiSpin() == 13) {
    outFile << "13/2"; 
  } else {
    outFile << particle->GetPDGiSpin()/2;
  }
  outFile << sSUP << sSYMBOL;
  if (particle->GetPDGiParity() == +1 ){
	outFile << "+";
  } else if (particle->GetPDGiParity() == -1 ){
	outFile << "-";
  } else {
    outFile << " ";
  }
  if (particle->GetPDGiConjugation() == +1 ){
	outFile << "+";
  } else if (particle->GetPDGiConjugation() == -1 ){
	outFile << "-";
  } else {
    outFile << " ";
  }
  outFile << eSYMBOL << eSUP;
  outFile <<  eTD << eTR << G4endl;
  // charge
  outFile << sTR << sTD << sB << "Charge" << eB << eTD;
  outFile << sTD <<  particle->GetPDGCharge()/eplus;
  outFile << eTD << eTR << G4endl;
  // Magnetic Moment
  outFile << sTR << sTD << sB << "Magnetic Moment" << eB << eTD;
  if  (particle->GetPDGMagneticMoment() != 0.0){
    outFile << sTD <<  particle->GetPDGMagneticMoment()/MeV*tesla;
    outFile << "[MeV/T]" << eTD << eTR << G4endl;
  } else {
    outFile << sTD <<  " not defined ";
    outFile << eTD << eTR << G4endl;
  }
  // life time
  outFile << sTR << sTD << sB << "Life Time" << eB << eTD;
  if ( particle->GetPDGLifeTime() >0.0 ) {
    outFile << sTD <<  particle->GetPDGLifeTime()/second;
    outFile << "[sec]" << eTD << G4endl;
  } else {
    if (particle->GetPDGStable()) {
      outFile << sTD << "stable" << eTD;
    } else if (particle->IsShortLived()) {
      outFile << sTD << "short-lived" << eTD;
    } else {
      outFile << sTD <<  "not Defined" << eTD;
    }
  }
  outFile << eTR << G4endl;

  outFile << eTABLE << G4endl;
  outFile << "<HR>" << G4endl;

  // Qurak content 
  outFile << "<H2>" << " Quark Content " << "</H2>" << G4endl;

  outFile << sTABLE << '"' << "60%" << '"' << " > " << G4endl;
 
  outFile << sTR;
  outFile << sTD << sB << "flavour " << eB << eTD ;
  outFile << sTD << sB << " quark " << eB << eTD; 
  outFile << sTD << sB << " anti-quark " << eB << eTD; 
  outFile << eTR;

  static const char* quarkName[6] = { "d", "u", "s", "c", "b", "t" };
  for (G4int flv = 0; flv <6 ; flv++ ){
    outFile << sTR;
    outFile << sTD << sB << quarkName[flv] << eB << eTD ;
    outFile << sTD << sB << particle->GetQuarkContent(flv+1) << eB << eTD ;
    outFile << sTD << sB << particle->GetAntiQuarkContent(flv+1) << eB << eTD ;
    outFile << eTR;
  }
  outFile << eTABLE << G4endl;
  outFile << "<HR>" << G4endl;

 // Decay Table  
  G4DecayTable* dcyTable = particle->GetDecayTable(); 
  if (dcyTable != 0) { 
    outFile << "<H2>" << " Decay Table " << "</H2>" << G4endl;

    outFile << sTABLE << '"' << "80%" << '"' << " > " << G4endl;
    
    outFile << sTR;
    outFile << sTD << sB << "BR" << eB << eTD ;
    outFile << sTD << sB << "kinematics" << eB << eTD; 
    outFile << eTR;

    for (G4int i=0; i< dcyTable->entries(); i++){
      G4VDecayChannel * channel = dcyTable->GetDecayChannel(i);
      outFile << sTR << G4endl;;
      // column 1  : BR
      outFile << sTD << channel->GetBR() << eTD;
      // column 2 : Kinematics
      outFile << sTD << channel->GetKinematicsName() << eTD;
      // column 3..  : daughters
      for (G4int j=0; j< channel->GetNumberOfDaughters(); j++){
        outFile << sTD << channel->GetDaughter(j)->GetParticleName() << eTD;
      }
      outFile << eTR << G4endl;
    }
    outFile << eTABLE << G4endl;
    outFile << "<HR>" << G4endl;
  }
  
  outFile << sB;
  outFile << "<A HREF=" << '"' <<  "index.html" << '"' << ">back to index</A>";
  outFile << eB << G4endl;

  PrintFooter(outFile); 
}

 void G4HtmlPPReporter::PrintHeader(std::ofstream& outFile)
{
  outFile << "<HTML>" << G4endl;
  outFile << "<HEAD>" << G4endl;
  outFile << " <META HTTP-EQUIV=" << "\"" << " Content-Type" << "\"";
  outFile << " CONTENT=" 
	  << "\"" << "text/html; charset=iso-8859-1" 
	  << "\"" << ">" << G4endl;
  outFile << " <TITLE>Geant4 Particle List </TITLE>" << G4endl;
  outFile << "</HEAD>" << G4endl;
  outFile << "<! -- Generated automatically by Geant4, "
	  << " -- !>" << G4endl;
  outFile << "<BODY>" << G4endl;
}

 void G4HtmlPPReporter::PrintFooter(std::ofstream& outFile)
{
  outFile << "<HR>" << G4endl;
  outFile << "</BODY>" << G4endl;
  outFile << "</HTML>" << G4endl;
}

const char*  G4HtmlPPReporter::sTABLE = "<TABLE WIDTH=";
const char*  G4HtmlPPReporter::eTABLE = "</TABLE>";
const char*  G4HtmlPPReporter::sTR = "<TR>";
const char*  G4HtmlPPReporter::eTR = "</TR>";
const char*  G4HtmlPPReporter::sTD = "<TD>";
const char*  G4HtmlPPReporter::eTD = "</TD>";
const char*  G4HtmlPPReporter::sB = "<B>";
const char*  G4HtmlPPReporter::eB = "</B>";
const char*  G4HtmlPPReporter::sLFONT = "<FONT SIZE = +1>";
const char*  G4HtmlPPReporter::eLFONT = "</FONT>";
const char*  G4HtmlPPReporter::sSYMBOL = "<FONT FACE = \"symbol\" >";
const char*  G4HtmlPPReporter::eSYMBOL = "</FONT>";
const char*  G4HtmlPPReporter::sSUP = "<SUP>";
const char*  G4HtmlPPReporter::eSUP = "</SUP>";
const char*  G4HtmlPPReporter::sSUB = "<SUB>";
const char*  G4HtmlPPReporter::eSUB = "</SUB>";
 
 
