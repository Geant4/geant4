// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: tst2HtmlReporter.cc,v 1.4 1999-12-15 14:51:20 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ---------------------------------------------------------------
#include "tst2HtmlReporter.hh"
#include "G4ios.hh"
#include "globals.hh"
#include "G4DecayTable.hh"  
#include "g4std/fstream"
#include "g4std/iomanip"
#include "g4rw/ctoken.h"

 tst2HtmlReporter::tst2HtmlReporter():tst2VParticleReporter()
{
 
}

 tst2HtmlReporter::~tst2HtmlReporter()
{
}    

 void tst2HtmlReporter::Print(const tst2ParticleContainer& container, 
                               const G4String& option)
{
  pList = &container;
  
  SparseOption( option );

  GenerateIndex();

  for (G4int i=0; i< entries(); i++){
    GeneratePropertyTable(GetParticle(i));
  }
}    


void tst2HtmlReporter::SparseOption(const G4String& option)
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

 void tst2HtmlReporter::GenerateIndex()
{
  //--- open index file -----
  G4String fileName = baseDir + "index.html";
  G4std::ofstream outFile(fileName, G4std::ios::out );
  outFile.setf( G4std::ios:: scientific, G4std::ios::floatfield );

  // header
  PrintHeader(outFile);
  
  // comment
  outFile << "<! -- " << comment << " --!> " << G4endl;
  outFile << G4endl;
 
  
  outFile << sTABLE << '"' << "80%" << '"' << " > " << G4endl;

  // Raw #1
  outFile << sTR;
  outFile << sTD << sLFONT << "Code" << eLFONT<< eTD; 
  outFile << sTD << sLFONT << "Name" << eLFONT<< eTD; 
  outFile << sTD << sLFONT << "Anti-Particle" << eLFONT<< eTD; 
  outFile << eTR << G4endl;;

  for (G4int i=0; i< entries(); i++){
	outFile << sTR << G4endl;;
	// column 1  : endcoding
    outFile << sTD << GetEncoding(i) << eTD << G4endl;; 
    // column 2  : name 
    G4String name = GetParticle(i)->GetParticleName();
    outFile << sTD;
	outFile << "<A HREF=" << '"' << name << ".html" << '"' << ">";
    outFile << name << "</A>" << eTD << G4endl;;    
	// column 3 AntiParticle
	if  ( (GetAntiParticle(i)!= 0) &&
          (GetParticle(i)->GetAntiPDGEncoding() != GetEncoding(i) ) ) {
      outFile << sTD <<  GetAntiParticle(i)->GetParticleName() << eTD << G4endl;;
    }
    outFile << eTR << G4endl;;
  }
  
  outFile << eTABLE << G4endl;

  // footer
  PrintFooter(outFile); 

}

 void  tst2HtmlReporter::GeneratePropertyTable(G4ParticleDefinition* particle)
{
  G4String name = particle->GetParticleName();
  //--- open index file -----
  G4String fileName = baseDir + name + ".html";
  G4std::ofstream outFile(fileName, G4std::ios::out );
  outFile.setf( G4std::ios:: scientific, G4std::ios::floatfield );
  outFile << G4std::setprecision(7) << G4endl;

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

 void tst2HtmlReporter::PrintHeader(G4std::ofstream& outFile)
{
  outFile << "<HTML>" << G4endl;
  outFile << "<HEAD>" << G4endl;
  outFile << " <META HTTP-EQUIV=" << "\"" << " Content-Type" << "\"";
  outFile << " CONTENT=" << "\"" << "text/html; charset=iso-8859-1" << "\"" << ">" << G4endl;
  outFile << " <TITLE>Geant4 Particle List </TITLE>" << G4endl;
  outFile << "</HEAD>" << G4endl;
  outFile << "<!-- Generated automatically by Geant4, ";
  outFile << " -- >" << G4endl;
  outFile << "<BODY>" << G4endl;
}

 void tst2HtmlReporter::PrintFooter(G4std::ofstream& outFile)
{
  outFile << "<HR>" << G4endl;
  outFile << "</BODY>" << G4endl;
  outFile << "</HTML>" << G4endl;
}

const char*  tst2HtmlReporter::sTABLE = "<TABLE WIDTH=";
const char*  tst2HtmlReporter::eTABLE = "</TABLE>";
const char*  tst2HtmlReporter::sTR = "<TR>";
const char*  tst2HtmlReporter::eTR = "</TR>";
const char*  tst2HtmlReporter::sTD = "<TD>";
const char*  tst2HtmlReporter::eTD = "</TD>";
const char*  tst2HtmlReporter::sB = "<B>";
const char*  tst2HtmlReporter::eB = "</B>";
const char*  tst2HtmlReporter::sLFONT = "<FONT SIZE = +1>";
const char*  tst2HtmlReporter::eLFONT = "</FONT>";
const char*  tst2HtmlReporter::sSYMBOL = "<FONT FACE = \"symbol\" >";
const char*  tst2HtmlReporter::eSYMBOL = "</FONT>";
const char*  tst2HtmlReporter::sSUP = "<SUP>";
const char*  tst2HtmlReporter::eSUP = "</SUP>";
const char*  tst2HtmlReporter::sSUB = "<SUB>";
const char*  tst2HtmlReporter::eSUB = "</SUB>";
 
 
