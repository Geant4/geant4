// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: tst2HtmlReporter.cc,v 1.1 1999-06-17 04:45:31 kurasige Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ---------------------------------------------------------------
#include "tst2HtmlReporter.hh"
#include "G4ios.hh"
#include "globals.hh"
#include "G4DecayTable.hh"  
#include <fstream.h>
#include <iomanip.h>
#include <rw/ctoken.h>

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
  RWCTokenizer savedToken( option );
  
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
  ofstream outFile(fileName, ios::out );
  outFile.setf( ios:: scientific, ios::floatfield );

  // header
  PrintHeader(outFile);
  
  // comment
  outFile << "<! -- " << comment << " --!> " << endl;
  outFile << endl;
 
  
  outFile << sTABLE << '"' << "80%" << '"' << " > " << endl;

  // Raw #1
  outFile << sTR;
  outFile << sTD << sLFONT << "Code" << eLFONT<< eTD; 
  outFile << sTD << sLFONT << "Name" << eLFONT<< eTD; 
  outFile << sTD << sLFONT << "Anti-Particle" << eLFONT<< eTD; 
  outFile << eTR << endl;;

  for (G4int i=0; i< entries(); i++){
	outFile << sTR << endl;;
	// column 1  : endcoding
    outFile << sTD << GetEncoding(i) << eTD << endl;; 
    // column 2  : name 
    G4String name = GetParticle(i)->GetParticleName();
    outFile << sTD;
	outFile << "<A HREF=" << '"' << name << ".html" << '"' << ">";
    outFile << name << "</A>" << eTD << endl;;    
	// column 3 AntiParticle
	if  ( (GetAntiParticle(i)!= 0) &&
          (GetParticle(i)->GetAntiPDGEncoding() != GetEncoding(i) ) ) {
      outFile << sTD <<  GetAntiParticle(i)->GetParticleName() << eTD << endl;;
    }
    outFile << eTR << endl;;
  }
  
  outFile << eTABLE << endl;

  // footer
  PrintFooter(outFile); 

}

 void  tst2HtmlReporter::GeneratePropertyTable(G4ParticleDefinition* particle)
{
  G4String name = particle->GetParticleName();
  //--- open index file -----
  G4String fileName = baseDir + name + ".html";
  ofstream outFile(fileName, ios::out );
  outFile.setf( ios:: scientific, ios::floatfield );
  outFile << setprecision(7) << endl;

  PrintHeader(outFile);
   
  // particle name
  outFile << "<H2>" << name << "</H2>" << endl;
  outFile << "<HR>" << endl;

  // encoding, type
  outFile << sTABLE << '"' << "40%" << '"' << " > " << endl;
  outFile << sTR << sTD << sB << "PDG encoding" << eB << eTD;
  outFile << sTD <<  particle->GetPDGEncoding() << eTD << eTR << endl;
  outFile << sTR << sTD << sB << "Type" << eB << eTD;
  outFile << sTD <<  particle->GetParticleType() << eTD << eTR << endl;
  outFile << eTABLE << endl;
  outFile << "<HR>" << endl;

  // Properties  
  outFile << sTABLE << '"' << "60%" << '"' << " > " << endl;
  // mass
  outFile << sTR << sTD << sB << "Mass" << eB << eTD;
  outFile << sTD <<  particle->GetPDGMass()/GeV;
  outFile << " [GeV/c" << sSUP << "2" << eSUP << "]" << eTD << eTR << endl;
  // width
  outFile << sTR << sTD << sB << "Width" << eB << eTD;
  outFile << sTD <<  particle->GetPDGWidth()/GeV;
  outFile << " [GeV/c" << sSUP << "2" << eSUP << "]" << eTD << eTR << endl;
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
  outFile <<  eTD << eTR << endl;
  // charge
  outFile << sTR << sTD << sB << "Charge" << eB << eTD;
  outFile << sTD <<  particle->GetPDGCharge()/eplus;
  outFile << eTD << eTR << endl;
  // life time
  outFile << sTR << sTD << sB << "Life Time" << eB << eTD;
  if ( particle->GetPDGLifeTime() >0.0 ) {
    outFile << sTD <<  particle->GetPDGLifeTime()/second;
    outFile << "[sec]" << eTD << endl;
  } else {
    if (particle->GetPDGStable()) {
      outFile << sTD << "stable" << eTD;
    } else if (particle->IsShortLived()) {
      outFile << sTD << "short-lived" << eTD;
    } else {
      outFile << sTD <<  "not Defined" << eTD;
    }
  }
  outFile << eTR << endl;

  outFile << eTABLE << endl;
  outFile << "<HR>" << endl;

  // Qurak content 
  outFile << "<H2>" << " Quark Content " << "</H2>" << endl;

  outFile << sTABLE << '"' << "60%" << '"' << " > " << endl;
 
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
  outFile << eTABLE << endl;
  outFile << "<HR>" << endl;

 // Decay Table  
  G4DecayTable* dcyTable = particle->GetDecayTable(); 
  if (dcyTable != 0) { 
    outFile << "<H2>" << " Decay Table " << "</H2>" << endl;

    outFile << sTABLE << '"' << "80%" << '"' << " > " << endl;
    
    outFile << sTR;
    outFile << sTD << sB << "BR" << eB << eTD ;
    outFile << sTD << sB << "kinematics" << eB << eTD; 
    outFile << eTR;

    for (G4int i=0; i< dcyTable->entries(); i++){
      G4VDecayChannel * channel = dcyTable->GetDecayChannel(i);
      outFile << sTR << endl;;
      // column 1  : BR
      outFile << sTD << channel->GetBR() << eTD;
      // column 2 : Kinematics
      outFile << sTD << channel->GetKinematicsName() << eTD;
      // column 3..  : daughters
      for (G4int j=0; j< channel->GetNumberOfDaughters(); j++){
        outFile << sTD << channel->GetDaughter(j)->GetParticleName() << eTD;
      }
      outFile << eTR << endl;
    }
    outFile << eTABLE << endl;
    outFile << "<HR>" << endl;
  }
  
  outFile << sB;
  outFile << "<A HREF=" << '"' <<  "index.html" << '"' << ">back to index</A>";
  outFile << eB << endl;

  PrintFooter(outFile); 
}

 void tst2HtmlReporter::PrintHeader(ofstream& outFile)
{
  outFile << "<HTML>" << endl;
  outFile << "<HEAD>" << endl;
  outFile << " <META HTTP-EQUIV=" << "\"" << " Content-Type" << "\"";
  outFile << " CONTENT=" << "\"" << "text/html; charset=iso-8859-1" << "\"" << ">" << endl;
  outFile << " <TITLE>Geant4 Particle List </TITLE>" << endl;
  outFile << "</HEAD>" << endl;
  outFile << "<!-- Generated automatically by Geant4, ";
  outFile << " -- >" << endl;
  outFile << "<BODY>" << endl;
}

 void tst2HtmlReporter::PrintFooter(ofstream& outFile)
{
  outFile << "<HR>" << endl;
  outFile << "</BODY>" << endl;
  outFile << "</HTML>" << endl;
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
 
 
