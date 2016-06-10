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
// $Id: G4LevelReader.cc 88407 2015-02-18 09:18:44Z vnivanch $
//
// -------------------------------------------------------------------
//
//      GEANT4 header file 
//
//      File name:     G4NucLevel
//
//      Author:        V.Ivanchenko (M.Kelsey reading method is used)
// 
//      Creation date: 4 January 2012
//
//      Modifications:
//
// -------------------------------------------------------------------

#include "G4LevelReader.hh"
#include "G4NucleiProperties.hh"
#include "G4NucLevel.hh"
#include "G4SystemOfUnits.hh"
#include "G4Pow.hh"
#include <vector>
#include <fstream>
#include <sstream>

G4String G4LevelReader::fTrans[] = {
  "1-", "1+", "2-", "2+", "3-", "3+", "4-", "4+", "5-", "5+"};

G4LevelReader::G4LevelReader() 
  : fMinProbability(1.e-8),fVerbose(0)
{
  fTimeFactor = CLHEP::second/G4Pow::GetInstance()->logZ(2);
  char* directory = getenv("G4LEVELGAMMADATA");
  if(directory) {
    fDirectory = directory;
  } else {
    G4Exception("G4LevelReader()","had0707",FatalException,
		"Environment variable G4LEVELGAMMADATA is not defined");
    fDirectory = "";
  } 
  fFile = fDirectory + "/z100.a200";
  fPol = "  ";
  for(G4int i=0; i<10; ++i) { fICC[i] = 0.0; }
  for(G4int i=0; i<20; ++i) { buffer[i] = ' '; }
  bufp[0] = bufp[1] = ' ';

  fEnergy = fCurrEnergy = fTrEnergy = fProb = fTime = fSpin = fAlpha 
    = fNorm1 = fNorm2 = fNorm3 = 0.0;

  size_t nn = 10;
  vTransEnergy.reserve(nn);
  vGammaCumProbability.reserve(nn);
  vGammaECumProbability.reserve(nn);
  vGammaProbability.reserve(nn);
  vShellProbability.reserve(nn);
  vTrans.reserve(nn);

  nn = 100;
  vEnergy.reserve(nn);
  vTime.reserve(nn);
  vTimeg.reserve(nn);
  vSpin.reserve(nn);
  vLevel.reserve(nn);
}

G4LevelReader::~G4LevelReader() 
{}

const G4LevelManager* 
G4LevelReader::CreateLevelManager(G4int Z, G4int A)
{
  std::ostringstream ss;
  ss << "/z" << Z << ".a" << A;
  G4String st =  G4String(ss.str());
  fFile = fDirectory + st;
  return MakeLevelManager(Z, A, fFile);
}

const G4LevelManager* 
G4LevelReader::MakeLevelManager(G4int Z, G4int A, const G4String& filename)
{
  vEnergy.clear();
  vTime.clear();
  vTimeg.clear();
  vSpin.clear();
  vLevel.clear();

  vEnergy.push_back(0.0f);
  vTime.push_back(FLT_MAX);
  vTimeg.push_back(FLT_MAX);
  vSpin.push_back(0);
  vLevel.push_back(0);

  std::ifstream infile(filename, std::ios::in);

  // file is not opened
  if (!infile.is_open()) {
    if (fVerbose > 0) {
      G4cout << " G4LevelReader: fail open file for Z= " 
	     << Z << " A= " << A 
	     << " <" << filename << ">" <<  G4endl;
    }

  } else {

    if (fVerbose > 0) {
      G4cout << "G4LevelReader: open file for Z= " 
	     << Z << " A= " << A 
	     << " <" << filename << ">" <<  G4endl;
    }
    // read line by line
    G4bool end = false;
    G4bool next = true;
    G4int nline = 0;
    fCurrEnergy = DBL_MAX;
    do {

      fNorm3 = 0.0;

      G4bool isOK = (ReadDataItem(infile,fEnergy)   &&
		     ReadDataItem(infile,fTrEnergy) &&
		     ReadDataItem(infile,fProb)     &&
		     ReadDataItem(infile,fPol)      &&
		     ReadDataItem(infile,fTime)     &&
		     ReadDataItem(infile,fSpin)     &&
		     ReadDataItem(infile,fAlpha));

      fEnergy   *= CLHEP::keV;
      fTrEnergy *= CLHEP::keV;

      if(isOK) {
	for(G4int i=0; i<10; ++i) { 
	  isOK = (isOK && (ReadDataItem(infile,fICC[i])));
	}
      }
      if(!isOK) {
	end = true; 
	next = false; 
      }
      if(!isOK && fVerbose > 1) {
	G4cout << "Line #" << nline << " in file <" << filename
	       << " is corrupted Z= " << Z << " A= " << A << G4endl; 
	G4cout << "#Line " << nline << " " << fEnergy << " " << fTrEnergy 
	       << " " << fProb << " " <<  fPol << " " << fTime << " " 
	       << fSpin << " " << fAlpha << G4endl;
	G4cout << "      ";
	for(G4int i=0; i<10; ++i) { G4cout << fICC[i] << " "; }
	G4cout << G4endl;
      }
      // end of nuclear level data
      if(end || fEnergy > fCurrEnergy) {
	size_t nn = vTransEnergy.size();
        if(fVerbose > 1) {
	  G4cout << "Reader: new level E= " << fCurrEnergy 
		 << " Ntransitions= " << nn << " fNorm1= " << fNorm1 
		 << " fNorm2= " << fNorm2 << G4endl;
	} 
	if(nn > 0) {
	  if(fNorm1 > 0.0) {
	    --nn;
	    vTimeg.push_back((G4float)(fTime*fNorm2*fTimeFactor/fNorm1));
	    fNorm1 = 1.0/fNorm1; 
	    fNorm2 = 1.0/fNorm2; 
	    for(size_t i=0; i<nn; ++i) {
	      vGammaCumProbability[i]  = (G4float)(vGammaCumProbability[i]*fNorm1);
	      vGammaECumProbability[i] = (G4float)(vGammaECumProbability[i]*fNorm2);
	      /*
              G4cout << "Probabilities[" << i << "]= " << vGammaCumProbability[i]
		     << "  " << vGammaECumProbability[i] 
		     << " Etran= " << vTransEnergy[i]
		     << G4endl;
	      */
	    }
	    vGammaCumProbability[nn]  = 1.0f;
	    vGammaECumProbability[nn] = 1.0f;
	    /*
	    G4cout << "Probabilities[" << nn << "]= " << vGammaCumProbability[nn]
		   << "  " << vGammaECumProbability[nn] 
		   << " Etran= " << vTransEnergy[nn]
		   << G4endl;
	    */
	    //case of X-level
	  } else {
	    vGammaCumProbability[0]  = 0.0f;
	    vGammaECumProbability[0] = 0.0f;
	  }
      
	  vLevel.push_back(new G4NucLevel(vTransEnergy,
					  vGammaCumProbability,
					  vGammaECumProbability,
					  vGammaProbability,
					  vTrans,
					  vShellProbability));
	}
        if(!end) { next = true; }
      }
      // begin nuclear level data
      if(next) {
        //G4cout << "== Reader: begin of new level E= " << fEnergy << G4endl; 
	fCurrEnergy = fEnergy;
	vEnergy.push_back((G4float)fEnergy);
	vTime.push_back((G4float)(fTime*fTimeFactor));
        if(fSpin > 20.0) { fSpin = 0.0; }
        fProb = std::max(fProb, fMinProbability);  
	vSpin.push_back(G4lrint(2*fSpin));
        fNorm1 = 0.0;
        fNorm2 = 0.0;
	vTransEnergy.clear();
	vGammaCumProbability.clear();
	vGammaECumProbability.clear();
	vGammaProbability.clear();
	vShellProbability.clear();
	vTrans.clear();
        next = false;
      } 
      // continue filling level data
      if(!end) {
        if(fProb > 0.0) {
	  // by default transition to a ground state
          G4float efinal = (G4float)(fEnergy - fTrEnergy);
          G4float elevel = 0.0f;
	  // do not check initial energy
	  G4int nn = vEnergy.size() - 1;
	  if(0 < nn) {
	    G4float ediffMin = std::abs(efinal);
	    for(G4int i=0; i<nn; ++i) {
              G4float ediff = std::abs(efinal - vEnergy[i]);
	      /*
	      G4cout << "Elevel[" << i << "]= " << vEnergy[i] 
		     << " Efinal= " << efinal
		     << " Ediff= " << ediff
		     << " EdiffMin= " << ediffMin << G4endl;
	      */
	      if(ediff < ediffMin) {
		ediffMin = ediff;
                elevel = vEnergy[i];
	      }
	    }
	  }
	  G4double x = 1.0 + fAlpha;
	  fNorm1 += fProb;
	  fNorm2 += fProb*x;
	  vTransEnergy.push_back(elevel);
	  vGammaCumProbability.push_back((G4float)(fNorm1));
	  vGammaECumProbability.push_back((G4float)(fNorm2));
	  vGammaProbability.push_back((G4float)(1.0/x));
	  G4int tnum = 0;
	  for(; tnum<10; ++tnum) { if(fTrans[tnum] == fPol) { break; } }
	  vTrans.push_back(tnum);
	  if(fAlpha > 0.0) {
	    vShellProbability.push_back(NormalizedICCProbability(Z));
	  } else {
	    vShellProbability.push_back(nullptr);
	  }
	}
      }
      ++nline;
      if(nline > 10000) {
	G4cout << "G4LevelReader: Line #" << nline << " Z= " << Z << " A= "
	       << " this file is too long - stop loading" << G4endl;
	end = true;
      }
      // Loop checking, 05-Aug-2015, Vladimir Ivanchenko
    } while (!end);
    infile.close();
  }

  const G4LevelManager* man = new G4LevelManager(vEnergy,
						 vTime,vTimeg,
						 vSpin,vLevel); 
  if(fVerbose > 0) {
    G4cout << "Reader: new manager for Z= " << Z << " A= " << A 
	   << " Nlevels= " << vEnergy.size() << " E[0]= " 
	   << vEnergy[0]/CLHEP::MeV << " MeV   Emax= " 
	   << man->MaxLevelEnergy()/CLHEP::MeV << " MeV "
	   << " " << fVerbose 
	   << G4endl;
  } 
  return man;
}

G4bool G4LevelReader::ReadDataItem(std::istream& dataFile, G4double& x)
{
  x = 0.0;
  for(G4int i=0; i<20; ++i) { buffer[i] = ' '; }
  G4bool okay = true;
  dataFile >> buffer;
  if(dataFile.fail()) { okay = false; }
  else { x = strtod(buffer, nullptr); }

  return okay;
}

G4bool G4LevelReader::ReadDataItem(std::istream& dataFile, G4String& x)
{
  G4bool okay = true;
  dataFile >> bufp;
  if(dataFile.fail()) { okay = false; }
  else { x = G4String(bufp, 2); }
  return okay;
}

const std::vector<G4float>* G4LevelReader::NormalizedICCProbability(G4int Z)
{
  std::vector<G4float>* vec = nullptr;
  G4int LL = 3;
  G4int M = 5;
  G4int N = 1;
  if(Z <= 4) {
    LL = 1;
    M = N = 0;
  } else if(Z <= 6) {
    LL = 2;
    M = N = 0;
  } else if(Z <= 10) {
    M = N = 0;
  } else if(Z <= 12) {
    M = 1;
    N = 0;
  } else if(Z <= 17) {
    M = 2;
    N = 0;
  } else if(Z == 18) {
    M = 3;
    N = 0;
  } else if(Z <= 20) {
    M = 3;
  } else if(Z <= 27) {
    M = 4;
  }
  G4double norm = 0.0;
  if(LL < 3) { for(G4int i=LL+1; i<=4; ++i) { fICC[i] = 0.0; } }
  if(M < 5) { for(G4int i=M+4; i<=8; ++i) { fICC[i] = 0.0; } }
  if(N < 1) { fICC[9] = 0.0; }
  for(G4int i=0; i<10; ++i) {
    norm += fICC[i]; 
    fICC[i] = norm;
  }
  if(norm > 0.0) {
    norm = 1.0/norm;
    vec = new std::vector<G4float>;
    vec->reserve(10);
    for(G4int i=0; i<9; ++i) {
      vec->push_back((G4float)(fICC[i]*norm)); 
    }
    vec->push_back(1.0f);
  }
  return vec;
}
