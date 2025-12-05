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
#include "G4NuclearLevelData.hh"
#include "G4DeexPrecoParameters.hh"
#include "G4SystemOfUnits.hh"
#include "G4Pow.hh"
#include <fstream>
#include <sstream>

namespace
{
  const G4int countmax = 4;
  const G4int nfloting = 13;
  const G4double eTolarence = 2*CLHEP::eV;
  const G4String fFloatingLevels[13] = {
  "-", "+X", "+Y", "+Z", "+U", "+V", "+W", "+R", "+S", "+T", "+A", "+B", "+C"};
}

G4LevelReader::G4LevelReader(G4NuclearLevelData* ptr) 
  : fData(ptr)
{
  fAlphaMax = (G4float)1.e15;
  fTimeFactor = CLHEP::second/G4Pow::GetInstance()->logZ(2);
  fDirectory = G4String(G4FindDataDir("G4LEVELGAMMADATA"));
  if (fDirectory.empty()) {
    G4Exception("G4LevelReader::G4LevelReader()", "had014", FatalException,
      "G4LEVELGAMMADATA environment variable not set");
  }

  vTrans.resize(fTransMax,0);
  vRatio.resize(fTransMax,0.0f);
  vGammaCumProbability.resize(fTransMax,0.0f);
  vGammaProbability.resize(fTransMax,0.0f);
  vShellProbability.resize(fTransMax,nullptr);

  vEnergy.resize(fLevelMax,0.0);
  vSpin.resize(fLevelMax,0);
  vLevel.resize(fLevelMax,nullptr);
}

G4bool G4LevelReader::ReadData(std::istringstream& stream, G4double& x)
{
  stream >> x;
  return !stream.fail();
}

G4bool G4LevelReader::ReadDataItem(std::istream& dataFile, G4double& x)
{
  x = 0.0;
  for(G4int i=0; i<nbufmax; ++i) { buffer[i] = ' '; }
  G4bool okay = true;
  dataFile >> buffer;
  if(dataFile.fail()) { okay = false; }
  else { x = std::strtod(buffer, 0); }
  return okay;
}

G4bool G4LevelReader::ReadDataItem(std::istream& dataFile, G4float& x)
{
  x = 0.0f;
  for(G4int i=0; i<nbuf1; ++i) { buff1[i] = ' '; }
  G4bool okay = true;
  dataFile >> buff1;
  if(dataFile.fail()) { okay = false; }
  else { x = std::atof(buff1); }

  return okay;
}

G4bool G4LevelReader::ReadDataItem(std::istream& dataFile, G4int& ix)
{
  ix = 0;
  for(G4int i=0; i<nbuf2; ++i) { buff2[i] = ' '; }
  G4bool okay = true;
  dataFile >> buff2;
  if(dataFile.fail()) { okay = false; }
  else { ix = std::atoi(buff2); }

  return okay;
}

const std::vector<G4float>* G4LevelReader::NormalizedICCProbability(G4int Z)
{
  std::vector<G4float>* vec = nullptr;
  G4int LL = 3;
  G4int M = 5;
  G4int N = 1;
  G4int Kmax = 9;
  if(Z <= 27) {
    M = N = 0;
    if(Z <= 4) {
      LL = 1;
      Kmax = 2;
    } else if(Z <= 6) {
      LL = 2;
      Kmax = 3;
    } else if(Z <= 10) {
      Kmax = 4;
    } else if(Z <= 12) {
      M = 1;
      Kmax = 8;
    } else if(Z <= 17) {
      M = 2;
      Kmax = 8;
    } else if(Z == 18) {
      M = 3;
      Kmax = 8;
    } else if(Z <= 20) {
      M = 3;
      N = 1;
    } else {
      M = 4;
      N = 1;
    }
    if(LL < 3) { for(G4int i=LL+1; i<=4; ++i) { fICC[i] = 0.0f; } }
    if(M < 5)  { for(G4int i=M+4;  i<=8; ++i) { fICC[i] = 0.0f; } }
    if(N < 1)  { fICC[9] = 0.0f; }
  }
  G4float norm = 0.0f;
  for (G4int i = 0; i <= Kmax; ++i) {
    norm += fICC[i]; 
    fICC[i] = norm;
  }
  if (norm > 0.0f) {
    norm = 1.0f/norm;
  }
  vec = new std::vector<G4float>(Kmax + 1, 0.0f);
  for (G4int i = 0; i < Kmax; ++i) {
    (*vec)[i] = fICC[i]*norm;
  }
  (*vec)[Kmax] = 1.0f;
  if (fVerbose > 3) {
    G4long prec = G4cout.precision(3);
    G4cout << "# InternalConv: ";
    for (G4int i = 0; i <= Kmax; ++i) { G4cout << " " << (*vec)[i]; }
    G4cout << G4endl;
    G4cout.precision(prec);
  }
  return vec;
}

const G4LevelManager* 
G4LevelReader::CreateLevelManager(G4int Z, G4int A)
{
  std::ostringstream ss;
  ss << fDirectory << "/z" << Z << ".a" << A;
  std::ifstream infile(ss.str(), std::ios::in);

  // file is not opened
  if (!infile.is_open()) {
    if(fVerbose > 1) {
      G4ExceptionDescription ed;
      ed << "Regular file " << ss.str() << " is not opened! Z="
         << Z << " A=" << A;  
      G4Exception("G4LevelReader::LevelManager(..)","had014",
		  JustWarning, ed, "Check file path");
    }
    return nullptr;
  }
  // file is opened
  if (fVerbose > 1) {
    G4cout << "G4LevelReader: open file " << ss.str() << " for Z= " 
	   << Z << " A= " << A <<  G4endl;
  }
  return LevelManager(Z, A, infile);
}

const G4LevelManager* 
G4LevelReader::MakeLevelManager(G4int Z, G4int A, const G4String& filename)
{
  std::ifstream infile(filename, std::ios::in);
  
  // file is not opened
  if (!infile.is_open()) {
    if(fVerbose > 1) {
      G4ExceptionDescription ed;
      ed << "External file " << filename << " is not opened! Z=" 
         << Z << " A=" << A;  
      G4Exception("G4LevelReader::LevelManager(..)","had014",
		  FatalException, ed, "Check file path");
    }
    return nullptr;
  }
  // file is opened
  if (fVerbose > 1) {
    G4cout << "G4LevelReader: open external file " << filename 
           << " for Z= " << Z << " A= " << A << G4endl;
  }
  return LevelManager(Z, A, infile);
}

const G4LevelManager* 
G4LevelReader::LevelManager(G4int Z, G4int A, std::ifstream& infile)
{
  G4bool allLevels = fData->GetParameters()->StoreICLevelData();
  fPol = "  ";
  G4int i = 0;
  for (;;) {
    infile >> i1 >> fPol;    // Level number and floating level
    // normal end of file
    if (infile.eof()) {
      if (fVerbose > 1) { 
	G4cout << "### End of file Z= " << Z << " A= " << A 
	       << " Nlevels= " << i << G4endl;
      }
      break;
    }
    // start reading new level data
#ifdef G4VERBOSE
    if(fVerbose > 2) { 
      G4cout << "New line: i1= " << i1 << "  fPol= <" << fPol << "> " << G4endl;
    }
#endif
    // read new level data
    if (!(ReadDataItem(infile, fEnergy) &&
 	  ReadDataItem(infile, fTime) &&
	  ReadDataItem(infile, fSpin) &&
	  ReadDataItem(infile, ntrans))) {
      if (fVerbose > 1) { 
	G4cout << "### Incomplete end of file Z= " << Z << " A= " << A 
	       << " Nlevels= " << i << G4endl;
      }
      break;
    }
    fEnergy *= CLHEP::keV;
    for (k=0; k<nfloting; ++k) {
      if (fPol == fFloatingLevels[k]) {
	break;
      }
    }
    // if a previous level has higher energy the current should be ignored
    // data with wrong level should red anyway
    if (0 < i) {
      if (fEnergy < vEnergy[i-1]) {
#ifdef G4VERBOSE
	++count1;
	if (count1 < countmax && fVerbose > 0) {
	  G4cout << "### G4LevelReader: broken level " << i
		 << " E(MeV)= " << fEnergy << " < " << vEnergy[i-1]
		 << " for isotope Z= " << Z << " A= " 
		 << A << " level energy increased" << G4endl;
	}
#endif
	// for any case
	fEnergy = vEnergy[i-1] + eTolarence;
      }
    }
    vEnergy[i] = fEnergy;
    if (fTime > 0.0)  { fTime *= fTimeFactor; }
    else if (fTime < 0.0) { fTime = DBL_MAX; }

    G4int twos = G4lrint(fSpin + fSpin);
    twos = std::max(twos, -100);
    vSpin[i] = 100 + twos + k*100000;
#ifdef G4VERBOSE
    if (fVerbose > 2) {
      G4cout << "   Level #" << i1 << " E(MeV)=" << fEnergy/CLHEP::MeV
	     << "  LTime(s)=" << fTime << " 2S=" << vSpin[i]
	     << "  meta=" << vSpin[i]/100000 << " idx=" << i  
	     << " ntr=" << ntrans << G4endl;
    }
#endif
    vLevel[i] = nullptr;
    if (ntrans == 0) {
      vLevel[i] = new G4NucLevel(0, fTime, vTrans,
				 vGammaCumProbability,
				 vGammaProbability,
				 vRatio,
				 vShellProbability);
    } else if (ntrans > 0) {

      G4bool isTransOK = true;
      if (ntrans > fTransMax) {
	fTransMax = ntrans;
	vTrans.resize(fTransMax);
	vRatio.resize(fTransMax);
	vGammaCumProbability.resize(fTransMax);
	vGammaProbability.resize(fTransMax);
	vShellProbability.resize(fTransMax);
      }
      fNorm1 = 0.0f;
      for (G4int j=0; j<ntrans; ++j) {
       
	if (!(ReadDataItem(infile, i2) &&
	      ReadDataItem(infile, fTransEnergy) &&
	      ReadDataItem(infile, fProb) &&
	      ReadDataItem(infile, tnum) &&
	      ReadDataItem(infile, vRatio[j]) &&
	      ReadDataItem(infile, fAlpha))) {
#ifdef G4VERBOSE
	  ++count2;
	  if (count2 < countmax && fVerbose > 0) {
	    G4cout << "### Fail to read transition j=" << j
		   << " j=" << j << " i2=" << i2 
		   << " Z=" << Z << " A=" << A << G4endl; 
	  }
#endif
	  isTransOK = false;
	}
        if (i2 >= i) {
#ifdef G4VERBOSE
	  ++count2;
	  if (count2 < countmax) {
	    G4cout << "### G4LevelReader: broken transition " << j 
		   << " from level " << i << " to " << i2
		   << " for isotope Z= " << Z << " A= " 
		   << A << "; the transition probability set to zero" << G4endl; 
	  }
#endif
	  isTransOK = false;
	  fProb = 0.0f;
	}
	vTrans[j] = i2*10000 + tnum;
        fAlpha = std::min(std::max(fAlpha,0.f), fAlphaMax);
	G4float x = 1.0f + fAlpha;
	fNorm1 += x*fProb;
	vGammaCumProbability[j] = fNorm1;
	vGammaProbability[j] = 1.0f/x;
	vShellProbability[j] = nullptr;
	if (fVerbose > 2) { 
	  G4long prec = G4cout.precision(4);
	  G4cout << "### Transition #" << j << " to level " << i2 
		 << " i2= " << i2 <<  " Etrans(MeV)= " << fTransEnergy*CLHEP::keV
		 << "  fProb= " << fProb << " MultiP= " << tnum
		 << "  fMpRatio= " << fRatio << " fAlpha= " << fAlpha 
		 << G4endl;
	  G4cout.precision(prec);
	}
	if (fAlpha > 0.0f) {
	  for (k=0; k<10; ++k) {
	    if (!ReadDataItem(infile,fICC[k])) {
	      isTransOK = false;
#ifdef G4VERBOSE
	      ++count2;
	      if (count2 < countmax) { 
		G4cout << "### G4LevelReader: fail to read conversion coeff k= " << k 
		       << " for transition j= " << j 
		       << "  Z= " << Z << " A= " << A << G4endl; 
	      }
#endif
	      for(kk=k; kk<10; ++kk) { fICC[kk] = 0.f; }
	    }
	  }
	  if (allLevels) { 
	    vShellProbability[j] = NormalizedICCProbability(Z);
	  }
	}
      }
      if (ntrans > 0) {
        G4int nt = ntrans - 1;      
        if (fVerbose > 2) {
          G4cout << "=== New G4NucLevel: Ntrans=" << ntrans 
	         << " Time(ns)=" << fTime
	         << " IdxTrans=" << vTrans[nt]/10000
		 << " isOK=" << isTransOK
	         << G4endl;
        }
	fNorm1 = (FLT_MIN < fNorm1) ? 1.0f/fNorm1 : 0.0f; 
	for (k=0; k<nt; ++k) {
	  vGammaCumProbability[k] *= fNorm1;
#ifdef G4VERBOSE
	  if (fVerbose > 3) {
	    G4cout << "Probabilities[" << k 
		   << "]= " << vGammaCumProbability[k]
		   << "  " << vGammaProbability[k]
		   << " idxTrans= " << vTrans[k]/10000
		   << G4endl;
	  }
#endif
	}
	vGammaCumProbability[nt] = 1.0f;
	vLevel[i] = new G4NucLevel((std::size_t)ntrans, fTime, vTrans,
				   vGammaCumProbability,
				   vGammaProbability,
				   vRatio,
				   vShellProbability);
      }
    }
    ++i;
    if (i == fLevelMax) {
      fLevelMax += 10;
      vEnergy.resize(fLevelMax, 0.0);
      vSpin.resize(fLevelMax, 0);
      vLevel.resize(fLevelMax, nullptr);
    }
  }
  G4LevelManager* lman = nullptr;
  if (1 <= i) {
    lman = new G4LevelManager(Z, A, (std::size_t)i, vEnergy, vSpin, vLevel);
    if (fVerbose > 1) {
      G4cout << "=== Reader: new manager for Z=" << Z << " A=" << A 
	     << " Nlevels=" << i << " E[0]=" 
	     << vEnergy[0]/CLHEP::MeV << " MeV   E1=" 
	     << vEnergy[i-1]/CLHEP::MeV << " MeV"
	     << " count1,2=" << count1 << ", " << count2 << G4endl;
    }
  }
  return lman;
}
