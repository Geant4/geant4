//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//

#include "SimpleAnalyzer.hh"


int SimpleAnalyzer::mOutputFileCnt = 0;


void SimpleAnalyzer::analyze(const G4UppTrackVector& allTracks) const
{
  cout << "Analyzer running..." << endl;
  mOutputFileCnt++;
  string output_file = mOutputFile + (char)('0'+mOutputFileCnt);
  cout << "type filename: " << output_file << endl;
  ofstream out(output_file.c_str());
  ofstream outd((output_file+'d').c_str());
  ofstream outp((output_file+'p').c_str());
  for (int i=0; i<allTracks.size(); i++) {
    G4ParticleDefinition* aDefinition = allTracks[i]->GetDefinition();
    G4ThreeVector aPos = allTracks[i]->GetPosition();
    G4int encoding = aDefinition->GetPDGEncoding();
    if ( (aDefinition->GetParticleName()=="proton") || 
	 (aDefinition->GetParticleName()=="neutron") ) {
      out 
	<< aPos.x()/fermi << " " 
	<< aPos.y()/fermi << " " 
	<< aPos.z()/fermi << " "
	<< aDefinition->GetParticleName() << endl;
    }
    else if ( (encoding==211) || (encoding==-211) || (encoding==111) ) {
      outp
	<< aPos.x()/fermi << " " 
	<< aPos.y()/fermi << " " 
	<< aPos.z()/fermi << " "
	<< aDefinition->GetParticleName() << endl;
    } else
      outd
	<< aPos.x()/fermi << " " 
	<< aPos.y()/fermi << " " 
	<< aPos.z()/fermi << " "
	<< aDefinition->GetParticleName() << endl;
  }
}
