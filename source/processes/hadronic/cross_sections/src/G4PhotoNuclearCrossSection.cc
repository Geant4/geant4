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
// * authors in the GEANT4 collaboration.                             *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4PhotoNuclearCrossSection.cc,v 1.2 2001-10-26 09:48:00 mkossov Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// G4 Physics class: PhotoNuclearCrossSection for gamma+A cross sections
// M.V. Kossov, ITEP(Moscow), 10-OCT-01
// 

#include "G4PhotoNuclearCrossSection.hh"

// The main member function giving the gamma-A cross section (E in GeV, CS in mb)
G4double G4PhotoNuclearCrossSection::GetCrossSection(const G4DynamicParticle* aPart,
													 const G4Element* anEle, G4double T)
{
  const G4double kinEnergy = aPart->GetKineticEnergy()/MeV;
  const G4double lE=log(kinEnergy);
  const G4int targetAtomicNumber = static_cast<int>(anEle->GetN()+.499); //@@ Nat mixture
  const G4int targZ = static_cast<int>(anEle->GetZ());
  const G4int targN = targetAtomicNumber-targZ;
  // Associative memory for acceleration
  static G4int    lastN;     // The last N of calculated nucleus
  static G4int    lastZ;     // The last Z of calculated nucleus
  static G4double lastHighE; // The last High Energy A-dependence
  static G4double lastGDRc1; // The last logAmplitude of the 1-st GDR maximum
  static G4double lastGDRp1; // The last A-power of the 1-st GDR maximum
  static G4double lastGDRt1; // The last Threshold of the 1-st GDR maximum
  static G4double lastGDRs1; // The last Slope of the 1-st GDR maximum
  static G4double lastGDRc2; // The last logAmplitude of the 2-nd GDR maximum
  static G4double lastGDRp2; // The last A-power of the 2-nd GDR maximum
  static G4double lastGDRt2; // The last Threshold of the 2-nd GDR maximum
  static G4double lastGDRs2; // The last Slope of the 2-nd GDR maximum
  static G4double lastQDAmp; // The last Amplitude of the QuasiDeuteron region [exp/(1+exp)]
  static G4double lastQDWid; // The last Width of the QuasiDeuteron region [.4] or SecRes (H1,H2)
  static G4double lastQDPos; // The last Position of the QuasiDeuteron region [3.8] or SecRs (H1,H2)
  static G4double lastDelAm; // The last Amplitude of the Delta Resonance [.41*(Z+N)]
  static G4double lastDelWd; // The last Width of the Delta Resonance [11.9-ln(A)*1.24]
  static G4double lastDelPs; // The last Position of the Delta Resonance [5.84-.09/(1+.003*A^2)]
  static G4double lastDelTh; // The last Threshold of the Delta Resonance [5.13-.00075*A]
  static G4double lastDelSl; // The last Slope of the Delta Resonance [0.04->0.09]
  static G4double lastRopAm; // The last Amplitude of the Roper Resonance [-2.+ln(A)*0.84]
  static G4double lastRopWd; // The last Width of the Roper Resonance [.1+1.65*ln(A)]
  static G4double lastRopPs; // The last Position of the Roper Resonance [6.46+.061*ln(A)]       
  static G4std::vector <G4int> colN;       // N of calculated nucleus
  static G4std::vector <G4int> colZ;       // Z of calculated nucleus
  static G4std::vector <G4double> HighE;   // High Energy A-dependence
  static G4std::vector <G4double> GDRc1;   // logAmplitude of the 1-st GDR maximum
  static G4std::vector <G4double> GDRp1;   // A-power of the 1-st GDR maximum
  static G4std::vector <G4double> GDRt1;   // Threshold of the 1-st GDR maximum
  static G4std::vector <G4double> GDRs1;   // Slope of the 1-st GDR maximum
  static G4std::vector <G4double> GDRc2;   // logAmplitude of the 2-nd GDR maximum
  static G4std::vector <G4double> GDRp2;   // A-power of the 2-nd GDR maximum
  static G4std::vector <G4double> GDRt2;   // Threshold of the 2-nd GDR maximum
  static G4std::vector <G4double> GDRs2;   // Slope of the 2-nd GDR maximum
  static G4std::vector <G4double> QDAmp;   // Amplitude of the QuasiDeuteron region [exp/(1+exp)]
  static G4std::vector <G4double> QDWid;   // Width of the QuasiDeuteron region [.4] or R(H1,H2)
  static G4std::vector <G4double> QDPos;   // Position of the QuasiDeuteron region [3.7] or R(H1,H2)
  static G4std::vector <G4double> DelAm;   // Amplitude of the Delta Resonance [.41*(Z+N)]
  static G4std::vector <G4double> DelWd;   // Width of the Delta Resonance [11.9-ln(A)*1.24]
  static G4std::vector <G4double> DelPs;   // Position of the Delta Resonance [5.84-.09/(1+.003*A2)]
  static G4std::vector <G4double> DelTh;   // Threshold of the Delta Resonance [5.13-.00075*A]
  static G4std::vector <G4double> DelSl;   // Slope of the Delta Resonance [0.04->0.09]
  static G4std::vector <G4double> RopAm;   // Amplitude of the Roper Resonance [-2.+ln(A)*0.84]
  static G4std::vector <G4double> RopWd;   // Width of the Roper Resonance [.1+1.65*ln(A)]
  static G4std::vector <G4double> RopPs;   // Position of the Roper Resonance [6.46+.061*ln(A)]
  G4double sigma=0.;
  if( aPart->GetDefinition()->GetPDGEncoding() == 22 &&
	  kinEnergy                                 > ThresholdEnergy(targZ, targN))
  {
    G4double A=targN+targZ;
    if(targN!=lastN || targZ!=lastZ)          // Otherwise the set of parameters is ready
	{
      G4int n=colN.size();
      G4bool in=false;
      if(n) for(G4int i=0; i<n; i++) if(colN[i]==targN && colZ[i]==targZ) // Calculated nucleus
	  { // @@ Parameters can be combined in a type (structure) to accelerate the retrieve process
        in=true;
        lastHighE=HighE[i];                   // High Energy A-dependence
        lastGDRc1=GDRc1[i];                   // logAmplitude of the 1-st GDR maximum
        lastGDRp1=GDRp1[i];                   // A-power of the 1-st GDR maximum
        lastGDRt1=GDRt1[i];                   // Threshold of the 1-st GDR maximum
        lastGDRs1=GDRs1[i];                   // Slope of the 1-st GDR maximum
        lastGDRc2=GDRc2[i];                   // logAmplitude of the 2-nd GDR maximum
        lastGDRp2=GDRp2[i];                   // A-power of the 2-nd GDR maximum
        lastGDRt2=GDRt2[i];                   // Threshold of the 2-nd GDR maximum
        lastGDRs2=GDRs2[i];                   // Slope of the 2-nd GDR maximum
        lastQDAmp=QDAmp[i];                   // Amplitude of the QuasiDeuteron region
        lastQDWid=QDWid[i];                   // Width of the QuasiDeuteron region or SecR(H1,H2)
        lastQDPos=QDPos[i];                   // Position of the QuasiDeuteron region or SecR(H1,H2)
        lastDelAm=DelAm[i];                   // Amplitude of the Delta Resonance
        lastDelWd=DelWd[i];                   // Width of the Delta Resonance
        lastDelPs=DelPs[i];                   // Position of the Delta Resonance
        lastDelTh=DelTh[i];                   // Threshold of the Delta Resonance
        lastDelSl=DelSl[i];                   // Slope of the Delta Resonance
        lastRopAm=RopAm[i];                   // Amplitude of the Roper Resonance
        lastRopWd=RopWd[i];                   // Width of the Roper Resonance
        lastRopPs=RopPs[i];                   // Position of the Roper Resonance       
	  }
	  if(!in)                                 // Fill the new set of parameters for the new nucleus
	  {
        lastN    = targN;                     // The last N of calculated nucleus
        lastZ    = targZ;                     // The last Z of calculated nucleus
        G4double lnA=log(A);
        if(A==1)
        {
          lastHighE=1.;                       // High Energy A-dependence
          lastDelAm=.55;                      // The last Amplitude of the Delta Resonance (not .41)
          lastQDAmp=.08;                      // The last Amplitude of the Third Resonance (like H2)
          lastQDWid=90.;                      // The last Width of the Third Resonance (like H2)
		  lastQDPos=6.93;                     // The last Position of the Third Resonance (like H2)
          lastDelWd=18.;                      // The last Width of the Delta Resonance (not 11.9)
          lastRopAm=.22;                      // The last Amplitude of the Roper Resonance (not .14)
          lastRopWd=20.;                      // The last Width of the Roper Resonance (close to H2)
          lastRopPs=6.57;                     // The last Position of the Roper Resonance (like H2)
		}
        else
        {
          lastHighE=exp(-lnA*(.115-.0048*lnA))*A; // High Energy A-dependence
          lastGDRc1=GetGDRc1(targZ, targN);   // The last logAmplitude of the 1-st GDR maximum
          lastGDRp1=GetGDRp1(targZ, targN);   // The last A-power of the 1-st GDR maximum
          lastGDRt1=GetGDRt1(targZ, targN);   // The last Threshold of the 1-st GDR maximum
          lastGDRs1=GetGDRs1(targZ, targN);   // The last Slope of the 1-st GDR maximum
          lastGDRc2=GetGDRc2(targZ, targN);   // The last logAmplitude of the 2-nd GDR maximum
          lastGDRp2=GetGDRp2(targZ, targN);   // The last A-power of the 2-nd GDR maximum
          lastGDRt2=GetGDRt2(targZ, targN);   // The last Threshold of the 2-nd GDR maximum
          lastGDRs2=GetGDRs2(targZ, targN);   // The last Slope of the 2-nd GDR maximum
          lastDelWd=GetDelWd(targZ, targN);   // The last Width of the Delta Resonance
          if(A==2)
		  {
            lastDelAm=.88;                    // The last Amplitude of the Delta Resonance (not .82)
            lastQDAmp=.078;                   // The last Amplitude of the Third Resonance
            lastQDWid=90.;                    // The last Width of the Third Resonance (like H1)
		    lastQDPos=6.93;                   // The last Position of the Third Resonance (like H1)
            lastRopAm=.34;                    // The last Amplitude of the Roper Resonance (not .14)
            lastRopWd=15.;                    // The last Width of the Roper Resonance (close to H1)
            lastRopPs=6.57;                   // The last Position of the Roper Resonance (like H1)
		  }
          else
          {
            lastDelAm=GetDelAm(targZ, targN); // The last Amplitude of the Delta Resonance
            lastQDAmp=GetQDAmp(targZ, targN); // The last Amplitude of the QuasiDeuteron region
            lastQDWid=.4;                     // The last Width of the QuasiDeuteron region
	    	lastQDPos=3.7;                    // The last Position of the QuasiDeuteron region
            lastRopAm=GetRopAm(targZ, targN); // The last Amplitude of the Roper Resonance
            lastRopWd=GetRopWd(targZ, targN); // The last Width of the Roper Resonance
            lastRopPs=GetRopPs(targZ, targN); // The last Position of the Roper Resonance       
		  }
		}
        lastDelPs=GetDelPs(targZ, targN);     // The last Position of the Delta Resonance
        lastDelTh=GetDelTh(targZ, targN);     // The last Threshold of the Delta Resonance
        lastDelSl=GetDelSl(targZ, targN);     // The last Slope of the Delta Resonance
        colN.push_back(targN);
        colZ.push_back(targZ);
        HighE.push_back(lastHighE);
        GDRc1.push_back(lastGDRc1);
        GDRp1.push_back(lastGDRp1);
        GDRt1.push_back(lastGDRt1);
        GDRs1.push_back(lastGDRs1);
        GDRc2.push_back(lastGDRc2);
        GDRp2.push_back(lastGDRp2);
        GDRt2.push_back(lastGDRt2);
        GDRs2.push_back(lastGDRs2);
        QDAmp.push_back(lastQDAmp);
        QDWid.push_back(lastQDWid);
        QDPos.push_back(lastQDPos);
        DelAm.push_back(lastDelAm);
        DelWd.push_back(lastDelWd);
        DelPs.push_back(lastDelPs);
        DelTh.push_back(lastDelTh);
        DelSl.push_back(lastDelSl);
        RopAm.push_back(lastRopAm);
        RopWd.push_back(lastRopWd);
        RopPs.push_back(lastRopPs);
	  } // End of creation of the new set of parameters
    } // End of parameters udate
    // ============================== NOW the Magic Formula =================================
    G4double qdeut=lE-lastQDPos;
    G4double delta=lE-lastDelPs;
    G4double roper=lE-lastRopPs;
    if(A>1) sigma+=exp(lastGDRc1-lE*lastGDRp1)/(1.+exp((lastGDRt1-lE)/lastGDRs1))+       // 1-st GDR
			       exp(lastGDRc2-lE*lastGDRp2)/(1.+exp((lastGDRt2-lE)/lastGDRs2));       // 2-nd GDR
	sigma+=lastQDAmp/(1.+lastQDWid*qdeut*qdeut) + lastRopAm/(1.+lastRopWd*roper*roper)+  // QD+Roper
	       lastDelAm/(1.+lastDelWd*delta*delta)/(1.+exp((lastDelTh-lE)/lastDelSl))+      // Delta
	       lastHighE*(0.0116*exp(lE*0.16)+.4*exp(-lE*0.2))/(1.+exp((7.-lE)/0.2));        // High E
  } // End of "sigma" calculation
  return sigma*millibarn;
}

// Correction function for Be,C @@ Move to header
G4double G4PhotoNuclearCrossSection::LinearFit(G4double X, G4int N, const G4double* XN,
											   const G4double* YN)
{
  G4double Xj=XN[0];
  G4double Xh=XN[N-1];
  if(X<=Xj) return Xj; //-----+
  else if(X>=Xh) return Xh;//-|
  G4double Xp=0.; //          |
  G4int j=0;   //             |
  while (X>Xj && j<N)//<------+
  {
    j++;
    Xp=Xj;
    Xj=XN[j];
  }
  return YN[j]-(Xj-X)*(YN[j]-YN[j-1])/(Xj-Xp);
}


