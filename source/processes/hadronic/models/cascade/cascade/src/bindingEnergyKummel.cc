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

G4double G4InuclSpecialFunctions::bindingEnergyKummel(G4double A, 
						      G4double Z) {
  G4int verboseLevel = 2;
  if (verboseLevel > 3) {
    G4cout << " >>> G4InuclSpecialFunctions::bindingEnergyKummel" << G4endl;
  }

  // calculates the nuclei binding energy using Kummel mass formula

// data for Kummel mass formula 
  const G4double ASH[6] = {20.0, 28.0, 50.0, 82.0, 126.0, 184.0};
  const G4double OMN[3] = {2.511464e0, 7.640294e-2, 2.618317e-2};
  const G4double OMP[3] = {2.511464e0, 7.640294e-2, 6.243892e-3};
  const G4double TMET[5] = {9.477193e0, 7.166809e-1, 1.2169e-1, 4.07179e-2, 2.927847e-2}; 
  const G4double TMET1[4] = {-4.775384e-2, 5.880883e-5, -1.98964e-3, 3.427961e-6}; 

  const G4double BET[5] = {1.535628e4, 1.502538e4, 1.481831e4, 1.459232e4, 1.40664e4};

  const G4double SHN[5] = {8.0, 22.0, 32.0, 44.0, 58.0};

  const G4double TNUN[5] = {8.728778e2, 8.45789e1, 9.210738e1, 5.98e1, 3.724e1}; 

  const G4double TNUP[5] = {8.728778e2, 8.45789e1, 3.74e1, 5.42e1, 3.724e1};

  const G4double TKSN[5] = {3.34582e1, 3.858216e0, -1.489218e-1, -9.2e-1, -6.36e-1};

  const G4double TKSP[5] = {3.34582e1, 3.858216e0, 7.2e-1, -8.2e-1, -6.36e-1};

  const G4double RON[5] = {-6.441735e3, -5.422058e3, -3.722645e3, -2.437172e4, -2.645867e4};

  const G4double ROP[5] = {-6.75601e3, -5.877588e3, -3.216382e4, -3.010688e4, -2.645867e4}; 

  const G4double SINGM[4] = {-3.578531e0, -2.6911e0, -7.487347e-1, 0.0};

  const G4double C[2] = {6.1e3, 8.31e3};

  const G4double AKU = 6.04122e5;

  const G4double US = 1.661835e4;

  const G4double UC = 6.218614e2;

  const G4double UT = 1.983151e4;

  const G4double SIPG = -2.067547e0;

  const G4double TAU = 2.231241e0;

  const G4double AL0 = 0.151;

  const G4double USB = 1.95114e4;

  const G4double UCB = 753.3;

  //  const G4double ALD = 15.4941;
  // const G4double ALD1 = 17.9439;

  // const G4double C3 = 0.7059;

  //const G4double C4 = 1.21129;

  //  const G4double PKLD = 1.7826;

  const G4double DMU = 3.132902e4;
        
  G4double DM;

  G4double AN = A - Z;

  G4int INS = 0;
  G4double ANE = 0.0;
  G4double HNE = 0.0;
  G4int i(0);

  for (i = 1; i < 6; i++) {

    if (AN <= ASH[i]) {
      ANE = AN - ASH[i - 1];
      HNE = ASH[i] - AN;
      INS = i - 1;

      break;
    };
  };
 
  G4int IPS = 0;

  G4double APR = 0.0;

  G4double HPR = 0.0;

  for (i = 1; i < 6; i++) {

    if (Z <= ASH[i]) {
      APR = Z - ASH[i - 1];
      HPR = ASH[i] - Z;
      IPS = i - 1;

      break;
    };
  };
 
  G4double PPHN = ANE * HNE;
  G4double PPHZ = APR * HPR;

  // omega terms
  G4double OMT = 0.0;

  if (INS <= 2) OMT += OMN[INS] * PPHN;

  if (IPS <= 2) OMT += OMP[IPS] * PPHZ;

  if (verboseLevel > 3) {
    G4cout << " OMT " << OMT << G4endl;
  }

  // theta term 1

  G4double TET = 0.0;

  if (verboseLevel> 3) {
    G4cout << " PPHN " << PPHN << " PPHZ " << PPHZ << G4endl;
    G4cout << " INS " << INS << " IPS " << IPS << G4endl;
  }

  if (PPHN > 0.5 && PPHZ > 0.5) {
    G4int IT = 0;

    if ((INS + 1) * (IPS + 1) > 0) {

      if (INS - IPS == 1) {

	if (INS < 5) {
	  IT = INS;

	} else {
	  IT = -1; 
	};		 

      } else {
	IT = -1;
      }; 
    };

    if (verboseLevel > 3){
      G4cout << " IT " << IT << G4endl; 
    }

    if (IT >= 0) TET = TMET[IT] * PPHN * PPHZ;
  };

  if (verboseLevel> 3) {
    G4cout << " TET " << TET << G4endl;
  }

  // theta term 2
  G4double TET1 = 0.0;

  if (IPS == 2) {
    TET1 += TMET1[0] * PPHZ * PPHZ;

    if (INS == 3) TET1 += TMET1[1] * PPHN * PPHZ * HNE * HPR;
  };

  if (INS == 3) {

    G4double TVSP = PPHN * PPHN * HNE;
    TET1 += TMET1[2] * TVSP;
    TET1 += TMET1[3] * TVSP * PPHN;
  };

  if (verboseLevel > 3) {
    G4cout << " TET1 " << TET1 << G4endl;
  }

  // betta, nu, ksy terms
  G4double TBET = 0.0;

  if (INS != 0) {

    for (G4int i = 0; i <= INS - 1; i++) TBET += BET[i] * SHN[i];
  };
  TBET += BET[INS] * ANE;
	
  if (IPS != 0) {

    for (G4int i = 0; i <= IPS - 1; i++) TBET += BET[i] * SHN[i];
  };
  TBET += BET[IPS] * APR;

  if (verboseLevel> 3) {
    G4cout << " TBET " << TBET << G4endl;
  }

  G4double TBET1 = 0.0;

  if (PPHN > 0.0 || PPHZ > 0.0) 
    TBET1 = 0.5 * ((TNUN[INS] + TKSN[INS] * ANE) * PPHN + 
		   (TNUP[IPS] + TKSP[IPS] * APR) * PPHZ);
  TBET -= TBET1;

  if (verboseLevel > 3) {
    G4cout << " TBET1 " << TBET1 << G4endl;
  }

  // deformation
  G4double TDEF = 0.0;
  G4double X = Z * Z / A;
  G4double X1 = pow(A, 0.3333333);
  G4double X2 = X1 * X1;

  if (IPS != INS && INS >= 3 && IPS >= 2) {
    G4double X3 = 2.0 * USB - UCB * X;
    G4double DNZ = 0.0;

    if (!(INS < 3 || (INS - IPS) != 1)) DNZ = C[INS - 3];
    DNZ = TBET1 - DNZ;
    G4double X4 = 0.2 * X2 * AL0 * AL0 * X3;
  
    if (DNZ > X4) {
      G4double X5 = USB + X * UCB;
      G4double X6 = log(DNZ / X4);
      G4double X7 = sqrt(X6);

      // G4double ALM = AL0 * (X7 + 0.143 * AL0 * X5 / X3);
      TDEF = -X4 * (X6 + 1.0) + DNZ + 0.038 * X2 * pow(AL0 * X7, 3.0) * X5;
    };
  };

  if (verboseLevel > 3) {
    G4cout << " TDEF " << TDEF << G4endl;
  }

  // pairing
  G4double TPE = 0.0; 
  G4double TPEN = 0.0;
  G4double TPEP = 0.0;
  G4double DN0 = 0.0;
  G4double DZ0 = 0.0;
  G4double AV = 2.0 * int(0.5 * AN + 0.1); 	

  if (AN > AV) DN0 = 1.0;

  AV = 2.0 * int(0.5 * Z + 0.1); 

  if (Z > AV) DZ0 = 1.0;

  if (DN0 * DZ0 > 0.0) TPE = DMU / A;

  if (DN0 + DZ0 > 0.0) {

    if (DN0 > 0.0) {
      TPEN = RON[INS] / X1;

      if (INS >= 1) {

	if (AN >= 82.0) {
	  TPEN /= X1;

	  if (INS == 3) {
	    G4double PN = 0.0;
	    G4double HN = 0.0;

	    if (AN > 90.0) PN = AN - 90.0;

	    if (AN < 116.0) HN = 116.0 - AN;
	    TPEN += TAU * PN * HN;
	  };
	  TPEN += SINGM[INS - 1] * PPHN;
	}; 
      };
    };

    if (DZ0 > 0.0) {
      TPEP = ROP[IPS] / X1;

      if (IPS == 1) TPEP += SIPG * PPHZ;

      if (IPS > 1) TPEP = TPEP / X1;

    };  
  };
  TPE += TPEN + TPEP;

  if (verboseLevel > 3) {
    G4cout << " TPE " << TPE << " TPEN " << TPEN << " TPEP " << TPEP << G4endl;
  }

  // collect everything
  DM = (AKU - US * X2 + Z * (Z - 1.0) / X1 * (OMT - UC) - UT / A *
	(AN - Z) * (AN - Z) + TET + TET1 + TBET + TDEF + TPE) * 0.001;	  

  if (verboseLevel > 3) {
    G4cout << " kummel " << G4endl; 
  }

  return DM;
}
