#include "globals.hh"
#include "Randomize.hh" // AH for test environment only: delete in production version

#include "G4Cascade.hh"

#include "g4std/iostream"

G4Cascade::G4Cascade():verboseLevel(0){};
G4Cascade::~G4Cascade(){};

G4int G4Cascade::go()
{   
if( verboseLevel > 1 ) G4cout << " Entering go" << G4endl;
   rout1(); // AH Testing
if( verboseLevel > 1 ) G4cout << " Leaving go" << G4endl;
return 0;
}
void stor( G4float p0,   // AH implemetn this
           G4int i, 
           G4float e, 
           G4float a, 
           G4float b, 
           G4float g,
           G4float w, 
           G4float j, 
           G4float erem, 
           G4float p, 
           G4float wmass, 
           G4float itype ){}; 


 void G4Cascade::pol1( G4double &cs, G4double &si )
  {
    G4double u;
    cs = dflran();
    u = dflran();
    if( u < 0.5 ) cs = -cs;
    si = sqrt( 1.0-(cs*cs) );
    return;
  }

 void G4Cascade::azio( G4double &s, G4double &c )
  {
    G4double r1, r2, r1sq, r2sq;
    G4double sum = 2.0;
    while( sum > 1.0 )
    {
      r1 = dflran();
      r2 = dflran();
      sum = r1*r1 + r2*r2;    // ysq
    }
    sum /= 2.0;               // (xsq+ysq)/2
    c = (sum-r1*r1)/sum;      // (ysq-xsq)/(xsq+ysq)
    s = (r1*r2)/sum;          // (2*x*y)/(xsq+ysq)
    r1 = dflran();
    if( r1 >= 0.5 )s = -s;
    return;
  }

 void G4Cascade::piAH( G4float po, 
		     G4int itip, 
		     G4float ener, 
		     G4float pmom, 
		     G4float a, 
		     G4float b, 
		     G4float g )
  {
    G4float ratt[12];
    G4float dum[18];
    G4double r2, csp, snp;
    const G4float ang[15] = {0.0, 0.0, 0.01745, 0.0349, 0.0698, 0.1396, 0.2094, 0.3491, 0.7845,
                             0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    const G4float pm[] = {4.5, 7.5, 10.5, 14.0, 22.0, 34.0, 50.0, 90.0, 170.0, 290.0, 490.0, 1000.0};
    const G4float dndpip[12][3] = {{-5.3242445e-03,    2.2346497e-02,    4.5267868e-01},
                                   {1.1044145e-03,   -3.2094955e-02,    5.6748962e-01},
                                   {8.9269876e-04,   -2.8545380e-02,    5.5277443e-01},
                                   {5.0222874e-04,   -2.0401001e-02,    5.1032925e-01},
                                   {2.0570680e-04,   -1.1921406e-02,    4.4970989e-01},
                                   {6.6768611e-05,   -5.8776736e-03,    3.8399315e-01},
                                   {2.2491440e-05,   -2.9177666e-03,    3.3454037e-01},
                                   {5.6584831e-06,   -1.2033470e-03,    2.9090124e-01},
                                   {9.4826100e-07,   -3.7099421e-04,    2.5414348e-01},
                                   {1.6872946e-07,   -1.1652336e-04,    2.3341179e-01},
                                   {3.4642653e-08,   -4.0553510e-05,    2.2265625e-01},
                                   {5.3050826e-09,   -1.1835014e-05,    2.1562809e-01}};
    const G4float dndpim[12][3] = {{-3.7355423e-03,    2.9491067e-02,    2.0408440e-01},
                                   {-1.3875589e-04,   -1.4313459e-03,    2.7039909e-01},
                                   {1.4499575e-04,   -5.3753257e-03,    2.8402138e-01},
                                   {1.0721385e-04,   -4.5639873e-03,    2.7966785e-01},
                                   {4.8756599e-05,   -2.8809309e-03,    2.6756477e-01},
                                   {1.6727019e-05,   -1.4842749e-03,    2.5234318e-01},
                                   {5.7825819e-06,   -7.5179338e-04,    2.4009037e-01},
                                   {1.4668331e-06,   -3.1205267e-04,    2.2888833e-01},
                                   {2.4605833e-07,   -9.6268952e-05,    2.1935964e-01},
                                   {4.3786713e-08,   -3.0238181e-05,    2.1398067e-01},
                                   {8.9794412e-09,   -1.0516495e-05,    2.1118736e-01},
                                   {1.3733370e-09,   -3.0666124e-06,    2.0936286e-01}};
    G4float a1 = 3.386;
    G4float a2 = 4.146;
    G4float a3 = 4.556;
    G4float ta3 = a3+a3;
    G4float a5 = 9.600;
    G4float a6 = 4.823;
    G4float topi = 6.28318;
    G4float xp   = 0.7853982;
    G4float xp2a3 = xp*xp*a3;
    G4float sxp  = sqrt(xp);
    G4float spo  = sqrt(po);
    G4float pos  = po*po;
    G4float r1= dflran();
G4float dmax;
G4float a4;
G4float rm;
    switch ( itip )
    {
     case 1:  // proton
     case 2:  // neutron
      a4 = 7.141;
      rm = 0.940;
      break;
     case 3:  // pi plus
     case 4:  // pi zero
       a4 = 7.141;
      rm = 0.139;
      break;
     case 5:  // pi minus
      a4 = 1.853;
       rm = 0.139;
    }
    //
    // sample uniformly from po
    //
    G4float ps = po*r1;
    //
    // find dndp(max) for ps
    //
G4int inx;
    if( po < pm[6] )
    {
      if( po < pm[4] )
      {
        if( po > pm[2] )
        {
          inx = 3;
          if( po > pm[3] )inx = 4;
        }
        else
        {
          inx = 1;
          if( po > pm[1] )inx = 2;
        }
        goto L50;
      }
      inx = 5;
      if( po > pm[5] )inx = 6;
    }
    else if( po < pm[9] )
    {
      inx = 7;
      if( po >= pm[7] )
      {
        inx = 8;
        if( po > pm[8] )inx = 9;
      }
    }
    else
    {
      inx = 10;
      if( po >= pm[10] )
      {
        inx = 11;
        if( po > pm[11] )inx = 12;
      }
    }
    goto L50;
  L50:
    //    goto (51, 51, 51, 51, 52), itip; // AH use in code
  L51:
    dmax = dndpip[1][inx]*po*po + dndpip[2][inx]*po + dndpip[3][inx];
    goto L53;
  L52:
    dmax = dndpim[1][inx]*po*po + dndpim[2][inx]*po + dndpim[3][inx];
  L53:
    if( ener != 0.0 )ps = pmom;
G4float pss = ps*ps;
G4float psa6 = ps*a6;
G4float sps = sqrt(ps);
    G4float f1 = a1*pss*exp(-a2*ps/spo);
    G4float f2 = (a4*pss/po)*exp(-a5*pss/pos);
    G4float d1 = ta3*ps*spo;
    G4float d2 = psa6*psa6;
    G4float e1 = xp2a3*ps*spo;
    if( e1 > 50.0 )e1 = 50.0;
    G4float e2 = xp*psa6;
    if( e2 > 50.0 )e2 = 50.0;
    G4float dpsmx = topi*f1/d1*(1.0-exp(-e1)) + topi*f2/d2*(1.0-exp(-e2)*100.0);
    G4float ratio = dpsmx/dmax;
    r2 = dflran();
    if( ener != 0.0 )goto L8;
    if( r2 <= ratio )goto L8;
    r1 = dflran();
    ps =  po*r1;
    goto L53;
  L8:
    //
    // for current ps and r2,  find angles
    //
    // set angle intervals - coarse-  test for ratio and r2
    //  0 -  1 degrees (.00000 - .01745 radians)
    //  1 -  2 degrees (.01745 - .03490 radians)
    //  2 -  4 degrees (.03490 - .06980 radians)
    //  4 -  8 degrees (.06980 - .13960 radians)
    //  8 - 12 degrees (.13960 - .20940 radians)
    // 12 - 20 degrees (.20940 - .34910 radians)
    // 20 - 45 degrees (.34910 - .78540 radians)
    G4float f3 = topi*f1;
    G4float f4 = topi*f2;
    G4float a7 = 1.0/d1;
    G4float a8 = 1.0/d2;
    G4float f3a7 = f3*a7;
    G4float f4a8 = f4*a8;
    a7 *= 2;
    r2 = dflran();
    G4float r3 = r2;
    G4float exl = -1.0;
    G4float f6l = -1.0;
    ratt[1] = 0.0;
    G4float angfrs = 0.0;
    G4float anglas = ang[9];
  L58:
    G4float tu = (angfrs+anglas)*0.5;
    G4float txu = tu*tu/a7;
    if( txu > 50.0 )txu = 50.0;
    G4float exu = -exp(-txu);
    G4float f6u = 0.0;
    if( psa6*tu > 50.0 )goto L59;
    f6u = exp(-psa6*tu)*(-psa6*tu-1.0);
  L59:
    G4float f5 = exu-exl;
    G4float f6 = f6u-f6l;
    G4float rat = f3a7*f5+f4a8*f6;
    ratt[2] = rat/dpsmx;
    if( abs(anglas-tu)/anglas < 1.0e-4 )goto L80;
    if( r3 <= ratt[2] )goto L68;
    angfrs = tu;
    goto L58;
  L68:
    anglas = tu;
    goto L58;
  L80:
    G4float tht = tu;
    //
    // find angles -
    //
    g = cos(tht);
    azio( csp, snp );
    G4float z = sqrt(1.0-g*g);
    a = z * csp;
    b = z * snp;
    //
    ener = sqrt(pss+rm*rm) - rm;
    if( ener < 0.0 )ener = 1.0e-6;
    pmom = ps;
    return;
  }

 void G4Cascade::nucnuc( G4float p0, 
			 G4int nofas,
			 G4int itype )
  {
    // no nucleon conservation
    G4float dum[18];
    G4double val, csa, sna, csp, snp;
    G4float sp[5];
    G4float count = 0.0;
    G4float witty = 0.03;
    G4float ddd = 0.15;
    G4float c45d = 0.707;
    G4int ntimes = 0;
    G4float erem = einc;
    st0r = einc - 0.139;
    nofas = 0;
    prob( itype, p0 );
    sp[1] = pb[1]/(pb[1]+pb[2]);
    sp[2] = 1.0;
    sp[3] = pb[3]/(pb[3]+pb[4]+pb[5]);
    sp[4] = pb[4]/(pb[3]+pb[4]+pb[5])+sp[3];
    sp[5] = 1.0;
    numnuc = 0;
  L8:
    val = dflran();
    G4bool test = true;
    G4int i;
    for ( i = 1; i < 6; ++i )
    {
      if( val < pb[i] )
      {
        test = false;
        break;
      }
    }
    if( test ) G4Exception("nucnuc2");
    G4int j = 1;
    G4int e;
    G4int a=0; // AH dummy value added
    G4int b=0; // AH dummy value added
    G4int g=0; // AH dummy value added
    G4int w;
    G4int p=0; // AH dummy value added
    G4int wmass;
    stor( p0, i, e, a, b, g, w, j, erem, p, wmass, itype );
    if( w != 0.0 )
    {
      erem -= e+wmass;
      if( erem < 0.0 )goto L10;
  L23:
      if( ++nofas > 60 )
      {
        //write(io, 8890);
        //8890 format('0', 10x , 'nofas greater than 60 in nucnuc');
        if( ++ntimes >= 10 )G4Exception("nucnuc1");
        erem = einc;
        st0r = einc - 0.139;
        nofas = 0;
        prob( itype, p0 );
        sp[1] = pb[1]/(pb[1]+pb[2]);
        sp[2] = 1.0;
        sp[3] = pb[3]/(pb[3]+pb[4]+pb[5]);
        sp[4] = pb[4]/(pb[3]+pb[4]+pb[5])+sp[3];
        sp[5] = 1.0;
        numnuc = 0;
      }
      else
      {
        if( i > 2 || ++numnuc < 3 )
        {
          efas[nofas] = e*1000.0;
          it[nofas] = i-1;
          alpfas[nofas] = a;
          betfas[nofas] = b;
          gamfas[nofas] = g;
          wtfas[nofas] = abs(w);
          if( w < 0.0 )return;
        }
        else
        {
          --nofas;
          j = 2;
          stor( p0, i, e, a, b, g, w, j, erem, p, wmass, itype );
          erem += e + wmass;
        }
      }
      goto L8;
    }
    w = 1.0;
    test = false;
   do
    {
      e = 0.0;
      if( i <= 2 )
      {
        prot( p0, e, p, a, b, g );
        if( e > st0r )test = true;
      }
      else
      {
        piAH( p0, i, e, p, a, b, g );
        if( e > st0r )test = true;
      }
    }
    while ( test );
    if( i > 2 )
    {
      j = 4;
      stor( p0, i, e, a, b, g, w, j, erem, p, wmass, itype );
      if( i == 0 )goto L8;
      wmass = 0.139;
      goto L5;
    }
    j = 4;
    stor( p0, i, e, a, b, g, w, j, erem, p, wmass, itype );
    if( i != 0 )
    {
      wmass = 0.0;
  L5:
      erem -= e+wmass;
      if( erem < 0.0 )goto L10;
      if( ++nofas > 60 )
      {
        //write(io, 8890);
        //8890 format('0', 10x , 'nofas greater than 60 in nucnuc');
        if( ++ntimes >= 10 )G4Exception("nucnuc1");
        erem = einc;
        st0r = einc - 0.139;
        nofas = 0;
        prob( itype, p0 );
        sp[1] = pb[1]/(pb[1]+pb[2]);
        sp[2] = 1.0;
        sp[3] = pb[3]/(pb[3]+pb[4]+pb[5]);
        sp[4] = pb[4]/(pb[3]+pb[4]+pb[5])+sp[3];
        sp[5] = 1.0;
        numnuc = 0;
      }
      else
      {
        if( i > 2 || ++numnuc < 3 )
        {
          efas[nofas] = e*1000.0;
          it[nofas] = i-1;
          alpfas[nofas] = a;
          betfas[nofas] = b;
          gamfas[nofas] = g;
          wtfas[nofas] = abs(w);
          if( w < 0.0 )return;
        }
        else
        {
          --nofas;
          j = 2;
          stor( p0, i, e, a, b, g, w, j, erem, p, wmass, itype );
          erem += e + wmass;
        }
      }
    }
    goto L8;
  L10:
    erem += e + wmass;
    j = 2;
    stor( p0, i, e, a, b, g, w, j, erem, p, wmass, itype );
    if( numnuc > 1 )
    {
      val  = dflran();
      i = 3;
      if( val > sp[3] )i = 4;
      if( val > sp[4] )i = 5;
      if( erem <= einc*ddd )
      {
        val = dflran();
        if( val < witty )
        {
          if( (i > 2) && (erem < 0.139) )
          {
            if( erem < 0.0 )G4Exception("nucnuc4");
            ++count;
            for ( G4int i = 1; i <= nofas; ++i )efas[i] += erem/nofas*1000.0;;
            return;
          }
          azio( csa, sna );
          do
          {
            pol1( csp, snp );
          }
          while( csp > c45d );
          a = snp*csa;
          b = snp*sna;
          g = csp;
          i > 2.0 ? wmass  = 0.139 : wmass = 0.0;
          e = erem - wmass;
          w = -1.0;
          goto L23;
        }
      }
      if( i > 2 )
      {
        e = erem - 0.139;
        if( e < 0.0 )
        {
          if( erem < 0.0 )G4Exception("nucnuc4");
          ++count;
          for ( G4int i = 1; i <= nofas; ++i )efas[i] += erem/nofas*1000.0;;
          return;
        }
        p = sqrt(e*e+0.278*e);
        w = 1.0;
        piAH( p0, i, e, p, a, b, g );
        j = 3;
        stor( p0, i, e, a, b, g, w, j, erem, p, wmass, itype );
        w = -1.0;
        goto L23;
      }
    }
    else if( numnuc > 0 )
    {
      val = dflran();
      i = 1;
      if( val > sp[1] )i = 2;
    }
    else
    {
      ++numnuc;
      val = dflran();
      val > sp[1] ? i = 2 : i = 1;
      val = dflran();
      e = erem*val;
      erem -= e;
      p = sqrt(e*(e+1.88));
      w = 1.0;
      prot( p0, e, p, a, b, g );
      j = 3;
      stor( p0, i, e, a, b, g, w, j, erem, p, wmass, itype );
      if( ++nofas > 60 )
      {
        //write(io, 8890);
        //8890 format('0', 10x , 'nofas greater than 60 in nucnuc');
        if( ++ntimes >= 10 )G4Exception("nucnuc1");
        erem = einc;
        st0r = einc - 0.139;
        nofas = 0;
        prob( itype, p0 );
        sp[1] = pb[1]/(pb[1]+pb[2]);
        sp[2] = 1.0;
        sp[3] = pb[3]/(pb[3]+pb[4]+pb[5]);
        sp[4] = pb[4]/(pb[3]+pb[4]+pb[5])+sp[3];
        sp[5] = 1.0;
        numnuc = 0;
        goto L8;
      }
      efas[nofas] = e*1000.0;
      it[nofas] = i-1;
      alpfas[nofas] = a;
      betfas[nofas] = b;
      gamfas[nofas] = g;
      wtfas[nofas] = abs(w);
      val = dflran();
      i = 1;
      if( val > sp[1] )i = 2;
    }
    e = erem;
    if( e < 0.0 )G4Exception("nucnuc3");
    p = sqrt(e*e+1.88*e);
    w = 1.0;
    prot( p0, e, p, a, b, g );
    j = 3;
    stor( p0, i, e, a, b, g, w, j, erem, p, wmass, itype );
    w = -1.0;
    goto L23;
  }
 
 
 void G4Cascade::pinuc( G4float p0, 
			G4int nofas, 
			G4int itype )
  {
    G4float dum[18];
    G4double val, csa, sna, csp, snp;
    G4float sp[5];
    G4float count = 0.0, witty = 0.03, ddd = 0.15, c45d = 0.707;
    G4int ntimes = 0;
    goto L8889;
  L8888:
    //write(io, 8890);
    //8890 format('0', 10x, 'nofas greater than 60 in pinuc');
    if( ++ntimes >= 10 )G4Exception("pinuc1");
  L8889:
    G4float erem = einc + 0.139;
    st0r = einc - 0.139;
    nofas = 0;
    prob( itype, p0 );
    G4float cons = pb[1]+pb[2];
    sp[1] = pb[1]/cons;
    sp[2] = 1.0;
    cons = pb[3]+pb[4]+pb[5];
    sp[3] = pb[3]/cons;
    sp[4] = pb[4]/cons+sp[3];
    sp[5] = 1.0;
    numnuc = 0;
  L8:
    val = dflran();
    for ( G4int i = 1; i < 6; i++ )
    {
      if( val < pb[i] )goto L2;
    }
    G4Exception("pinuc2");
  L2:
    G4int j = 1;
G4float e;
G4float a;
G4float b;
G4float g;
G4float w;
G4float p;
G4float wmass;
    stor( p0, i, e, a, b, g, w, j, erem, p, wmass, itype );
    if( w != 0.0 )goto L5;
    w = 1.0;
    wmass = 0.0;
    if( i > 2 )wmass = 0.139;
  L101:
    e = 0.0;
    piAH( p0, i, e, p, a, b, g );
    if( e > st0r )goto L101;
    j = 4;
    stor( p0, i, e, a, b, g, w, j, erem, p, wmass, itype );
    if( i == 0 )goto L8;
  L5:
    erem = erem - e - wmass;
    if( erem < 0.0 )goto L10;
  L23:
    nofas = nofas + 1;
    if( nofas > 60 )goto L8888;
    if( i > 2 )goto L24;
    numnuc = numnuc+1;
    if( numnuc < 2 )goto L24;
    nofas = nofas-1;
    j = 2;
    stor( p0, i, e, a, b, g, w, j, erem, p, wmass, itype );
    erem = erem + e + wmass;
    goto L8;
  L24:
    efas[nofas] = e*1000.0;
    it[nofas] = i-1;
    alpfas[nofas] = a;
    betfas[nofas] = b;
    gamfas[nofas] = g;
    wtfas[nofas] = abs(w);
    if( w < 0.0 )return;
    goto L8;
  L10:
    erem = erem + e + wmass;
    j = 2;
    stor( p0, i, e, a, b, g, w, j, erem, p, wmass, itype );
    if( numnuc > 0 )goto L330;
    val = dflran();
    i = 1;
    if( val > sp[1] )i = 2;
    e = erem;
    p = sqrt(e*(e+1.88));
    goto L663;
  L330:
    val = dflran();
    i = 3;
    if( val > sp[3] )i = 4;
    if( val > sp[4] )i = 5;
    if( erem > (einc*ddd) )goto L162;
    val = dflran();
    if( val < witty ) goto L200;
  L162:
    if( erem < 0.139 && i > 2 )goto L25;
    wmass = 0.0;
    if( i > 2 )wmass = 0.139;
    e = erem - wmass;
    if( e < 0.0 )G4Exception("wmass");
    if( i > 2 )goto L662;
    p = sqrt(e*e + 1.88*e);
    goto L663;
  L662:
    p = sqrt(e*e + 0.278*e);
  L663:
    w = 1.0;
    piAH(p0, i, e, p, a, b, g);
  L105:
    j = 3;
    stor( p0, i, e, a, b, g, w, j, erem, p, wmass, itype );
    w =  -1.0;
    goto L23;
  L25:
    if( erem < 0.0 )G4Exception("pinuc3");
    count = count + 1.0;
    { 
    G4float temp = (erem/nofas)*1.e3;
    for ( i = 1; i <= nofas; i++ ) efas[i] += temp;
    }

    return;
  L200:
    if( i > 2 && erem < 0.139 )goto L25;
    azio( csa, sna );
  L201:
    pol1( csp, snp );
    if( csp > c45d )goto L201;
    a = snp*csa;
    b = snp*sna;
    g = csp;
    wmass = 0.0;
    if( i > 2 )wmass = 0.139;
    e = erem - wmass;
    w = -1.0;
    goto L23;
  }
 
 void G4Cascade::prob( G4int itype, G4float pinc )
  {
    const G4float pz[12] = {4.5, 7.5, 10.5, 14.0, 22.0, 34.0, 50.0, 90.0, 170.0, 290.0, 490.0, 1000.0};
    const G4float sumas = 1.079965;
    const G4float tpmas = 1.880150;
    const G4float cpn = 1.75;
    const G4float cpp = 0.875;
    const G4float witty = 2.0;  // witty used in p + p to get nu, n
    const G4float cnunn = 2.0;  // cnunp and cnunn used in n + p to get nu, p and nu, n
    const G4float cnunp = 1.0;
    const G4float abc = 0.5;
    const G4float witt2 = 1.0;  // witt2, abc, and ab2 used in pi+ + p to get nu, p, nu, n, and nu, pi+
    const G4float ab2 = 1.0;
    const G4float aaa = 0.5;
    const G4float wit = 1.0;    // wit,  aaa,  abk used in pi- + p to get nu, p,  nu, n,  and nu, pi-
    const G4float abk = 1.0;
    const G4float cc3 = 0.0;    // cc3 used in pi+ + p to get nu, pi-
    const G4float cc4 = 0.0;    // cc4 used in pi- + p to get nu, pi+
G4float cnu[5];
G4float et[5];
    const G4float cpnu[12][3] = {{-1.3039768e-02,    1.3335037e-01,    6.7083645e-01},
                          {-2.2712350e-03,    4.0141106e-02,    8.7221718e-01}, 
                          {-6.1076880e-04,    1.6270638e-02,    9.5784855e-01}, 
                          {-2.3752451e-04,    8.5802078e-03,    9.9742126e-01}, 
                          {-7.3343515e-05,    3.9319992e-03,    1.0303373    }, 
                          {-1.8943101e-05,    1.5830994e-03,    1.0556831    }, 
                          {-5.4761767e-06,    6.8706274e-04,    1.0705557    }, 
                          {-1.2256205e-06,    2.5695562e-04,    1.0814734    }, 
                          {-1.9255094e-07,    7.5042248e-05,    1.0894747    }, 
                          {-3.3469405e-08,    2.3089349e-05,    1.0937023    }, 
                          {-6.9267116e-09,    8.0280006e-06,    1.0958271    }, 
                          {-1.0377335e-09,    2.3245811e-06,    1.0972347    }};
    const G4float cpk[12][3] = {{-1.4113426e-02,    5.3243160e-01,    3.4706497e-01},
                         {-3.8183928e-03,    4.4261360e-01,    5.4272461e-01}, 
                         {-1.6031265e-03,    4.1049194e-01,    6.5904236e-01}, 
                         {-8.7887049e-04,    3.9549255e-01,    7.3689270e-01}, 
                         {-4.1687489e-04,    3.8226700e-01,    8.3128357e-01}, 
                         {-1.7768145e-04,    3.7179470e-01,    9.4566345e-01}, 
                         {-8.0645084e-05,    3.6527634e-01,    1.0555573    }, 
                         {-3.0577183e-05,    3.6011600e-01,    1.1883698    }, 
                         {-9.1493130e-06,    3.5625267e-01,    1.3625183    }, 
                         {-2.9318035e-06,    3.5418224e-01,    1.5346680    }, 
                         {-1.0281801e-06,    3.5308743e-01,    1.6921387    }, 
                         {-2.9010698e-07,    3.5234547e-01,    1.8784180    }};
    const G4float cpipnu[12][3] = {{-1.3925135e-02,    1.9439220e-01,    2.6330948e-02},
                            {-4.9353242e-03,    1.1407471e-01,    2.0571613e-01}, 
                            {-1.9524097e-03,    7.0723534e-02,    3.6305332e-01}, 
                            {-9.7447634e-04,    5.0484657e-02,    4.6774292e-01}, 
                            {-4.1649491e-04,    3.4583092e-02,    5.8101559e-01}, 
                            {-1.6780943e-04,    2.3748994e-02,    6.9900227e-01}, 
                            {-7.8253448e-05,    1.7732978e-02,    7.9998779e-01}, 
                            {-3.3177203e-05,    1.3085365e-02,    9.1971684e-01}, 
                            {-1.2346776e-05,    9.3094110e-03,    1.0908327    }, 
                            {-5.0680246e-06,    6.8699718e-03,    1.2951612    }, 
                            {-2.2703025e-06,    5.2463412e-03,    1.5307465    }, 
                            {-8.7339140e-07,    3.8301349e-03,    1.8892488    }};
    const G4float cpipk[12][3] = {{-5.6201220e-04,    1.8275261e-01,   -1.3130379e-01},
                           {-9.0676546e-04,    1.8500519e-01,   -1.3445663e-01}, 
                           {-4.1097403e-04,    1.7770004e-01,   -1.0755920e-01}, 
                           {-1.9192696e-04,    1.7316246e-01,   -8.4091187e-02}, 
                           {-6.6518784e-05,    1.6959763e-01,   -5.8746338e-02}, 
                           {-1.8179417e-05,    1.6750813e-01,   -3.6132813e-02}, 
                           {-5.3644180e-06,    1.6665268e-01,   -2.1896362e-02}, 
                           {-1.1995435e-06,    1.6622925e-01,   -1.1184692e-02}, 
                           {-1.8626451e-07,    1.6604996e-01,   -3.3721924e-03}, 
                           {-2.9802322e-08,    1.6599941e-01,    4.8828125e-04}, 
                           {-3.7252903e-09,    1.6598701e-01,    2.6855469e-03}, 
                           {-1.1641532e-09,    1.6598129e-01,    3.6621094e-03}};
    const G4float cpimnu[12][3] = {{-6.8021417e-03,    1.0507488e-01,    2.5606155e-02},
                            {-2.4644732e-03,    6.6451073e-02,    1.1157703e-01}, 
                            {-1.0727644e-03,    4.6208382e-02,    1.8512058e-01}, 
                            {-5.9020519e-04,    3.6210060e-02,    2.3689651e-01}, 
                            {-2.8741732e-04,    2.7553260e-02,    2.9872513e-01}, 
                            {-1.3258681e-04,    2.0779014e-02,    3.7283039e-01}, 
                            {-6.7789108e-05,    1.6413033e-02,    4.4637489e-01}, 
                            {-3.0794181e-05,    1.2584507e-02,    5.4531479e-01}, 
                            {-1.1968659e-05,    9.1618299e-03,    7.0087528e-01}, 
                            {-5.0016679e-06,    6.8241954e-03,    8.9691067e-01}, 
                            {-2.2565946e-06,    5.2303076e-03,    1.1282806    }, 
                            {-8.7132503e-07,    3.8254932e-03,    1.4840240    }};
    const G4float cpimk[12][3] = {{-9.9664927e-04,    1.0531235e-01,   -6.1017036e-02},
                           {-4.5859814e-04,    1.0030270e-01,   -4.9365997e-02}, 
                           {-1.6081333e-04,    9.5970154e-02,   -3.3605576e-02}, 
                           {-6.6816807e-05,    9.4027519e-02,   -2.3590088e-02}, 
                           {-2.1457672e-05,    9.2741013e-02,   -1.4461517e-02}, 
                           {-5.5730343e-06,    9.2055321e-02,   -7.0800781e-03}, 
                           {-1.6018748e-06,    9.1790199e-02,   -2.6702881e-03}, 
                           {-3.5762787e-07,    9.1663361e-02,    5.1879883e-04}, 
                           {-5.2154064e-08,    9.1610909e-02,    2.8533936e-03}, 
                           {-5.8207661e-09,    9.1595650e-02,    4.2114258e-03}, 
                           {-1.8626451e-09,    9.1592789e-02,    4.7149658e-03}, 
                           {-4.6566129e-10,    9.1590405e-02,    4.7454834e-03}};
    G4float hold;
    G4float fectr;
    G4float factr;

    // initialize pb array
    //
    for ( G4int i = 1; i < 6; i++ ) pb[i] = 0.0;
    //
    // find pinc position in pz array to get all coefficients
    //
    G4int inx;
    if( pinc < pz[6] )goto L35;
    if( pinc < pz[9] )goto L30;
    inx = 10;
    if( pinc < pz[10] )goto L50;
    inx = 11;
    if( pinc > pz[11] )inx = 12;
    goto L50;
  L30:
    inx = 7;
    if( pinc < pz[7] )goto L50;
    inx = 8;
    if( pinc > pz[8] )inx = 9;
    goto L50;
  L35:
    if( pinc < pz[4] )goto L40;
    inx = 5;
    if( pinc > pz[5] )inx = 6;
    goto L50;
  L40:
    if( pinc > pz[2] )goto L45;
    inx = 1;
    if( pinc > pz[1] )inx = 2;
    goto L50;
  L45:
    inx = 3;
    if( pinc > pz[3] )inx = 4;
  L50:
    G4float pinc2 = pinc*pinc;
    // goto (100, 200, 300, 400, 500), itype; // AH add functionality
    // p + p collisions;
  L100:
    G4float etot = einc+tpmas;
    cnu[1] = pinc2*cpnu[1][inx]+pinc*cpnu[2][inx]+cpnu[3][inx];
    cnu[3] = pinc2*cpipnu[1][inx]+pinc*cpipnu[2][inx]+cpipnu[3][inx];
    cnu[5] = pinc2*cpimnu[1][inx]+pinc*cpimnu[2][inx]+cpimnu[3][inx];
    et[1] = pinc2*cpk[1][inx]+pinc*cpk[2][inx]+cpk[3][inx];
    et[3] = pinc2*cpipk[1][inx]+pinc*cpipk[2][inx]+cpipk[3][inx];
    et[5] = pinc2*cpimk[1][inx]+pinc*cpimk[2][inx]+cpimk[3][inx];
    G4float dnu1 = cnu[1];
    G4float dnu2 = witty-cnu[1];
    G4float et2 = et[1]*dnu2/cnu[1];
    G4float et1 = et[1];
    G4float dene;
    G4float fune;

    if( pinc < 18.0 )goto L220;
    cnu[2] = dnu2;
    et[2] = et2;
    fune = etot-et[1]-et[2]-et[3]-et[5];
    dene = et[3]/cnu[3];
    goto L600;
    //
    // n + p collisions
    //
  L200:
    etot = einc+tpmas;
    cnu[1] = pinc2*cpnu[1][inx]+pinc*cpnu[2][inx]+cpnu[3][inx];
    cnu[3] = pinc2*cpipnu[1][inx]+pinc*cpipnu[2][inx]+cpipnu[3][inx];
    cnu[5] = pinc2*cpimnu[1][inx]+pinc*cpimnu[2][inx]+cpimnu[3][inx];
    et[1] = pinc2*cpk[1][inx]+pinc*cpk[2][inx]+cpk[3][inx];
    et[3] = pinc2*cpipk[1][inx]+pinc*cpipk[2][inx]+cpipk[3][inx];
    et[5] = pinc2*cpimk[1][inx]+pinc*cpimk[2][inx]+cpimk[3][inx];
    dnu2 = cnunn-cnunp;
    et2 = et[1]*dnu2/cnu[1];
    et1 = et[1]*cnunp/cnu[1];
    dnu1 = cnunp;
    if( pinc < 18.0 )goto L220;
    cnu[2] = dnu2;
    et[2] = et2;
    et[1] = et1;
    cnu[1] = dnu1;
  L210:
    dene = et[3]/cnu[3];
    fune = etot-et[1]-et[2]-et[3]-et[5];
    goto L600;
  L220:
    cnu[2] = cpn-cpp;
    et[2] = et[1]*cnu[2]/cnu[1];
    et[1] = et[1]*cpp/cnu[1];
    cnu[1] = cpp;
    fectr = (18.0-pinc)/15.5;
    factr = 1.00001-fectr;
    cnu[2] = cnu[2]*fectr+dnu2*factr;
    cnu[1] = cnu[1]*fectr+dnu1*factr;
    et[2] = et[2]*fectr+et2*factr;
    et[1] = et[1]*fectr+et1*factr;
    goto L210;
    //
    // pi-plus + p collisions
    //
  L300:
    etot = einc+sumas;
    cnu[1] = abc;
    cnu[2] = witt2-abc;
    cnu[3] = pinc2*cpipnu[1][inx]+pinc*cpipnu[2][inx]+cpipnu[3][inx];
    cnu[5] = pinc2*cpimnu[1][inx]+pinc*cpimnu[2][inx]+cpimnu[3][inx];
    et[3] = pinc2*cpipk[1][inx]+pinc*cpipk[2][inx]+cpipk[3][inx];
    et[5] = pinc2*cpimk[1][inx]+pinc*cpimk[2][inx]+cpimk[3][inx];
    dene = et[3]/cnu[3];
    cnu[3] = cnu[3]+ab2;
    et[1] = dene*cnu[1];
    et[2] = dene*cnu[2];
    et[3] = dene*cnu[3];
    et[5] = et[5]*(1.0+cc3/cnu[5]);
    cnu[5] = cnu[5]+cc3;
    fune = etot-et[1]-et[2]-et[3]-et[5];
    goto L600;
    //
    // pi-zero + p collisions
    // you cant get here from there
    //
  L400:
    G4Exception("prob");
    //
    // pi-minus + p collisions
    //
  L500:
    etot = einc+sumas;
    cnu[1] = aaa;
    cnu[2] = wit-aaa;
    cnu[3] = pinc2*cpipnu[1][inx]+pinc*cpipnu[2][inx]+cpipnu[3][inx];
    cnu[5] = pinc2*cpimnu[1][inx]+pinc*cpimnu[2][inx]+cpimnu[3][inx];
    et[3] = pinc2*cpipk[1][inx]+pinc*cpipk[2][inx]+cpipk[3][inx];
    et[5] = pinc2*cpimk[1][inx]+pinc*cpimk[2][inx]+cpimk[3][inx];
    dene = et[3]/cnu[3];
    et[1] = dene*cnu[1];
    et[2] = dene*cnu[2];
    hold = et[5]/cnu[5];
    cnu[5] = cnu[5]+abk;
    et[5] = hold*cnu[5];
    et[3] = et[3]*(1.0+cc4/cnu[3]);
    cnu[3] = cnu[3]+cc4;
    fune = etot-et[1]-et[2]-et[3]-et[5];
  L600:
    cnu[4] = fune/dene;
    pb[1] = cnu[1];
    pb[2] = cnu[2]+pb[1];
    pb[3] = cnu[3]+pb[2];
    pb[4] = cnu[4]+pb[3];
    pb[5] = cnu[5]+pb[4];
    for ( i = 1; i < 6; i++ ) pb[i] /= pb[5];
    return;
  }
 
 void G4Cascade::prot(    G4float p0, 
			  G4float e, 
			  G4float p, 
			  G4float x, 
			  G4float y, 
			  G4float z )
  {
    //
    // selects from the ranft distribution over all momenta up to p;
    // over angles from 0 to 45-deg.0;
    // to be used for p+p and n+p inelastic events;
    //
    // p0 = momentum of incident particle (gev/c);
    // e  = energy of secondary nucleon (gev);
    // p  = momentum of secondary nucleon (gev/c);
    // x  = x-direction cosine;
    // y  = y-direction cosine;
    // z  = z-direction cosine;
    //
    G4float dum[18];
    G4double dndpm;   // is the maximum of dndp = a5 + a6/(pow(p0,a7))
    G4double dndp, b1, b2, b3, b4, b5, b6, b7, r1, r2, r3, r4, r5, csp, snp, r6;
    const G4double a1 = 0.885;           // a1, a2, and a3 are the ranft coefficients
    const G4double a2 = 0.101;
    const G4double a3 = 4.256;
    const G4double mp = 0.940075;        // mp is mass of proton
    const G4double mp2 = 0.883741005625; // mp2 is mp*mp
    const G4double a4 = 0.616849233;     // a4 is 45-deg. in rads squared
    const G4double a5 = 0.216470;
    const G4double a6 = 0.672377;
    const G4double a7 = 1.1004872;
    if( e != 0.0 )goto L4;
    r1 = p0;
    b1 = r1*r1;
    b2 = sqrt(1.0 + b1/mp2); // AH original dsqrt
    dndpm = a5 + a6/pow(r1,a7);
  L1:
    r2 = dflran();
    r2 = r1*r2;          // r2 = chosen p --- rejection follows
    b3 = r2*r2;
    b4 = sqrt(1.0 + b3/mp2);
    b5 = abs(r2 + r2*b2 - r1*b4);
    b6 = 1.0 + b2 - r1*r2/(mp2*b4);
    b7 = -a3*b3*a4;
    if( b7 < -50.0 )goto L2;
    b7 = exp(b7);
    b7 = 1.0 - b7;
    goto L3;
  L2:
    b7 = 1.0;
  L3:
    dndp = (pi/a3)*(a1/r1 +(a2/b1)*b5)*b6*b7;
    r3 = dndp/dndpm;
    r4 = dflran();
    if( r4 > r3 )goto L1;
    //
    // rejection completed---accept p = r2
    // next follows z
    //
  L6:
    r4 = dflran();
    r5 = (-1.0/(a3*b3))*log(1.0-r4*b7);
    z = cos(sqrt(r5));
    azio( csp, snp );
    r6 = sqrt(1.0 - z*z);
    x = r6*csp;
    y = r6*snp;
    p = r2;
    e = sqrt(b3 + mp2) - mp;
    return;
  L4:
    r2 = p;
    b3 = r2*r2;
    b7 = -a3*b3*a4;
    if( b7 < -50.0 )goto L5;
    b7 = exp(b7);
    b7 = 1.0 - b7;
    goto L6;
  L5:
    b7 = 1.0;
    goto L6;
  }
 
 void G4Cascade::rout1()
  {
if( verboseLevel > 1 ) G4cout << "   Entering rout1" << G4endl;
    out[11] = zee;
    G4float value2 = pow(zee,0.66666667);

    for ( G4int i = 5; i < 8; i++ )  // AH original G4int i = 5; 7; i++  -> add 1 to all limits
    {
      space[i+2] = out[i]*out[11];
      space[i+5] = out[i+3]*value2 + 7.0;
    }

 
    //
    // scaled protons per cc and potential proton well depth
    // (mev )in each region
    //
    out[12] = amasno-out[11];
    //
    // no. of neutrons n,  stored
    //
    value2 = pow(out[12],0.66666667);
    //
    // n 2/3
    //
    for ( i = 5; i < 8; i++ )
    {
      space[i-4] = out[i]*out[12];
      space[i-1] = out[i+3]*value2 + 7.0;
    }
    //
    // scaled neuts. per cc and pot. neut. well depth (mev)
    // in each region
    //
    for ( i = 1; i < 4; i++ )
    {
      hvn[i] = 0.5*space[i+3];
      hvp[i] = 0.5*space[i+9];
      awd[i] = hvn[i]+hvp[i];
      fvnp[i] = 0.5*awd[i];
      vnvp[i] = space[i+3]-space[i+9];
      pmac[i] = vnvp[i]-hvn[i];
      ppan[i] = -vnvp[i]-hvp[i];
      thpn[i] = hvp[i]-vnvp[i];
      ffptfn[i] = -vnvp[i]+fvnp[i];
      tffn[i] = space[i+9]-fvnp[i];
      tffp[i] = vnvp[i]+tffn[i];
    }

 

    G4float pppda = (2.0*zee)/(zee+amasno-1.0);
    G4float ppmda = (2.0*out[12])/(amasno+out[12]-1.0);
    G4float ppnda = (2.0*zee*out[12])/(amasno*amasno-amasno);
    G4float ppnna = (out[12]*out[12]-out[12])/(out[12]*out[12]+zee*zee-amasno);
    //
    // pion absorption calc. for each reg.  1/2 nwd,  (-pnan,  -ppac)
    // 1/2 pwd,  (-pnap,  -pnac )av. well depth,  1/4 av. well depth, 
    // (n-p)well depth,  (-vpvn),  1/2nwd -pwd,  (pmap,  -vphn),  1/2pwd -nwd, 
    // (-vnhp),  3/2pwd -nwd,  5/4pwd -3/4nwd,  3/4pwd -1/4nwd, 
    // 3/4nwd -1/4pwd,  prob. pip deut abs,  prob pim deut abs
    // prob pin deut abs,  prob pin nn abs rather than pp
    //
    G4int k = 15;
    for ( i = 4; i < 7; i++ )
    {
      out[k] = space[i+6] + einc;
      out[k+3] = space[i] + einc;
      --k;
    }
    //
    // total k.e. in mev incident proton(neutron)particle in each region
    //
    out[30] = zee/out[4]*1.4412 * 10E-13;
    //
    // coulomb potential at surface in mev.  conversion
    // factor = mev-cm per proton    (bg33)
    //
    // if( ctofe )90, 85, 90; // AH do functionality, ctofe ?
    G4int ctofe = 0; // AH dummy
  L85:
    ctofe = out[30];
    //
    // if ctofe = 0, then equate it to1/2 potential energy at surface
    //
  L90:
    for( i = 1; i < 4; i++ )
    {
      cfepn[i+3] = space[i+3] + ctofen;
      cfepn[i] = space[i+9] + ctofe;
    }
    //
    // bg33p--cutoff energies in each region for neutrons(protons)
    //
    G4int in = 0;
    G4float value1 = 6.28318531E+10 * pow(1.19366207E-1, 0.33333333);
    //
    // calc. of fermi momenta per cm.  pf equ. 2pi*((3/8pi)to 1/3)*e10
    //
    for ( i = 1; i < 4; i++ )
    {
      fmpn[i] = value1 * pow(space[i+6],0.33333333);
      fmpn[i+3] = value1 * pow(space[i],0.33333333);
    }
    //
    // fermi momenta per cm. of protons(neutrons)
    //
    for ( i = 1; i < 5; i++ ) rands[i] = randi[i];
    if( verboseLevel > 1 ) G4cout << "   Leaving rout1 " << G4endl;
    return;
  }
 

/* // AH rest of file commented for testing 
 void G4Cascade::rout2( G4double *t )
  {
    // dimension t(19)
    i1 = 1;
    value2 = einc+space[12];
    bovera( value2, pnms, ans );
    space[14] = 0.20d-24*ans;
    // pip+p;
    space[15] = 0.023d-24*ans;
    // (pim+p)el;
    space[16] = 0.0451d-24*ans;
    // (pim+p)ex;
    if( value1-100.0 )125, 125, 130;
  L125:
    fmax[1] = space[14];
    fmax[2] = space[15];
    fmax[3] = space[16];
    s[1] = 0.0;
    s[2] = 0.0;
  L12500:
    crdet( 1, t[1], einc );
    // (pip-p )absorption cross section;
    space[17] = crdt[1];
    if( i1 )200, 200, 12501;
  L12501:
    fmax[4] = space(17);
  L200:
    return;
  L130:
    if( value2-2600.0 )140, 140, 135;
    // value2 is greater than 2600;
  L135:
    i1 = 0;
    goto L200;
  L140:
    space[17] = 0.0;
    if( value2-220.0 )141, 141, 142;
  L141:
    s[1] = 0.75d-27*ans;
    // (pip-p)s.p. 400mev;
    s[2] = 4.7d-27*ans;
    // (pim-p)s.p. 400mev;
    goto L146;
  L142:
    if( value2-400.0 )145, 145, 160;
  L145:
    space[14] = 0.20d-24;
    space[16] = 0.0451d-24;
    s[1] = 7.8d-27*ans;
    s[2] = 21.8d-27*ans;
    // 660 mev;
  L146:
    if( einc-360.0 )12500, 200, 200;
  L160:
    if( value2-500.0)165, 165, 175;
  L165:
    space[14] = 0.113d-24;
    // 250 mev;
    space[15] = 20.5d-27*ans;
    // 620 mev;
    space[16] = 27.7d-27;
    // 250 mev;
    s[1] = 13.8d-27*ans;
    s[2] = 24.4d-27*ans;
    // 800 mev;
    i1 = -1;
    goto L146;
  L175:
    s[2] = 30.4d-27*ans;
    space[15]= 26.3d-27*ans;
    // 900;
    if( value2-600.0 )180, 180, 185;
    // 940, 325;
  L180:
    space[14] = 53.0d-27;
    // 325;
    // value2 lte 600;
    // 325;
    space[16] = 16.2d-27;
    s[1] = 15.2d-27*ans;
    // 940;
    goto L200;
  L185:
    if( value2-800.0 )190, 190, 195;
    // 1200, 400;
  L190:
    space[14] = 33.0d-27;
    // 400;
    space[16]= 12.0d-27*ans;
    // 400;
    s[1] = 20.9d-27*ans;
    // 1200;
    goto L200;
  L195:
    space[14] = 19.3d-27*ans;
    // 1300;
    space[16] = 8.2d-27*ans;
    // 540;
    s[1] = 23.3d-27*ans;
    // 1400;
    s[2] = 30.4d-27*ans;
    // sigma(a )900+2mb for future correction;
    goto L200;
    // values of fi for both incident and cascade particles for;
    // charged pions(pi +or- ), single production.  s[1] = n)s.p.,  s[2] = p);
    // s.p.,  space[14] = n)s,  space[15]= p)d.s.,  space[16]= p)abs,  space(17);
    //  = p)abs  (ppac);
  }
 
 void G4Cascade::rout3()
  {
    if( no-4 )205, 210, 210;
  L205:
    isw[11] = 1;
    goto L215;
  L210:
    isw[11] = 0;
  L215:
    undis();
    inc = 1;
    return;
  }
 
 void G4Cascade::rout4()
  {
    geo();
    if( i1 )10, 5, 5;
  L5:
    curr[3] = pnms;
    // pi+ or -mass/cm;
    curr[1] = no;
    partin();
    spac32( 32 );
  L10:
    return;
  }
 
 void G4Cascade::rout5( G4double *t, G4double *b, G4double *r )
  {
    //real*8 t(126), b(126), r(126)
    if( not-2 )240, 245, 250;
  L240:
    crjab( 1, t[1] );
    goto L254;
    // (pip-p)elastic scattering crs.0;
  L245:
    crjab( 1, b[1] );
    goto L254;
    // (pim-p)direct scattering crs.0;
  L250:
    crjab( 1, r[1] );
  L254:
    return;
  }
 
 void G4Cascade::rout6()
  {
    if( i3 )315, 345, 345;
  L315:
    abz = 1.0;
    med = clsm;
    knot = not;
    value1 = dflran();
    if( isw[11] )325, 320, 325;
  L320:
    if( value1-ppmda )330, 365, 365;
    // prob. pim-deut abs.0;
  L136:
    return;
  L365:
    i3 = 1;
    goto L136;
  L325:
    if( value1-pppda )330, 365, 365;
    // prob. pip-deut abs.0;
  L330:
    if( isw[11] )340, 335, 340;
  L335:
    it = 13;
    absec = pmac[med];
    goto L345;
  L340:
    it = 14;
    absec = -hvn[med];
  L345:
    strkp = -1.0;
    i1 = 0;
    i2 = med;
    bb();
    i3 = 0;
    goto L136;
  }
 
 void G4Cascade::rout6a()
  {
  L350:
    i1 = 0;
    spisom();
    strkp = -2.0;
    i1 = 1;
    spisom();
    strkp = -1.0;
    com = (awd[med]-7.0)*2.0*rcpmv;
    if( com-e[2] )350, 350, 355;
  L355:
    pm[2] = 2.0*dncms;
    pm[3] = dncms;
    e[2] = pm[2]+e[2];
    return;
  }
 
 void G4Cascade::rout7()
  {
    i3 = 0;
    if( curr[1]-3.0 )135, 380, 370;
  L370: if( curr[1]-5.0)385, 375, 135;
    // proton,  neutron not permitted;
  L135:
    i3 = -1;
  L136:
    return;
  L375:
    it = 7;
    ifca = 5;
    // pi meson - [5];
    absec = pmac(med);
    // pim +pp abs.  tyor = pmapp(20021);
    goto L400;
  L380:
    it = 10;
    ifca = 3;
    // tyor = ppan(20004 ) pip-nn abs. energy correction pimeson +;
    absec = ppan(med);
    goto L405;
  L385:
    value1 = dflran();
    if( value1-ppnna)390, 395, 395;
  L390:
    it = 8;
    ifca = 4;
    // pnann(20015 )= tyor  pin-nn abs  pimeson 0;
    absec = -hvn(med);
    goto L405;
  L395:
    it = 9;
    ifca = 2;
    // pnapp(20011 )= tyor  pin+pp abs.  pimeson 0;
    absec = -hvp(med);
  L400:
    strkp = -1.0;
    e[1] = wkrpn(med)*rcpmv+pm[1];
    goto L410;
  L405:
    strkp = -2.0;
    e[1] = wkrpn(med+3)*rcpmv+pm[1];
  L410:
    if( inc)420, 415, 420;
  L415:
    p1clc;
    goto L136;
  L420:
    p1cli;
    goto L136;
  }
 
 void G4Cascade::rout7a()
  {
    include 'COM.F';
  L425:
    i1 = -1;
    spisom;
    goto L(430, 435, 430, 430, 435), ifca;
  L430:
    value1 = space(med+3)-7.0;
    goto L440;
  L435:
    value1 = space(med+9)-7.0;
  L440:
    if( value1)135, 445, 445;
  L135:
    i3 = 2;
  L136:
    return;
  L445:
    if( (value1*2.0*rcpmv)-e[2])425, 425, 450;
  L450:
    pm[3] = dncms;
    pm[2] = 2.0*dncms;
    e[2] = pm[2]+e[2];
    value1 = ex;
    if( med-2)535, 510, 455;
  L535:
    i3 = 3;
    goto L136;
  L510:
    i3 = 4;
    goto L136;
  L455:
    if( inc)460, 470, 460;
  L460:
    if( isw[1])465, 269, 465;
  L269:
    i3 = 5;
    goto L136;
  L465:
    if( isw[2])289, 284, 289;
  L289:
    i3 = 6;
    goto L136;
  L284:
    i3 = 9;
    goto L136;
  L470:
    if( isw[1])475, 485, 475;
  L485:
    i3 = 7;
    goto L136;
  L475:
    if( isw[2])480, 490, 480;
  L480:
    i3 = 1;
    goto L136;
  L490:
    i3 = 8;
    goto L136;
  }
 
 void G4Cascade::rout8()
  {
    include 'COM3.F';
    i3 = 1;
    if( iv)585, 580, 580;
  L580:
    if( value1-value2)585, 585, 600;
  L585:
    if( isw[3])595, 590, 595;
  L590:
    ifc = 7+ifcc;
    // 7 = bg6e(2461 ) 8 = bg6ia(4026 ) ntnt(21626 ) bg48x(12762 )= 19;
    if( in)4640, 591, 4640;
  L4640:
    i3 = 2;
  L4641:
    return;
  L591:
    c[3] = d[2];
    goto L270;
  L270:
    i3 = 3;
    goto L4641;
  L595:
    ifc = 8+ifcc;
    if( in)530, 596, 530;
  L530:
    i3 = 4;
    goto L4641;
  L596:
    c[3] = d[2]+d[3]+d[4];
    goto L270;
  L600:
    signex;
    goto L4641;
  }
 
 void G4Cascade::rou10()
  {
    include 'COM3.F';
    i3 = 0;
    if( ex-d[4])635, 635, 630;
  L630:
    spac32(31);
    goto L605;
  L605:
    i3 = -1;
  L606:
    return;
  L635:
    curr[2] = out[15];
    wkrpn[1] = out[15];
    wkrpn[4] = out(18);
    // k.e. for protons and neutrons region 1;
    goto L606;
  }
 
 void G4Cascade::rou12()
  {
    include 'COM3.F';
    common/run/ke;
    i3 = 0;
    azio( sopc, sops);
    coll(0);
    if( col[15])135, 682, 135;
  L135:
    i3 = -1;
  L136:
    return;
  L682:
    if( ke)135, 683, 681;
  L681:
    com = (e[4]-dncms)/rcpmv;
    goto L684;
  L683:
    cole4;
  L684:
    i1 =  -1;
    value1 = com;
    if( pt[14]-2.0)710, 136, 136;
  L710:
    i3 = 1;
    goto L136;
  }
 
 void G4Cascade::rou13()
  {
    include 'COM3.F';
    i3 = 0;
    if( iv)695, 705, 705;
  L695:
    if( abz)700, 705, 700;
  L700:
    if( ifca-2)701, 703, 702;
  L701:
    in = 0;
  L350:
    i3 = 1;
    goto L426;
  L705:
    signex;
    if( ifc-12)706, 706, 707;
  L702:
    if( ifca-6)703, 701, 710;
  L703:
    in = 0;
  L425:
    i3 = -1;
  L426:
    return;
  L710:
    if( ifca-8)701, 704, 704;
  L704:
    in = 1;
    goto L350;
  L706:
    in = 0;
    goto L426;
  L707:
    if( ifc-18)708, 708, 709;
  L708:
    in = -1;
    goto L426;
  L709:
    in = 1;
    goto L426;
  }
 
 void G4Cascade::rou14()
  {
    include 'COM3.F';
    common/run/ke;
    common/numprt/numm, num, numm2, nhist, nodat, noevt, numa, numb;
    if( i3)861, 850, 720;
  L720:
    if( i1)725, 740, 740;
  L725:
    i1 = 0;
    value1 = (e[3]-pm[3])/rcpmv;
    if( abz)730, 735, 730;
  L730:
    value1 = value1+absec;
  L735:
    if( pt[2]-2.0)710, 685, 740;
    // pt[2] = 1 = proton   pt[2] = 2 = neutron   pt[2] = 3, 4, 5 = pion;
  L685:
    i3 = 2;
  L686:
    goto L1000;
  L710:
    i3 = 1;
    goto L686;
  L740:
    pinst;
    if( i1 != 0 )write(6, 136)i1, i3;
  L136:
    format(1h , 2x, 'rou14:just called pinst:i1 = ', i5, 2x, 'i3 = ', i5/);
    if( i1)135, 800, 135;
  L135:
    i3 = 3;
    goto L686;
  L800:
    i1 = 0;
    m = pt[2];
    value2 = value1;
    if( m-3)805, 820, 820;
  L805:
    if( eco(m)-value2)815, 806, 810;
  L806:
    if( eco(m))135, 135, 810;
  L810:
    pt(i1+3 )= 0.0;
    pnbc(m )= pnbc(m)+1.0;
    goto L834;
  L815:
    pt(i1+3 )= value2;
    if( i1)860, 835, 860;
  L820:
    ccofe = clcfe;
    if( m-4)821, 822, 822;
  L821:
    if( strkp+2.0)826, 824, 826;
  L822:
    if( strkp+2.0)826, 826, 823;
  L823:
    ccofe = clcfe-ctofe+ctofen;
    goto L826;
  L824:
    ccofe = clcfe+ctofe-ctofen;
  L826:
    if( value2-ccofe)810, 810, 825;
  L825:
    if( strkp+2.0)830, 830, 815;
  L830:
    pt[3] = value1-space(med+3)+space(med+9);
  L834:
    if( i1)845, 835, 845;
  L835:
    m = pt[14];
    if( m-3)840, 135, 135;
  L840:
    value2 = com;
      i1 = 12;
      goto L805;
  L845:
      if( pt[3])860, 850, 860;
  L850:
      punp;
      //if( i1 < 0)write(6, 851)i1;
      //851 format(1h , 2x, 'rou14:just called punp:i1 = ', i5/);
      if( i1)135, 1340, 4415;
// -,   = error  0 = end of record  + = piscc(6607);
 L1340:
      i3 = 4;
      goto L686;
 L4415:
      i3 = 5;
      goto L686;
  L860:
      collm(-1);
  L861:
      if( ke > 0)goto L1000;
       stpr;
      if( i1)135, 850, 135;
  L1000:
      return;
  }
 
 void G4Cascade::rou17( t, b, r, w, g )
  {
      dimension t(117), b(101), r(117), w(234), g(234);
      include 'COM3.F';
      real *8 t, b, r, w, g;
      common/numprt/numm, num, numm2, nhist, nodat, noevt, numa, numb;
      if( i3)2001, 2001, 2040;
 L2001:
      pt(38 )= 0.0;
// 1-alpha  part.6+1;
      value1 = rlke-180.0;
      crdet(1, t[1], value1);
      com2 = crdt[1];
      ftr = dncms*rlke*2.0*rcpmv+2.9877156d27;
// e**2 = min**2+ncnms*rlke*2*rcpmv;
      univer = sqrt(ftr);  // AH orignal dsqrt
// e;
 L2005:
      value2 = dflran();
      com = value2*com2;
// r-prime;
      gene(b[1]);
      com1 = (com*com+ftr-.501264d26)/(2.0*univer);
// m1r prime)**2+e**2-2(pnms)/2e = e alpha;
      a = com1*com1-com*com;
      if( a)2006, 2009, 2009;
 L2006:
      pacnt = pacnt+1.0;
      goto L2005;
 L2009:
      unive = ((univer-com1)*com1/univer)*sqrt(a);
// ((e beta*e alpha*p alpha)/e )= f(m, tr);
      crdet(1, r[1], value1);
// (pi-nuc)fmax(rlke)isobar sampling s.p.0;
      com1 = dflran();
      if( (unive/crdt[1])-com1)2005, 2010, 2010;
// random no. less or equal than f(m, tr)/fmax(tr);
 L2010:
      angid;
      pm[3] = com;
      pm[4] = poms;
      pt[2] = 3.0;
      pt[4] = poms;
      pt[14] = 3.0;
      pt[16]= poms;
      pt(26 )= 1.0;
      pt(28 )= dncms;
      if( isw[9])2020, 2015, 2020;
 L2015:
      if( isw[10])2030, 2025, 2030;
 L2020:
      if( isw[10])2035, 135, 2035;
  L135:
      i3 = -1;
  L136:
      return;
 L2025:
      value1 = 0.4;
      value2 = 6.6666667d-1;
      value3 = 0.0;
      goto L2037;
 L2030:
      crdet(2, w[1], value1);
// (pich-p)fract. fin.sta.with recl. pi1 pi0 l.e.0;
      value3 = 3.3333333d-1;
      goto L2036;
 L2035:
      crdet(2, g[1], value1);
      value3 = strkp;
// (pin-p)fract.fin.sta.with recl.pi1 pio l.e.0;
 L2036:
      value1 = crdt[1];
      value2 = crdt[2];
 L2037:
      alpha;
 L2040:
      ecpl;
      if( i1)2045, 2045, 135;
 L2045:
      coll(-1);
      if( col[15])135, 2050, 135;
 L2050:
      if( pt(38))2084, 2055, 2084;
 L2084:
      i3 = 0;
      goto L136;
 L2055:
      pt(39 )= 0.0;
      pt[3] = ((e[4]-pm[4])/rcpmv)+pt[3];
      i3 = 1;
      goto L136;
  }
 
 void G4Cascade::rou18()
  {
      include 'COM3.F';
      goto L(4010, 4015, 4035, 4095, 4071, 4080), i3;
 L4010:
      i = 3;
      col[15]= 1.0;
      k = 27;
      goto L4020;
 L4015:
      i = 3;
      col[15]= 4.0;
      k = 15;
 L4020:
      pnidk[1] = pm[i];
      j = i;
      for ( G4int l = 2; i < 5; i++ )
      {
        pnidk(l )= pxyz(j);
        j += 4;
      }
      pnidk[5] = e[i];
      pnidk(6 )= pt(k-11);
      idk;
      if( k-27)4031, 4030, 4031;
 L4030:
      pt[15]= pt[15]+((pnidk[12]-pnidk( 6))/rcpmv);
 L4031:
      pt(k )= pt(k)+((pnidk(13)-dncms)/rcpmv);
      i3 = 1;
 L2057:
      iv = k;
      return;
 L4035:
      k = 3;
      col[15]= 2.0;
      if( pt[2]-3.0)4071, 4039, 4039;
 L4039:
      if( pt(k)-2500.0)4040, 4040, 4055;
 L4055:
      i3 = 5;
      goto L2057;
 L4040:
      if( pt(k))4050, 4050, 4045;
 L4045:
      ccofe = eco[1];
      if( pt(k-1)-4.0 )4047, 4046, 4046;
 L4046:
      ccofe = ccofe - ctofe + ctofen;
 L4047:
      if( pt(k )- ccofe  )4050, 4050, 4070;
 L4050:
      m = pt(k-1);
      pnbc(m )= pnbc(m)+1.0;
      pt(k )= 0.0;
      i3 = 3;
      goto L2057;
 L4070:
      if( k-3)4071, 4071, 4095;
 L4071:
      col[15]= 3.0;
      k = 15;
      if( pt[14]-2.0)135, 135, 4039;
  L135:
      i3 = 2;
      goto L2057;
 L4080:
      l = 14;
      for ( G4int m = 5; i < 8; i++ )
      {
        pt(m )= pnidk(l);
        pt(m+12 )= pnidk(l+3);
        ++l;
      }
      pt[11] = pnidk[12];
      pt(12 )= pnidk[6];
      i = 4;
      k = 39;
      col[15]= 5.0;
      goto L4020;
 L4095:
      i1 = 3;
 L4100:
      k = 12*i1-33;
      if( i1-4)4110, 4115, 4120;
 L4110:
      i2 = -1;
      goto L4130;
 L4115:
      i2 = 0;
      goto L4130;
 L4120:
      if( i1-5)4115, 4125, 4235;
 L4235:
      i3 = 4;
      goto L2057;
 L4125:
      i2 = 1;
 L4130:
      if( pt(k))4135, 4140, 4135;
 L4135:
      pstor;
 L4140:
      i1 = i1+1;
      goto L4100;
  }
 
 void G4Cascade::rou19()
  {
      include 'COM3.F';
      pt[3] = pt[3]+((pt[11]-pt[12])/rcpmv);
// collision allowed;
      k = 3;
 L4135:
      if( pt(k)-2500.0)4150, 4150, 4140;
 L4140:
      i3 = 1;
      goto L136;
 L4150:
      if( pt(k))4160, 4160, 4155;
 L4155:
      ccofe = eco[1];
      if( pt(k-1)-4.0 )4157, 4156, 4156;
 L4156:
      ccofe = ccofe - ctofe + ctofen;
 L4157:
      if( pt(k )- ccofe )4160, 4160, 4200;
 L4160:
      pt(k )= 0.0;
      if( pt(k-1)-3.0)135, 4170, 4165;
  L135:
      i3 = -1;
  L136:
      return;
 L4165:
      if( pt(k-1)-5.0)4170, 4170, 135;
 L4170:
      m = pt(k-1);
      pnbc(m )= pnbc(m)+1.0;
      goto L4185;
 L4175:
      i2 = 2;
 L4176:
      i1 = (k/12)+3;
      pstor;
 L4185:
      if( k-15)4190, 4195, 4210;
 L4190:
      k = 15;
      if( pt[15])4195, 4195, 4175;
 L4195:
      k = 27;
      pt(27 )= pt(27)+((pnidk[12]-pt(k+1))/rcpmv);
      goto L4135;
 L4200:
      if( k-15)4175, 4205, 4205;
 L4205:
      i2 = 0;
      goto L4176;
 L4210:
      if( k-27)135, 4215, 4235;
 L4215:
      if( pt(39))4235, 4235, 4220;
 L4235:
      i3 = 0;
      goto L136;
 L4220:
      i2 = 1;
      k = 39;
      goto L4176;
  }
 
 void G4Cascade::rou21( v, w, x, y, z )
  {
    dimensionv(161), w(101), x(161), y(130), z(176);
      include 'COM3.F';
      real *8 v, w, x, y, z;
      value2 = rlke*4.81633308d24+9.0554256d27;
// e(tr)**2 = rlke*rcpmv*2*ncms+4*ncms**2;
      value3 = sqrt(value2);
      goto L(4330, 4360, 4365, 4410), i3;
 L4330:
      isw(12 )= 0;
 L4331:
      pt(38 )= 0.0;
      i1 = 0;
      ans = rlke;
 L4333:
      value1 = ans-300.0;
      crdet(1, v[1], value1);
// (nuc-nuc )f(tr )isobar sampling;
      ftr = crdt[1];
 L4335:
      sn = dflran();
      com = sn*ftr;
// r prime = f(tr)*random;
      gene(w[1]);
// (nuc-nuc)mass of isobar s.p.    m(r prime);
      if( i1)4370, 4336, 4375;
 L4336:
      com1 = (com*com-sqnm+value2)/(2.0*value3);
// e gamma;
      a = com1*com1-com*com;
      if( a)4337, 4338, 4338;
 L4337:
      pgcnt = pgcnt+1.0;
      goto L4335;
 L4338:    univer = sqrt(a)*com1*(value3-com1)/value3;
// f(m, tr )= p gamma*e gamma*e delta/e;
      crdet(1, x[1], value1);
// (nuc-nuc)fmax(tr )isobar sampling s.p.0;
      com1 = dflran();
      if( com1-(univer/crdt[1]))4340, 4340, 4335;
 L4340:
      pm[4] = dncms;
      pm[3] = com;
      angid;
      pt[4] = dncms;
      pt(28 )= dncms;
      alp19;
 L2040:
      return;
 L4360:
      isw(12 )= 2;
      goto L4331;
 L4365:
      isw(13 )= 0;
 L4366:
      i1 = -1;
      ans = ((value3-pnms)**2-9.0554256d27)/4.81633308d24;
      goto L4333;
// tr prime     com1 = rlke prime;
 L4370:
      com1 = ((value3+dncms-com)**2-9.0554256d27)/4.81633308d24;
      com2 = com;
      ans = com1;
      com4 = ftr;
      i1 = 1;
      goto L4333;
 L4375:
      com1 = (com2*com2-com*com+value2)/(2.0*value3);
// e epsilon;
      a = com1*com1-com2*com2;
      if( a)4376, 4377, 4377;
 L4376:
      pecnt = pecnt+1.0;
      goto L4380;
// f(m1, m2, tr )= p epsilon*e epsilon*e zeta/e;
 L4377:
      univer = sqrt(a)*com1*(value3-com1)/value3;
      value1 = rlke-920.0;
      crdet(1, y[1], value1);
// (nuc-nuc)fmax(tr )isobar sampling d.p.  fmax(m1, m2, tr);
      value1 = dflran();
      if( sn-(univer*ftr/(crdt[1]*com4)))4385, 4385, 4380;
 L4380:
      ftr = com4;
      i1 = -1;
      goto L4335;
 L4385:
      value1 = dflran();
      if( value1-.5)4390, 4390, 4395;
 L4390:
      pm[3] = com2;
      pm[4] = com;
      goto L4400;
 L4395:
      pm[3] = com;
      pm[4] = com2;
 L4400:
      angid;
      pt[16]= dncms;
      pt(40 )= dncms;
      if( isw(13))4401, 4405, 4401;
 L4401:
     crdet(1, z[1], rlke);
      value1 = crdt[1];
// (n-p)fract.fin.sta.3/2 l.e.0;
 L4405:
      pt[2] = 3.0;
      pt[4] = poms;
      pt[14] = 1.0;
      pt(26 )= 3.0;
      pt(28 )= poms;
      pt(38 )= 1.0;
      alp28;
      goto L2040;
 L4410:
      isw(13 )= 2;
      goto L4366;
  }
 
 void G4Cascade::rou22( G4double *v, G4double *w, G4double *x, G4double *y, G4double *z )
  {
    // dimension v(19), w(19), x(126), y(126), z(126)
    io = 0;
    i5 = curr[1];
    goto (4480, 4830, 4500, 4515, 4535, 4825, 4640, 4840, 4855, 4620, 4885, 4611, 
          4544, 4638, 4955, 4415, 4696, 4670, 4479, 4650, 4610, 4870, 5005), iv;
  L4415:
    for ( G4int i = 1; i < 4; i++ )
    {
      xi[i] = curr[i+3];
      dcos[i] = curr[i+6];
    }
    in = -1;
    med = curr[10];
    i5 = curr[1];
    geo();
    if( i1)5046, 4424, 4424;
 L4424:
    if( curr[1]-2.0)4430, 4690, 4425;
 L4425:
    if( curr[1]-4.0)4695, 4696, 5045;
 L5045:
    iv = 1;
 L5046:
    return;
 L4430:
    isw[4] = 1;
// proton;
 L4431:
    abz = 0.0;
      i4 = 1;
      i2 = med;
      if( i2-1)135, 4436, 4470;
 L4436:
      i4 = 6;
      goto L4470;
 L4438:
      isw(6 )= 1;
      isw[5] = 1;
      i3 = 0;
      if( com-3600.0)4440, 4440, 135;
  L135:
      iv = 2;
      goto L5046;
// com = greatest energy inside nucleus for nucleons only;
 L4440:
      if( com-560.0)4445, 4445, 4450;
 L4445:
      if( com-160.0)4460, 4460, 4455;
 L4450:
      dfmax;
// dfmax fills out fmax(1-6);
      i1 = 6;
 L4451:
      if( curr[1]-2.0)4452, 4457, 4457;
 L4457:
      i3 = 1;
 L4452:
      store;
      ex = 0.0;
      signex;
      goto L(4479, 4535, 4610, 4630, 4544, 4665, 4810, 4825, 4855, 605, 4955), i4;
  L605:
      iv = 3;
      goto L5046;
 L4455:
      isw(6 )= 0;
      pfmax;
// pfmax fills out fmax(1-4);
      i1 = 4;
      goto L4451;
 L4460:
      isw[5] = 0;
      isw(6 )= 0;
      nn;
// nn fills out fmax(1-2), fmax(3-6 )= 0;
      i1 = 2;
      goto L4451;
 L4470:
      isw[1] = 0;
      isw[2] = 0;
      isw[3] = 0;
 L4471:
      m = med+ int(15.0-6.0*sngl(curr[1]));
      a = curr[2]-space(m);
      4472 for ( G4int i = 1; i < 4; i++ )
      {
        wkrpn[i] = a+space[i+9];
        wkrpn(i+3 )= a+space[i+3];
      }
 L4476:
      m = 4-3*isw[4];
      com = wkrpn(m);
      goto L4438;
 L4479:
      goto L(4890, 4890, 4480), i2;
 L4480:
      if( ex-d[2])4495, 4495, 4485;
// entry point for cascade charged and neutral pions after crjab;
// rejection, prod.reaction lt 180 fermi rejection for non-absorption;
// reaction also all cascade nuclear rejections including crjab or;
// rlke lt prod.threshold;
 L4485:
      if( d[3])4660, 4490, 4660;
 L4490:
      ccpes;
      if( i1)850, 850, 135;
// goto L Lpunp;
  L850:
      iv = 4;
      goto L5046;
 L4495:
      if( in )4951,  4951, 4965;
 L4951:
      bg6ca(3, 0);
      ifcc = 12;
 L4496:
      med = clsm;
      iv = 5;
      goto L5046;
 L4500:
      value1 = ex;
      iv = 6;
      goto L5046;
 L4515:
      value1 = ex+d[3];
  L500:
      iv = 7;
      goto L5046;
// e.p.for rejection following crjab or rlke lt 180 in prod.reactions;
// and for fermi rejection in scattering and production reactions;
// for cascade charged and neutral pions;
 L4535:
      if( ex-d[6])4495, 4495, 4490;
 L4540:
      isw[1] = 1;
 L4543:
      goto L(4544, 4544, 4544, 4825), i5;
// e.p.for all rejections, crjab, production, and fermi for cascade;
// nucleons;
 L4544:
      if( ex-d[3])4635, 4635, 4545;
 L4545:
      i2 = 3;
      i4 = 2;
      i3 = 1;
      i1 = 2;
      if( d[4])4550, 4476, 4550;
 L4550:
      isw[2] = 1;
      isw[3] = 1;
      i4 = 3;
      i2 = 1;
      goto L4476;
 L4610:
      goto L(4611, 4611, 4611, 4855), i5;
 L4611:
      if( ex-d[4])4615, 4615, 4625;
 L4615:
      bg6ca(1, 0);
      ifcc = 7;
      goto L4496;
 L4620:
      value1 = ex;
      iv = 8;
      goto L5046;
 L4625:
      i4 = 4;
      i1 = 2;
      i2 = 2;
      i3 = 1;
      goto L4476;
// e.p.for cascade nucleons and pi0, after crjab and prod.threshold;
// rejections;
 L4630:
      if( i5-2)4638, 4638, 4634;
 L4634:
      if( i5-4)135, 4955, 135;
// e.p.for cascade nucleons after all fermi rejections;
 L4638:
      if( ex-d[5])4635, 4635, 4655;
 L4635:
      bg6ca(2, 0);
      ifcc = 10;
      goto L4496;
 L4640:
      value1 = ex;
      goto L500;
 L4650:
      if( isw[3])4630, 4543, 4630;
 L4655:
      i4 = 2;
      i2 = 3;
 L4656:
      i3 = 1;
      i1 = 2;
      goto L4476;
 L4660:
      isw[1] = 1;
      if( curr[1]-3.0)4661, 4815, 4930;
 L4661:
      i4 = 5;
      i2 = 2;
      goto L4656;
// d.p.0;
 L4665:
      isw[1] = 1;
      isw[2] = 1;
      isw[3] = 1;
      goto L4479;
 L4670:
      any = fmax(not);
      if( not-5)4684, 1290, 4680;
 L1290:
      iv = 9;
      goto L5046;
 L4680:
      if( i5-4)1290, 4685, 1290;
 L4684:
      if( knot-15)4685, 5000, 5000;
 L4685:
      goto L(4686, 4686, 4686, 4985), i5;
 L4686:
      if( not-2)1161, 1157, 1270;
 L1161:
      iv = 10;
      goto L5046;
 L1157:
      iv = 11;
      goto L5046;
 L1270:
      iv = 12;
      goto L5046;
 L4690:
      isw[4] = 0;
// neutron;
      goto L4431;
 L4695:
      isw[11] = 1;
// curr[1] = 3 = pi+ meson(11575)--pi 0 = 4(15100)--pi -(14646);
 L4696:
      in = 1;
      isw[1] = 0;
      isw[2] = 0;
      isw[3] = 0;
      isw[5] = 1;
      isw(6 )= 0;
      isw(7 )= 0;
      isw(8 )= 1;
      i6 = i5-2;
      i2 = med;
      com = curr[2]-space(med+9);
      goto L(4706, 4705, 4707), med;
// e.p.for pi0 after fermi rejection;
 L4705:
      isw(8 )= 0;
 L4706:
      isw(7 )= 1;
 L4707:
      for ( G4int i = 1; i < 4; i++ )
 {
   wkrpn[i] = com+space[i+9];
   wkrpn(i+3 )= com+space[i+3];
 }
      com = com+space[4];
 L4710:
      if( com-2600.0)4715, 4715, 135;
 L4715:
      if( com-100.0)4716, 4716, 4717;
 L4716:
      lg = 4;
      goto L4910;
 L4717:
      lg = 6;
 L4720:
      spcn;
      if( value1)4725, 4725, 4734;
 L4725:
      com = curr[2]-space(med+9);
      if( com-360.0)4730, 4734, 4734;
 L4730:
      goto L(4731, 4915, 4731), i6;
 L4731:
      crdet(1, v[1], com);
      fmax[4] = crdt[1];
 L4734:
      goto L(4735, 4920, 4735), i6;
 L4735:
      i1 = 6;
 L4736:
      i4 = 7;
      i2 = 1;
      i3 = 0;
      goto L(4755, 4760, 4755), i6;
 L4755:
      if( isw[11])4760, 4790, 4760;
// pi+;
 L4760:
      if( isw[7])4765, 4785, 4765;
// pi-;
 L4765:
      if( isw[8])4452, 4780, 4452;
 L4780:
      i2 = 2;
      goto L4452;
 L4785:
      i2 = 3;
      goto L4452;
 L4790:
      i3 = 1;
      goto L4760;
 L4810:
      if( curr[1]-2.0)4811, 4811, 4479;
 L4811:
      a = curr[2]-space(med+9);
      goto L4472;
 L4815:
      i4 = 8;
 L4816:
      i2 = 2;
 L4818:
      i3 = 0;
      i1 = lg;
      if( curr[1]-4.0)4819, 4827, 4819;
 L4819:
      if( isw[11])4821, 4820, 4821;
 L4820:
      i3 = 1;
 L4821:
      if( curr[1]-3.0)4822, 4823, 4823;
 L4822:
      iv = 23;
      goto L5046;
 L4823:
      if( lg-4)4822, 4452, 4824;
 L4824:
      m = 5-iabs(i5-4);
      univer = fmax(m);
      spcn;
      fmax(m )= univer;
      goto L4452;
 L4827:
      i1 = i1+1;
      goto L4823;
// e.p.for charged and neutral pions after crjab rejection also, after;
// rlke <= 180 in prod.reactions and after fermi rejection in scatter;
// ing and prod.reactions;
 L4825:
      if( ex-d[3])4826, 4826, 4845;
 L4826:
      goto L(576, 4940, 576), i6;
  L576:
      iv = 13;
      goto L5046;
 L4830:
      any = fmax(not);
      if( curr[1]-3.0)4831, 4832, 4832;
 L4831:
      ifc = 12;
L23503:
      iv = 14;
      goto L5046;
 L4832:
      ifcc = (clsm-2.0)*((clsm*5.5)-8.5)+12.05;
      goto L23503;
 L4840:
      if( i4-10)4841, 4865, 4841;
// e.p.for escape prior to choosing reactions--cascade charged pion;
 L4841:
      if( clsm-2.0)4842, 4865, 4842;
 L4842:
      i4 = 2;
      goto L4816;
 L4845:
      i4 = 2;
// e.p.when cascade particle escapes from region 2;
      i2 = 3;
      if( d[4])4850, 4846, 4850;
 L4846:
      goto L(4818, 4960, 4818), i6;
 L4850:
      isw[2] = 1;
      isw[3] = 1;
      i4 = 9;
      i2 = 1;
      goto L4846;
// e.p.after all rejections except fermi in abs.reactions for;
// charged and neutral cc pions--escape from region 1 prior to;
// choosing reaction for cc charged pion;
 L4855:
      if( ex-d[4])4856, 4856, 4860;
 L4856:
      goto L(636, 4970, 636), i6;
  L636:
      iv = 15;
      goto L5046;
 L4860:
      i4 = 10;
      goto L(4816, 4950, 4816), i6;
 L4865:
      i4 = 2;
      i2 = 3;
      goto L4818;
 L4870:
      if( in)4875, 356, 4875;
  L356:
      iv = 16;
      goto L5046;
 L4875:
      ifca = 8*iabs(i6-2)-11*(i6-1)*(i6-3);
      if( isw[1])4880, 4500, 4880;
 L4880:
      if( isw[2])480, 4515, 480;
  L480:
      iv = 17;
      goto L5046;
 L4885:
      ifca = 10*iabs(i6-2)+12*(i6-1)*(3-i6);
      if( isw[3])530, 4640, 530;
  L530:
      iv = 18;
      goto L5046;
 L4890:
      if( curr[1]-3.0)4895, 4900, 4900;
 L4895:
      goto L(4610, 4540), med;
 L4900:
      isw[1] = 1;
      goto L(4905, 4825), med;
 L4905:
      isw[2] = 1;
      isw[3] = 1;
      goto L4855;
 L4910:
      isw[5] = 0;
      goto L(4911, 4980, 4911), i6;
 L4911:
      bovera(curr[2], pnms, ans);
      fmax[1] = 0.20 d-24*ans;
// (pi+p)scattering;
      fmax[2] = 23.0d-27*ans;
// (pim+p)scattering;
      fmax[3] = 45.1d-27*ans;
// (pim+p)exchange;
      com = curr[2]-space(med+9);
// (k.e.of pions outside nucleus;
      crdet(1, v[1], com);
      fmax[4] = crdt[1];
// c(pip+p)abs.0;
      i1 = 4;
      goto L4736;
 L4915:
      crdet(1, w[1], com);
      fmax[5] = crdt[1];
 L4920:
      i1 = 7;
      goto L4736;
 L4925:
      bg6ca(3, 4);
      ifcc = 24;
L11410:
      ka = 7;
      med = clsm;
      iv = 19;
      goto L5046;
 L4930:
      i4 = 8;
 L4931:
      i2 = 2;
      goto L4818;
 L4940:
      bg6ca(2, 3);
 L4945:
      ifcc = 21;
      goto L11410;
 L4950:
      i4 = 11;
      goto L4931;
 L4955:
      if( ex-d[5])4940, 4940, 4975;
 L4960:
      i1 = 9-i5;
      goto L4818;
 L4965:
      goto L(231, 4925, 231), i6;
  L231:
      iv = 20;
      goto L5046;
 L4970:
      bg6ca(1, 2);
      goto L4945;
 L4975:
      i1 = 5;
      goto L4865;
 L4980:
      bovera(curr[2], poms, ans);
      fmax[1] = 89.2d-27*ans;
// (pi0+p)elast.scat.0;
      fmax[2] = 45.1d-27*ans;
// (pi0+p)ex.scat.0;
      fmax[3] = fmax[1];
      space(48 )= fmax[1];
      fmax[4] = fmax[2];
      space(49 )= fmax[2];
      com = curr[2]-space(med+9);
      crdet(1, w[1], com);
// (pin-p)abs. crs.sec.0;
      fmax[5] = crdt[1];
// (pi0+p)abs.0;
      space(50 )= fmax[5];
      i1 = 5;
      goto L4736;
 L4985:
      if( not-2)4990, 4995, 5035;
 L5035:
      iv = 21;
      goto L5046;
 L4990:
      crjab(1, x[1]);
// (pin-p)direct scat.crs.0;
 L1170:
      iv = 22;
      goto L5046;
 L4995:
      crjab(1, y[1]);
// (pim-p)xch.scat.crs.0;
      goto L1170;
 L5000:
      if( knot-16)5005, 4995, 4995;
 L5005:
      crjab(1, z[1]);
// (pin-n)drct.scat.crs.0;
      goto L1170;
  }

  */
