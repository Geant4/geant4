#ifndef G4_HETC
#define G4_HETC

 class G4Cascade
 {
 public:

   G4Cascade();
   ~G4Cascade();

   inline  void setVerboseLevel(G4int level){verboseLevel = level;}
   G4int go(void);

 private:

   G4int verboseLevel;
   G4float efas[80];
   G4float alpfas[80];
   G4float betfas[80];
   G4float gamfas[80]; 
   G4float wtfas[80];
   G4float pb[5];
   G4float st0r;
   G4int it[80];
   G4int numnuc;

   void nucnuc(
	       G4float p0, 
	       G4int nofas, 
	       G4int itype
	       );

   void pinuc(  
	      G4float p0, 
	      G4int nofas, 
	      G4int itype 
	      );

   void prob(
	     G4int itype, 
	     G4float p0 
	     );

   void prot( 
	     G4float p0, 
	     G4float e, 
	     G4float p, 
	     G4float x, 
	     G4float y, 
	     G4float z );


   G4float out[18]; // AH Parameters for rout1
   G4float zee;
   G4float space[12];
   G4float amasno;
   G4float hvn[3];
   G4float hpw[3];
   G4float hvp[3];
   G4float awd[3];
   G4float fvnp[3];
   G4float vnvp[3];
   G4float pmac[3];
   G4float ppan[3];
   G4float thpn[3];
   G4float ffptfn[3];
   G4float tffn[3];
   G4float tffp[3];
   G4float einc;
   G4float cfepn[6];
   G4float ctofen;
   G4float fmpn[6];
   G4float rands[4];
   G4float randi[4];

   void rout1();
   /* // AH commented for testing

      void rout2();
      void rout3();
      void rout4();
      void rout5();
      void rout6();
      void rout6a();
      void rout7();
      void rout7a();
      void rout8();
      void rou10();
      void rou12();
      void rou13();
      void rou14();
      void rou17();
      void rou18();
      void rou19();
      void rou21();
      void rou22();

      */

   void pol1( G4double &cs, G4double &si );

   void azio( G4double &s, 
	      G4double &c );

   void piAH(G4float po, 
	     G4int itip, 
	     G4float ener, 
	     G4float pmom, 
	     G4float a,   
	     G4float b, 
	     G4float g );

   inline G4double dflran()
     {
         G4double value = G4UniformRand();
        while( value == 0.0 || value == 1.0 ) value = G4UniformRand();
       return value;
     }
 };
#endif
