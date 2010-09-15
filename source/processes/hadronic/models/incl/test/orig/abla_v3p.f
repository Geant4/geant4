c********************************************************************  
C Version modifiee par A.B. pour entrer dans le calcul dynamique
C     Cugnon(2 ou 3) + KHS (ici V3).
C Modifications pour le tirage (maxwell) de l'energie dans les 
C     petites cibles. (Sub Direct) Juin/2002 C.V. A.B. 
C *******************************************************************                                                                      
C                                                                       
C      SUBROUTINE ABLAINIT(STATUS,TSTAT,NAME,FPATH)
C********************************************************************                      

      SUBROUTINE INIT_EVAPORA(RACINE)
                      
C********************************************************************                                                                       
C     ON INPUT:  INPUT PARAMETERS FROM FILE                             
C---------------------------------------------------------------------  
C     ON OUTPUT:                                                        
C     STATUS - FLAG FOR END OF INPUT FILE                               
C     TSTAT  - FLAG FOR NTUPLE-OUTPUT                                   
C     NAME   - NAME FOR ISOTOPIC PRODUCTION CROSS SECTION FILES         
C     FPATH  - PATH FOR  "          "        "        "    "            
C---------------------------------------------------------------------
C
C     Modification 5-january-2000 by KHS and BJ
C
C     New treatment of dissipation. 
C     See report of Beatriz Jurado, Jan. 2000
C
C---------------------------------------------------------------------
C   
C     MODIFICATION 6-aug-1999 by JB and MVR
C
C Some problems arised from an uncorrect evaluation of the fission barrier
C { 1) shell correction ( ECGNZ(J,K) ) was not subctracted (20-jul-99)
C   2) fiss. barrier EF was calc. before ang. mom. correct. (6-aug-99) }
C
C---------------------------------------------------------------------  
C     PROJECTILE AND TARGET PARAMETERS + CROSS SECTIONS                 
C     COMMON /ABLAMAIN/ AP,ZP,AT,ZT,EAP,BETA,BMAXNUC,CRTOT,CRNUC,       
C                       R_0,R_P,R_T, IMAX,IRNDM,PI,                     
C                       BFPRO,SNPRO,SPPRO,SHELL                         
C                                                                       
C     AP,ZP,AT,ZT   - PROJECTILE AND TARGET MASSES                      
C     EAP,BETA      - BEAM ENERGY PER NUCLEON, V/C                      
C     BMAXNUC       - MAX. IMPACT PARAMETER FOR NUCL. REAC.             
C     CRTOT,CRNUC   - TOTAL AND NUCLEAR REACTION CROSS SECTION          
C     R_0,R_P,R_T,  - RADIUS PARAMETER, PROJECTILE+ TARGET RADII        
C     IMAX,IRNDM,PI - MAXIMUM NUMBER OF EVENTS, DUMMY, 3.141...         
C     BFPRO         - FISSION BARRIER OF THE PROJECTILE                 
C     SNPRO         - NEUTRON SEPARATION ENERGY OF THE PROJECTILE       
C     SPPRO         - PROTON    "           "   "    "   "              
C     SHELL         - GROUND STATE SHELL CORRECTION                     
C---------------------------------------------------------------------  
C                                                                       
C     ENERGIES WIDTHS AND CROSS SECTIONS FOR EM EXCITATION              
C     COMMON /EMDPAR/ EGDR,EGQR,FWHMGDR,FWHMGQR,CREMDE1,CREMDE2,        
C                     AE1,BE1,CE1,AE2,BE2,CE2,SR1,SR2,XR                
C                                                                       
C     EGDR,EGQR       - MEAN ENERGY OF GDR AND GQR                      
C     FWHMGDR,FWHMGQR - FWHM OF GDR, GQR                                
C     CREMDE1,CREMDE2 - EM CROSS SECTION FOR E1 AND E2                  
C     AE1,BE1,CE1     - ARRAYS TO CALCULATE                             
C     AE2,BE2,CE2     - THE EXCITATION ENERGY AFTER E.M. EXC.           
C     SR1,SR2,XR      - WITH MONTE CARLO                                
C---------------------------------------------------------------------  
C                                                                       
C     DEFORMATIONS AND G.S. SHELL EFFECTS                               
C     COMMON /ECLD/   ECGNZ,ECFNZ,VGSLD,ALPHA                           
C                                                                       
C     ECGNZ - GROUND STATE SHELL CORR. FRLDM FOR A SPHERICAL G.S.       
C     ECFNZ - SHELL CORRECTION FOR THE SADDLE POINT (NOW: == 0)         
C     VGSLD - DIFFERENCE BETWEEN DEFORMED G.S. AND LDM VALUE            
C     ALPHA - ALPHA GROUND STATE DEFORMATION (THIS IS NOT BETA2!)       
C             BETA2 = SQRT(5/(4PI)) * ALPHA                             
C---------------------------------------------------------------------  
C                                                                       
C     ARRAYS FOR EXCITATION ENERGY BY STATISTICAL HOLE ENERY MODEL      
C     COMMON /EENUC/  SHE, XHE                                          
C                                                                       
C     SHE, XHE - ARRAYS TO CALCULATE THE EXC. ENERGY AFTER              
C                ABRASION BY THE STATISTICAL HOLE ENERGY MODEL          
C---------------------------------------------------------------------  
C                                                                       
C     G.S. SHELL EFFECT                                                 
C     COMMON /EC2SUB/ ECNZ                                              
C                                                                       
C     ECNZ G.S. SHELL EFFECT FOR THE MASSES (IDENTICAL TO ECGNZ)        
C---------------------------------------------------------------------  
C                                                                       
C     OPTIONS AND PARAMETERS FOR FISSION CHANNEL                        
C     COMMON /FISS/    AKAP,BET,HOMEGA,KOEFF,IFIS,                       
C                            OPTSHP,OPTXFIS,OPTLES,OPTCOL               
C                                                                       
C     AKAP   - HBAR**2/(2* MN * R_0**2) = 10 MEV                        
C     BET    - REDUCED NUCLEAR FRICTION COEFFICIENT IN (10**21 S**-1)   
C     HOMEGA - CURVATURE OF THE FISSION BARRIER = 1 MEV                 
C     KOEFF  - COEFFICIENT FOR THE LD FISSION BARRIER == 1.0            
C     IFIS   - 0/1 FISSION CHANNEL OFF/ON                               
C     OPTSHP - INTEGER SWITCH FOR SHELL CORRECTION IN MASSES/ENERGY     
C            = 0 NO MICROSCOPIC CORRECTIONS IN MASSES AND ENERGY        
C            = 1 SHELL ,  NO PAIRING                                    
C            = 2 PAIRING, NO SHELL                                      
C            = 3 SHELL AND PAIRING                                      
C     OPTCOL - 0/1 COLLECTIVE ENHANCEMENT SWITCHED ON/OFF               
C     OPTXFIS- 0,1,2 FOR MYERS & SWIATECKI, DAHLINGER, ANDREYEV         
C              FISSILITY PARAMETER.                                     
C     OPTLES - CONSTANT TEMPERATURE LEVEL DENSITY FOR A,Z > TH-224      
C     OPTCOL - 0/1 COLLECTIVE ENHANCEMENT OFF/ON                        
C---------------------------------------------------------------------  
C                                                                       
C     OPTIONS                                                           
C     COMMON /OPT/    OPTEMD,OPTCHA,EEFAC                               
C                                                                       
C     OPTEMD - 0/1  NO EMD / INCL. EMD                                  
C     OPTCHA - 0/1  0 GDR / 1 HYPERGEOMETRICAL PREFRAGMENT-CHARGE-DIST. 
C              ***  RECOMMENDED IS OPTCHA = 1 ***                       
C     EEFAC  - EXCITATION ENERGY FACTOR, 2.0 RECOMMENDED                
C---------------------------------------------------------------------  
C                                                                       
C     FISSION BARRIERS                                                  
C     COMMON /FB/     EFA                                               
C     EFA    - ARRAY OF FISSION BARRIERS                                
C---------------------------------------------------------------------  
C                                                                       
C     LEVEL DENSITY PARAMETERS                                          
C     COMMON /ALD/    AV,AS,AK,OPTAFAN                                  
C     AV,AS,AK - VOLUME,SURFACE,CURVATURE DEPENDENCE OF THE             
C                LEVEL DENSITY PARAMETER                                
C     OPTAFAN - 0/1  AF/AN >=1 OR AF/AN ==1                             
C               RECOMMENDED IS OPTAFAN = 0                              
C---------------------------------------------------------------------  
C   ____________________________________________________________________
C  /                                                                    
C  /  INITIALIZES PARAMETERS IN COMMON /ABRAMAIN/, /EMDPAR/, /ECLD/ ... 
C  /  PROJECTILE PARAMETERS, EMD PARAMETERS, SHELL CORRECTION TABLES.   
C  /  CALCULATES MAXIMUM IMPACT PARAMETER FOR NUCLEAR COLLISIONS AND    
C  /  TOTAL GEOMETRICAL CROSS SECTION + EMD CROSS SECTIONS              
C   ____________________________________________________________________
C                                                                       
C                                                                       
      IMPLICIT NONE                                                     
      COMMON /ABLAMAIN/ AP,ZP,AT,ZT,EAP,BETA,BMAXNUC,CRTOT,CRNUC,R_0,   
     +                  R_P,R_T,IMAX,IRNDM,PI,BFPRO,SNPRO,SPPRO,SHELL   
      REAL*8 AP,ZP,AT,ZT,EAP,BETA,BMAXNUC,CRTOT,CRNUC,R_0,R_P,R_T,PI,   
     +       BFPRO,SNPRO,SPPRO,SHELL                                    
      INTEGER*4 IMAX,IRNDM,STATUS,SEED1,SEED2                           
C                                                                       
      COMMON /EMDPAR/ EGDR,EGQR,FWHMGDR,FWHMGQR,CREMDE1,CREMDE2,        
     &  AE1,BE1,CE1,AE2,BE2,CE2,SRE1,SRE2,XRE1,XRE2,DS1,DS2             
      INTEGER NQUAD                                                     
      PARAMETER(NQUAD=1000)                                             
      REAL*8 EGDR,EGQR,FWHMGDR,FWHMGQR,CREMDE1,CREMDE2                  
      REAL*8 AE1(NQUAD+1),BE1(NQUAD+1),CE1(NQUAD+1),AE2(NQUAD+1),       
     &  BE2(NQUAD+1),CE2(NQUAD+1),SRE1(NQUAD+1),SRE2(NQUAD+1),          
     &  XRE1(NQUAD+1),XRE2(NQUAD+1),DS1,DS2                             
C                                                                       
      REAL*8 DUQUAD,ELIM,EI,XR(NQUAD+1),SR1(NQUAD+1),SR2(NQUAD+1),      
     &       U1,U2                                                      
C                                                                       
      REAL*8 PR1(NQUAD+1),PR2(NQUAD+1)                                  
C                                                                       
      INTEGER I,NNUC,NEE,L,long                                              
C      PARAMETER(NEE=2001,NNUC=50)                                      
      REAL*8 XEE(1:2001),EEC(50,1:2001),SEE(50,1:2001),YE               
      REAL*8 SE1(1:2001), EPS, M                                        
C                                                                       
C XHE,SHE ARE COMMON TO BE USED IN REANUCL FOR THE LINEAR               
C INTERPOLATION                                                         
C                                                                       
      REAL*8 SHE(1:2001),XHE(50,1:2001)                                 
      COMMON /EENUC/ SHE,XHE                                            
C                                                                       
      COMMON /EC2SUB/ ECNZ                                              
      REAL*8 ECNZ                                                       
      DIMENSION ECNZ(0:153,0:98)                                        
C                                                                       
      INTEGER*4 Z,N,AP1,ZP1,AT1,ZT1,TSTAT                               
      REAL*8 R13,DNDE1,DNDE2,BP,SPMIN,SIGE1,SIGE2                       
      REAL*8 Y,BS,BK,AF,AN,BIPOL,SPDEF                                  
      INTEGER*4 INOW,OPTLES,OPTAFAN                                     
      DIMENSION INOW(1:2)                                               
      CHARACTER*30 DUM(25)                                              
      CHARACTER*50 NAME,FPATH                                           
      CHARACTER*80 TUPLE,ECSERN,VGSTAB,DEFTAB,FILENAME,RACINE,FILEDAT                           
      CHARACTER*100 FINA1,FINA2                                         
C                                                                       
      INTEGER*4 OPTEMD,OPTCHA,IFIS,OPTSHP,J,K,SEED,OPTXFIS,IO,OPTCOL    
      REAL*8 AKAP,BET,HOMEGA,KOEFF,MAZ,MA1Z,MA1Z1,EEFAC                 
      REAL*8 ECGNZ(0:153,0:98),ECFNZ(0:153,0:98),EFA(0:100,0:160),      
     &       VGSLD(0:153,0:98),ALPHA(0:153,0:98),AV,AK,AS,FISSILITY     
      COMMON /FISS/ AKAP,BET,HOMEGA,KOEFF,IFIS,                          
     &        OPTSHP,OPTXFIS,OPTLES,OPTCOL                              
      COMMON /ECLD/ ECGNZ,ECFNZ,VGSLD,ALPHA                             
      COMMON /OPT/ OPTEMD,OPTCHA,EEFAC                                  
      COMMON /FB/ EFA                                                   
      COMMON /ALD/AV,AS,AK,OPTAFAN                                      
C                                                                       
      REAL*4 SEGS,SELMAX,SBFIS                                          
      REAL*8 SEVAL                                                      
      PI = 3.141592653589793D0                                          
      R13 = 1.D0/3.D0                                                   
C                                                                       
C---------- SET INPUT VALUES                                            
C                                                                       
C *** INPUT FROM UNIT 10 IN THE FOLLOWING SEQUENCE !                    
C     AP1 =    INTEGER  !                                               
C     ZP1 =    INTEGER  !                                               
C     AT1 =    INTEGER  !                                               
C     ZT1 =    INTEGER  !                                               
C     EAP =    REAL     !                                               
C     IMAX =   INTEGER  !                                               
C     IFIS =   INTEGER SWITCH FOR FISSION                               
C     OPTSHP = INTEGER SWITCH FOR SHELL CORRECTION IN MASSES/ENERGY     
C            =0 NO MICROSCOPIC CORRECTIONS IN MASSES AND ENERGY         
C            =1 SHELL , NO PAIRING CORRECTION                           
C            =2 PAIRING, NO SHELL CORRECTION                            
C            =3 SHELL AND PAIRING CORRECTION IN MASSES AND ENERGY       
C     OPTEMD =0,1  0 NO EMD, 1 INCL. EMD                                
C               ELECTROMAGNETIC DISSOZIATION IS CALCULATED AS WELL.     
C     OPTCHA =0,1  0 GDR- , 1 HYPERGEOMETRICAL PREFRAGMENT-CHARGE-DIST. 
C               RECOMMENDED IS OPTCHA=1                                 
C     OPTCOL =0,1 COLLECTIVE ENHANCEMENT SWITCHED ON 1 OR OFF 0 IN DENSN
C     OPTAFAN=0,1 SWITCH FOR AF/AN = 1 IN DENSNIV 0 AF/AN>1 1 AF/AN=1   
C     AKAP =  REAL    ALWAYS EQUALS 10                                  
C     BET  =  REAL    REDUCED FRICTION COEFFICIENT / 10**(+21) S**(-1)  
C     HOMEGA = REAL   CURVATURE / MEV RECOMMENDED = 1. MEV              
C     KOEFF  = REAL   COEFFICIENT FOR FISSION BARRIER                   
C     OPTXFIS= INTEGER 0,1,2 FOR MYERS & SWIATECKI, DAHLINGER, ANDREYEV 
C              FISSILITY PARAMETER.                                     
C     EEFAC  = REAL EMPIRICAL FACTOR FOR THE EXCITATION ENERGY          
C                   RECOMMENDED 2.D0, STATISTICAL ABRASION MODELL 1.D0  
C     AV     = REAL KOEFFICIENTS FOR CALCULATION OF A(TILDE)            
C     AS     = REAL LEVEL DENSITY PARAMETER                             
C     AK     = REAL                                                     
C                                                                       
C This following inputs will be initialized in the main through the 
C         common /ABLAMAIN/  (A.B.)

C      READ(10,*,IOSTAT=IO) DUM(1),AP1                                   
C      READ(10,*,IOSTAT=IO) DUM(2),ZP1                                   
C      READ(10,*,IOSTAT=IO) DUM(3),AT1                                   
C      READ(10,*,IOSTAT=IO) DUM(4),ZT1                                   
C      READ(10,*,IOSTAT=IO) DUM(5),EAP                                   
C      READ(10,*,IOSTAT=IO) DUM(6),IMAX                                  
C      READ(10,*,IOSTAT=IO) DUM(7),IFIS                                  
C*** SWITCH-FISSION.1=ON.0=OFF
      IFIS = 1
C      READ(10,*,IOSTAT=IO) DUM(8),OPTSHP                                
C*** SHELL+PAIRING.0-1-2-3
      OPTSHP = 0
C      READ(10,*,IOSTAT=IO) DUM(9),OPTEMD                                
C*** OPTEMD =0,1  0 NO EMD, 1 INCL. EMD                                
      OPTEMD = 1
C      READ(10,*,IOSTAT=IO) DUM(10),OPTCHA                               
      OPTCHA = 1
C      READ(10,*,IOSTAT=IO) DUM(11),AKAP                                 
C*** NOT.TO.BE.CHANGED.(AKAP)
      AKAP = 10.
C      READ(10,*,IOSTAT=IO) DUM(12),BET                                  
C*** NUCLEAR.VISCOSITY.(BETA)
      BET = 1.5
C      READ(10,*,IOSTAT=IO) DUM(13),HOMEGA                               
C*** POTENTIAL-CURVATURE
      HOMEGA = 1.
C      READ(10,*,IOSTAT=IO) DUM(14),KOEFF                                
C*** FISSION-BARRIER-COEFFICIENT
      KOEFF = 1.
C      READ(10,*,IOSTAT=IO) DUM(15),OPTCOL                               
C*** COLLECTIVE ENHANCEMENT SWITCHED ON 1 OR OFF 0 IN DENSN (QR=val or =1.)
      OPTCOL = 0
C      READ(10,*,IOSTAT=IO) DUM(24),OPTLES                               
C*** SWITCH-FOR-LOW-ENERGY-SYS
      OPTLES = 0                                                                 
C      READ(10,*,IOSTAT=IO) DUM(16),EEFAC                                
      EEFAC = 2.
C      READ(10,*,IOSTAT=IO) DUM(20),OPTAFAN                              
      OPTAFAN = 0
C      READ(10,*,IOSTAT=IO) DUM(17),AV                                   
C      READ(10,*,IOSTAT=IO) DUM(18),AS                                   
C      READ(10,*,IOSTAT=IO) DUM(19),AK                                   
C*** LEVEL DENSITY PARAMETER
      AV = 0.073D0                                   
      AS = 0.095D0                                  
      AK = 0.0D0                                  
C      READ(10,FMT='(A80)',IOSTAT=IO) TUPLE                              
C      READ(10,FMT='(A80)',IOSTAT=IO) FINA1                              
C      READ(10,FMT='(A80)',IOSTAT=IO) FINA2                              
C      READ(10,FMT='(A80)',IOSTAT=IO) ECSERN                             
C      READ(10,FMT='(A80)',IOSTAT=IO) VGSTAB                             
C      READ(10,FMT='(A80)',IOSTAT=IO) DEFTAB                             
C      READ(10,*,IOSTAT=IO) DUM(23)                                      
C      STATUS=IO                                                         
C                                                                       
C------------ CONTROL OUTPUT TO LOG-FILE                                
C                                                                       
      WRITE(6,*) 'IFIS',IFIS                                            
      WRITE(6,*) 'OPTSHP',OPTSHP                                          
      WRITE(6,*) 'OPTEMD',OPTEMD                                          
      WRITE(6,*) 'OPTCHA',OPTCHA                                         
      WRITE(6,*) 'AKAP',AKAP                                           
      WRITE(6,*) 'BET',BET                                            
      WRITE(6,*) 'HOMEGA',HOMEGA                                         
      WRITE(6,*) 'KOEFF',KOEFF                                          
      WRITE(6,*) 'OPTCOL',OPTCOL                                         
      WRITE(6,*) 'OPTLES',OPTLES                                         
      WRITE(6,*) 'EEFAC',EEFAC                                          
      WRITE(6,*) 'OPTAFAN',OPTAFAN                                        
      WRITE(6,*) 'AV',AV                                             
      WRITE(6,*) 'AS',AS                                             
      WRITE(6,*) 'AK',AK                                             
C      WRITE(6,*) TUPLE(2:30),TUPLE(33:80)                               
C      WRITE(6,*) FINA1(2:30),FINA1(33:100)                              
C      WRITE(6,*) FINA1(2:30),FINA2(33:100)                              
C      WRITE(6,*) ECSERN(2:30),ECSERN(33:80)                             
C      WRITE(6,*) VGSTAB(2:30),VGSTAB(33:80)                             
C      WRITE(6,*) DEFTAB(2:30),DEFTAB(33:80)                             
C      TSTAT = 1                                                         
C      IF (TUPLE(33:36).EQ.'NONE') TSTAT = 0                             
C                                                                       
C----------- COPY IN REAL VARIABLES                                     
C                                                                       
C      AP = DFLOAT(AP1)                                                  
C      ZP = DFLOAT(ZP1)                                                  
C      AT = DFLOAT(AT1)                                                  
C      ZT = DFLOAT(ZT1)                                                  
C                                                                       
C      NAME  = FINA2(33:82)                                              
C      FPATH = FINA1(33:82)                                              
      OPTXFIS = 1                                                       
C                                                                       
C      IF (TSTAT.EQ.1) THEN                                              
C         OPEN(UNIT=21,STATUS='UNKNOWN',FILE=TUPLE(33:80))               
C      END IF                                                            
C                                                                       
C----------- READ SHELL CORRECTION TABLES                               
C                                                                       
C      OPEN(UNIT = 9,STATUS='UNKNOWN',SHARED,READONLY,FILE=ECSERN(33:80))
C                                                                       
C      OPEN(UNIT =11,STATUS='UNKNOWN',SHARED,READONLY,FILE=VGSTAB(33:80))
C                                                                       
C      OPEN(UNIT =13,STATUS='UNKNOWN',SHARED,READONLY,FILE=DEFTAB(33:80))
C                                                                       
C Verif par DIFF file1 file2: pas de differences dans les tables...(A.B.)                                          
      FILEDAT = 'frldm.tab'
      long=0
      i=0
      DO WHILE(long.EQ.0)
      	i=i+1
	IF(RACINE(i:i).EQ.' ') long=i-1
      END DO
      FILENAME=RACINE(1:long)//FILEDAT                                                         
      OPEN(UNIT = 9,STATUS='UNKNOWN',FILE=FILENAME)
C                                                                       
      FILEDAT = 'vgsld.tab'                                                   
      FILENAME=RACINE(1:long)//FILEDAT                                                                                                      
      OPEN(UNIT =11,STATUS='UNKNOWN',FILE=FILENAME)
C                                                                       
                                                 
      FILEDAT = 'flalpha.tab'                                                   
      FILENAME=RACINE(1:long)//FILEDAT                                                                                                      
      OPEN(UNIT =13,STATUS='UNKNOWN',FILE=FILENAME)    
C
C
C----------- INITIALIZE RANDOM NUMBERS                                  
C                                                                       
C
C--- This part is suppressed in WINDOWS version (KHS)
C                                                                       
c      CALL DATIME(INOW(1),INOW(2))                                      
C      IRNDM=INOW(2)                                                    
c      WRITE(6,*)'IRNDM=',IRNDM                                          
c      CALL RDMOUT(SEED)                                                 
c      WRITE(6,*)'INITIAL SEED ',SEED                                    
c      CALL RDMIN(INOW(2))                                               
c      CALL RDMOUT(SEED)                                                 
c      WRITE(6,*) 'SEED CHANGED ',SEED                                   
c      CALL RDMIN(12345)                                                 
c      CALL RDMOUT(SEED)                                                 
c      WRITE(6,*) 'SEED CHANGED ',SEED                                   
c      CALL NORRUT(SEED1,SEED2)                                          
c      WRITE(6,*) 'SEEDS FOR NORMAL DIST.',SEED1,SEED2                   
c      CALL NORRIN(12345,1073)                                           
c      CALL NORRUT(SEED1,SEED2)                                          
c      WRITE(6,*) 'SEEDS CHANGED',SEED1,SEED2                            
C                                                                       
C---------- INITIALIZE TABLES OF NUCLEAR SHELL EFFECTS.                 
C                                                                       
                                                                        
      DO 15  Z = 0,98,1                                                 
          READ(UNIT = 9,FMT=21)(ECNZ(N,Z),N=0,153)                      
 21       FORMAT (7D11.3)                                               
 15   CONTINUE                                                          
      CLOSE(UNIT = 9)                                                   
C                                                                       
      DO 32  Z = 0,98,1                                                 
        DO 33  N = 0,153,1                                              
             ECGNZ(N,Z) = ECNZ(N,Z)                                     
 33     CONTINUE                                                        
 32   CONTINUE                                                          
C                                                                       
      DO 16  Z = 0,98,1                                                 
          READ(UNIT = 11,FMT=21)(VGSLD(N,Z),N=0,153)                    
 16   CONTINUE                                                          
      CLOSE(UNIT = 11)                                                  
C                                                                       
      DO 17  Z = 0,98,1                                                 
          READ(UNIT = 13,FMT=21)(ALPHA(N,Z),N=0,153)                    
 17   CONTINUE                                                          
      CLOSE(UNIT = 13)                                                  
C                                                                       
      DO 30  Z = 0,98,1                                                 
        DO 31  N = 0,153,1                                              
         ECFNZ(N,Z) = 0.D0                                              
 31     CONTINUE                                                        
 30   CONTINUE                                                          
C                                                                       
      DO 73 K = 1,98,1                                                  
        DO 74 J = 1,153,1                                               
C                                                                       
C FISSION BARRIERS FROM SIERK (BARFIT) WITH 0 ANG. MOMENTUM             
C SEGS AND SELMAX ARE NOT 2001DED HERE.                                 
C A.J. 16.7.96                                                          
C           CALL BARFIT(K,(K+J),0,SBFIS,SEGS,SELMAX)                     
C          EFA(K,J) = DBLE(SBFIS) * KOEFF                               
C                                                                       
C ADD SHELL CORRECTION. PAIRING CANCELS.                                
C                                                                       
C          IF ((OPTSHP.EQ.1).OR.(OPTSHP.EQ.3)) THEN                     
C            EFA(K,J) = EFA(K,J) -  ECGNZ(J,K)                          
C          END IF                                                       
C                                                                       
C TO AVOID NEGATIVE VALUES FOR IMPOSSIBLE NUCLEI                        
C THE FISSION BARRIER IS SET TO ZERO IF SMALLER THAN ZERO.              
C                                                                       
C          IF (EFA(K,J).LT.0.D0) EFA(K,J)=0.D0                          
   74    CONTINUE                                                       
   73 CONTINUE                                                          
                                                                        

      RETURN                                                            
      END                                                               

C                                                                       
C******************************************************************     
C                                                                       
      DOUBLE PRECISION FUNCTION SPDEF(A,Z,OPTXFIS)                      
      INTEGER A,Z,OPTXFIS,INDEX                                         
      REAL*8 ALPHA2(36),X,FISSILITY,V,DX                                
C                                                                       
C INPUT:  A,Z,OPTXFIS MASS AND CHARGE OF A NUCLEUS,                     
C         OPTION FOR FISSILITY                                          
C OUTPUT: SPDEF                                                         
C                                                                       
C ALPHA2 SADDLE POINT DEF. COHEN&SWIATECKI ANN.PHYS. 22 (1963) 406      
C RANGING FROM FISSILITY X=0.30 TO X=1.00 IN STEPS OF 0.02              
C                                                                       
      DATA ALPHA2                                                       
     & / 2.5464D0, 2.4944D0, 2.4410D0, 2.3915D0, 2.3482D0, 2.3014D0,    
     &   2.2479D0, 2.1982D0, 2.1432D0, 2.0807D0, 2.0142D0, 1.9419D0,    
     &   1.8714D0, 1.8010D0, 1.7272D0, 1.6473D0, 1.5601D0, 1.4526D0,    
     &   1.3164D0, 1.1391D0, 0.9662D0, 0.8295D0, 0.7231D0, 0.6360D0,    
     &   0.5615D0, 0.4953D0, 0.4354D0, 0.3799D0, 0.3274D0, 0.2779D0,    
     &   0.2298D0, 0.1827D0, 0.1373D0, 0.0901D0, 0.0430D0, 0.0000D0/    
      DX = 0.02D0                                                       
      X  = FISSILITY(A,Z,OPTXFIS)                                       
      IF (X.GT.1.D0) X = 1.D0                                           
      IF (X.LT.0.D0) X = 0.D0                                           
      V  = (X-0.3D0)/DX + 1.D0                                          
      INDEX = IDNINT(V)                                                 
      IF (INDEX.LT.1) THEN                                              
         SPDEF = ALPHA2(1)                                              
         RETURN                                                         
      END IF                                                            
      IF (INDEX.EQ.36) THEN                                             
         SPDEF = ALPHA2(36)                                             
      ELSE                                                              
         SPDEF =  ALPHA2(INDEX) +                                       
     &    ( ALPHA2(INDEX+1) - ALPHA2(INDEX) ) / DX *                    
     &    ( X - (0.3D0 + DX*(INDEX-1) ) )                               
      END IF                                                            
      RETURN                                                            
      END                                                               
C                                                                       
C******************************************************************     
C                                                                       
      DOUBLE PRECISION FUNCTION FISSILITY(A,Z,OPTXFIS)                  
C                                                                       
C    CALCULATION OF FISSILITY PARAMETER                                 
C                                                                       
C    INPUT: A,Z INTEGER MASS & CHARGE OF NUCLEUS                        
C           OPTXFIS = 0 : MYERS, SWIATECKI                              
C                     1 : DAHLINGER                                     
C                     2 : ANDREYEV                                      
C                                                                       
      INTEGER A,Z,OPTXFIS                                               
      REAL*8 AA,ZZ,I                                                    
      AA = DFLOAT(A)                                                    
      ZZ = DFLOAT(Z)                                                    
      I  = DFLOAT(A-2*Z) / AA                                           
C                                                                       
C-------------  MYERS & SWIATECKI DROPLET MODELL                        
C                                                                       
      IF (OPTXFIS.EQ.0) THEN                                            
        FISSILITY = ZZ**2 / AA /50.8830D0 /                             
     &         (1.0D0 - 1.7826D0 * I**2)                                
      END IF                                                            
C                                                                       
      IF (OPTXFIS.EQ.1) THEN                                            
C                                                                       
C-------------- DAHLINGER FIT:                                          
C                                                                       
         FISSILITY = ZZ**2 / AA *                                       
     &   (49.22D0*(1.D0 - 0.3803D0*I**2 - 20.489D0*I**4))**(-1)         
      END IF                                                            
C                                                                       
       IF (OPTXFIS.EQ.2) THEN                                           
C                                                                       
C-------------- DUBNA FIT:                                              
C                                                                       
         FISSILITY = ZZ**2 / AA  /(48.D0*(1.D0 - 17.22D0*I**4))         
       END IF                                                           
C                                                                       
      RETURN                                                            
      END                                                               

C                                                                       
C********************************************************************** 
C                                                                       
      SUBROUTINE EVAPORA(ZPRF,APRF,EE,JPRF,ZF,AF,MTOTA,PLEVA,PXEVA,     
     &                   PYEVA,FF,INTTYPE,INUM)                               
C                                                                       
C     INPUT:                                                            
C                                                                       
C     ZPRF, APRF, EE(EE IS MODIFIED!), JPRF                             
C                                                                       
C     PROJECTILE AND TARGET PARAMETERS + CROSS SECTIONS                 
C     COMMON /ABRAMAIN/ AP,ZP,AT,ZT,EAP,BETA,BMAXNUC,CRTOT,CRNUC,       
C                       R_0,R_P,R_T, IMAX,IRNDM,PI,                     
C                       BFPRO,SNPRO,SPPRO,SHELL                         
C                                                                       
C     AP,ZP,AT,ZT   - PROJECTILE AND TARGET MASSES                      
C     EAP,BETA      - BEAM ENERGY PER NUCLEON, V/C                      
C     BMAXNUC       - MAX. IMPACT PARAMETER FOR NUCL. REAC.             
C     CRTOT,CRNUC   - TOTAL AND NUCLEAR REACTION CROSS SECTION          
C     R_0,R_P,R_T,  - RADIUS PARAMETER, PROJECTILE+ TARGET RADII        
C     IMAX,IRNDM,PI - MAXIMUM NUMBER OF EVENTS, DUMMY, 3.141...         
C     BFPRO         - FISSION BARRIER OF THE PROJECTILE                 
C     SNPRO         - NEUTRON SEPARATION ENERGY OF THE PROJECTILE       
C     SPPRO         - PROTON    "           "   "    "   "              
C     SHELL         - GROUND STATE SHELL CORRECTION                     
C                                                                       
C---------------------------------------------------------------------  
C     FISSION BARRIERS                                                  
C     COMMON /FB/     EFA                                               
C     EFA    - ARRAY OF FISSION BARRIERS                                
C---------------------------------------------------------------------  
C     OUTPUT:                                                           
C              ZF, AF, MTOTA, PLEVA, PTEVA, FF, INTTYPE, INUM           
C                                                                       
C     ZF,AF - CHARGE AND MASS OF FINAL FRAGMENT AFTER EVAPORATION       
C     MTOTA _ NUMBER OF EVAPORATED ALPHAS                               
C     PLEVA,PXEVA,PYEVA - MOMENTUM RECOIL BY EVAPORATION               
C     INTTYPE - TYPE OF REACTION 0/1 NUCLEAR OR ELECTROMAGNETIC         
C     FF      - 0/1 NO FISSION / FISSION EVENT                          
C     INUM    - EVENTNUMBER                                             
C   ____________________________________________________________________
C  /                                                                    
C  /  CALCUL DE LA MASSE ET CHARGE FINALES D'UNE CHAINE D'EVAPORATION   
C  /                                                                    
C  /  PROCEDURE FOR CALCULATING THE FINAL MASS AND CHARGE VALUES OF A   
C  /  SPECIFIC EVAPORATION CHAIN, STARTING POINT DEFINED BY (APRF, ZPRF,
C  /  EE)                                                               
C  /  On ajoute les 3 composantes de l'impulsion (PXEVA,PYEVA,PLEVA)
C  /    (actuellement PTEVA n'est pas correct; mauvaise norme...)                                                               
C  /____________________________________________________________________
C                                                                       
      SAVE                                                              
      COMMON /ABLAMAIN/ AP,ZP,AT,ZT,EAP,BETA,BMAXNUC,CRTOT,CRNUC,R_0,   
     +                  R_P,R_T,IMAX,IRNDM,PI,BFPRO,SNPRO,SPPRO,SHELL   
      REAL*8 AP,ZP,AT,ZT,EAP,BETA,BMAXNUC,CRTOT,CRNUC,R_0,R_P,R_T,PI,   
     &       BFPRO,SNPRO,SPPRO,SHELL                                   
      INTEGER*4 IMAX,IRNDM,INUM                                         
C                                                                       
      INTEGER*4 SORTIE                                                  
      REAL*8 ZF,AF,ZPRF,APRF,EPSILN,ALPHA,EE,PROBP,PROBN,PROBA,PTOTL,E  
      REAL*8 SN,SBP,SBA,MTOTA,X,AMOINS,ZMOINS,ECN,ECP,ECA,BP,BA         
      REAL*8 PLEVA,PTEVA,RNDX,RNDY,RNDZ,RNDN,JPRF,PXEVA,PYEVA                       
C                                                                       
      INTEGER*4 J,K,FF,INTTYPE,ITEST
      REAL*8 EFA,PROBF,EF,PC,MALPHA                                            
      REAL*8 CTET1,STET1,PHI1                                            
      REAL*4 RNDM,SBFIS,RND                                  
      REAL*8 ECGNZ(0:153,0:98),ECFNZ(0:153,0:98),
     &       VGSLD(0:153,0:98),ALPHAP(0:153,0:98)     
      COMMON /ECLD/ECGNZ,ECFNZ,VGSLD,ALPHAP                              

      INTEGER*4 IFIS,OPTSHP,OPTXFIS,OPTLES,OPTCOL    
      REAL*8 AKAP,BET,HOMEGA,KOEFF                 
      COMMON /FISS/AKAP,BET,HOMEGA,KOEFF,IFIS,OPTSHP,                    
     &        OPTXFIS,OPTLES,OPTCOL                                     
      COMMON /FB/ EFA(0:100,0:160)                                      

      DIMENSION IY(19)
C ial generateur pour le cascade (et les IY pour eviter les correlations)
       common/hazard/ial,IY

      real*4 acv,zpcv,pcv,xcv,ycv,zcv
      common/volant/acv(300),zpcv(300),pcv(300),xcv(300),
     s              ycv(300),zcv(300),iv 
C                                                                       
C-----------------------------------------------------------------------
C     IRNDM             DUMMY ARGUMENT FOR RANDOM-NUMBER FUNCTION       
C     SORTIE   LOCAL    HELP VARIABLE TO END THE EVAPORATION CHAIN      
C     ZF                NUCLEAR CHARGE OF THE FRAGMENT                  
C     ZPRF              NUCLEAR CHARGE OF THE PREFRAGMENT               
C     AF                MASS NUMBER OF THE FRAGMENT                     
C     APRF              MASS NUMBER OF THE PREFRAGMENT                  
C     EPSILN            ENERGY BURNED IN EACH EVAPORATION STEP          
C     MALPHA   LOCAL    MASS CONTRIBUTION TO MTOTA IN EACH EVAPORATION  
C                        STEP                                           
C     EE                EXCITATION ENERGY (VARIABLE)                    
C     PROBP             PROTON EMISSION PROBABILITY                     
C     PROBN             NEUTRON EMISSION PROBABILITY                    
C     PROBA             ALPHA-PARTICLE EMISSION PROBABILITY             
C     PTOTL             TOTAL EMISSION PROBABILITY                      
C     E                 LOWEST PARTICLE-THRESHOLD ENERGY                
C     SN                NEUTRON SEPARATION ENERGY                       
C     SBP               PROTON SEPARATION ENERGY PLUS EFFECTIVE COULOMB 
C                        BARRIER                                        
C     SBA               ALPHA-PARTICLE SEPARATION ENERGY PLUS EFFECTIVE 
C                        COULOMB BARRIER                                
C     BP                EFFECTIVE PROTON COULOMB BARRIER                
C     BA                EFFECTIVE ALPHA COULOMB BARRIER                 
C     MTOTA             TOTAL MASS OF THE EVAPORATED ALPHA PARTICLES    
C     X                 UNIFORM RANDOM NUMBER FOR NUCLEAR CHARGE        
C     AMOINS   LOCAL    MASS NUMBER OF EVAPORATED PARTICLE              
C     ZMOINS   LOCAL    NUCLEAR CHARGE OF EVAPORATED PARTICLE           
C     ECP               KINETIC ENERGY OF PROTON WITHOUT COULOMB        
C                        REPULSION                                      
C     ECN               KINETIC ENERGY OF NEUTRON                       
C     ECA               KINETIC ENERGY OF ALPHA PARTICLE WITHOUT COULOMB
C                        REPULSION                                      
C     PLEVA             TRANSVERSAL RECOIL MOMENTUM OF EVAPORATION      
C     PTEVA             LONGITUDINAL RECOIL MOMENTUM OF EVAPORATION     
C     FF                FISSION FLAG                                    
C     INTTYPE           INTERACTION TYPE FLAG                           
C     RNDX              RECOIL MOMENTUM IN X-DIRECTION IN A SINGLE STEP 
C     RNDY              RECOIL MOMENTUM IN Y-DIRECTION IN A SINGLE STEP 
C     RNDZ              RECOIL MOMENTUM IN Z-DIRECTION IN A SINGLE STEP 
C     RNDN              NORMALIZATION OF RECOIL MOMENTUM FOR EACH STEP  
C-----------------------------------------------------------------------
C                                                                       
      ZF = ZPRF                                                         
      AF = APRF
      PLEVA = 0.D0                                                      
      PTEVA = 0.D0                                                      
      PXEVA = 0.D0
      PYEVA = 0.D0

      SORTIE = 0                                                        
      FF = 0          

      itest = 0
c       if (ZF.eq.92) itest = 1 
      if (itest.eq.1) write(6,*) '***************************'
                                                  
 10   CONTINUE                                                          

      if (itest.eq.1) write(6,*) '------ZF,AF,EE------',
     &            IDNINT(ZF),IDNINT(AF),EE 

C      IF( (ZF.LE.0.D0).OR.(AF.LE.0.D0).OR.(AF.LT.ZF)) THEN             
C        WRITE(6,*) ZPRF,APRF                                           
C      END IF                                                           
C                                                                       
C CALCULATION OF THE PROBABILITIES FOR THE DIFFERENT DECAY CHANNELS     
C PLUS SEPARATION ENERGIES AND KINETIC ENERGIES OF THE PARTICLES        
C             
      CALL DIRECT(ZF,AF,EE,JPRF,PROBP,PROBN,PROBA,PROBF,PTOTL,          
     +            SN,SBP,SBA,ECN,ECP,ECA,BP,BA,INTTYPE,INUM,itest)  
C 
C modif A.B. incorporee OK!
      K = IDNINT(ZF)
      J = IDNINT(AF-ZF)

C     Now EF is calculated from EFA that depends on the subroutine
C     BARFIT which takes into account the modification on the ang. mom.
C     JB MVR 6-aug-1999
C     NOTE *** shell correction! (ECGNZ)  JB MVR 20-7-1999
C                                                                       
       IL  = IDNINT(JPRF)                                     
       CALL BARFIT(K,K+J,IL,SBFIS,SEGS,SELMAX)
       IF ((OPTSHP.EQ.1).OR.(OPTSHP.EQ.3)) THEN                     
            EFA(K,J) = DBLE(sbfis) -  ECGNZ(J,K)
C           efa(k,j) = DBLE(bfkhs(real(k),real(k+j))) - ECGNZ(j,k)
       ELSE                
            EFA(K,J) = DBLE(sbfis)
C           efa(k,j) = DBLE(bfkhs(real(k),real(k+j)))
       END IF 
      EF = EFA(K,J)                                                   
C                                                                       
C HERE THE FINAL STEPS OF THE EVAPORATION ARE CALCULATED                
C                                                                       
      IF ((SORTIE.EQ.1).OR.(PTOTL.EQ.0.D0)) THEN                        
C                                                                       
        E = DMIN1(SN,SBP,SBA)                                           
        IF (E.GT.1.D30) THEN                                            
         WRITE(6,*)'ERREUR A LA SORTIE EVAPORA,E>1.D30,AF=',AF,'ZF=',ZF 
        END IF                                                          
        IF (ZF.LE.6.) GO TO 100                                         
        IF (E.LT.0.D0) THEN                                             
          IF (SN.EQ.E) THEN                                             
            AF = AF-1.D0                                                
            ELSE IF (SBP.EQ.E) THEN                                     
              AF = AF-1.D0                                              
              ZF = ZF-1.D0                                              
            ELSE IF (SBA.EQ.E) THEN                                     
              AF = AF-4.D0                                              
              ZF = ZF-2.D0                                              
          END IF                                                        
          IF (AF.LT.2.5) GO TO 100                                      
          GO TO 10                                                      
        END IF                                                          
        GO TO 100                                                       
      END IF                                                            
 30   CONTINUE
      IRNDM = IRNDM+1                         
C                                                                       
C HERE THE NORMAL EVAPORATION CASCADE STARTS                            
C                                                                             
C RANDOM NUMBER FOR THE EVAPORATION                                     
C                                                                       
      X = DBLE(RNDM(IRNDM))*PTOTL                                       
C                                                                       
      IF (X.LT.PROBA) THEN                                              

C                                                                       
C ALPHA EVAPORATION                                                     
C                                                                       
        if (itest.eq.1) write(6,*) '< alpha evaporation >'
        AMOINS = 4.D0                                                   
        ZMOINS = 2.D0                                                   
        EPSILN = SBA+ECA
        PC = DSQRT((1.D0+(ECA+BA)/3.72834D3)**2-1.D0) * 3.72834D3       
        MALPHA = 4.D0                                                    
C Volant:
	iv = iv + 1
	acv(iv) = 4.
	zpcv(iv) = 2.
	pcv(iv) = PC

         ELSE IF (X.LT.PROBA+PROBP) THEN                                 
C                                                                       
C PROTON EVAPORATION                                                    
C                                                                       
          if (itest.eq.1) write(6,*) '< proton evaporation >'
          AMOINS = 1.D0                                                 
          ZMOINS = 1.D0                                                 
          EPSILN = SBP+ECP                                              
          PC = DSQRT((1.D0+(ECP+BP)/9.3827D2)**2-1.D0) * 9.3827D2       
          MALPHA = 0.D0                                                  
C Volant:
	iv = iv + 1
	acv(iv) = 1.
	zpcv(iv) = 1.
	pcv(iv) = PC
	                                                    
        ELSEIF (X.LT.PROBA+PROBP+PROBN) THEN                            
C                                                                       
C NEUTRON EVAPORATION                                                   
C                                                                       
          if (itest.eq.1) write(6,*) '< neutron evaporation >'
          AMOINS = 1.D0                                                 
          ZMOINS = 0.D0                                                 
          EPSILN = SN+ECN           
          PC = DSQRT((1.D0+(ECN)/9.3956D2)**2-1.D0) * 9.3956D2          
          MALPHA = 0.D0                                                  
C Volant:
	iv = iv + 1
	acv(iv) = 1.
	zpcv(iv) = 0.
	pcv(iv) = PC
	                                                    
        ELSE                                                            
C                                                                       
C FISSION                                                               
C                                                                       
C IN CASE OF FISSION-EVENTS THE FRAGMENT NUCLEUS IS THE MOTHER NUCLEUS  
C BEFORE FISSION OCCURS WITH EXCITATION ENERGY ABOVE THE FIS.- BARRIER. 
C FISSION FRAGMENT MASS DISTRIBUTION IS CALULATED IN SUBROUTINE FISDIS  
C                                                                       
          if (itest.eq.1) write(6,*) '< fission >'
          AMOINS = 0.D0                                                 
          ZMOINS = 0.D0                                                 
          EPSILN = EF                                                   
C                                                                       
          MALPHA = 0.D0                                                  
          PC = 0.D0                                                     
          FF = 1                                                        
C          IF (INTTYPE.EQ.1) THEN                                       
C             ITEC = ITEC + 1                                           
C             WRITE(6,*) 'X,EE,PROBF,PROBN,ITEC',X,EE,PROBF,PROBN,      
C     &       ITEC                                                      
C          END IF                                                       
      END IF                                                            
C      IF (EE.LT.EPSILN) THEN                                           
C        IF ((EE.GE.SBA+ECA).AND.(PROBA.GT.1.D-3)) GO TO 30             
C        IF ((EE.GE.SBP+ECP).AND.(PROBP.GT.1.D-3)) GO TO 30             
C        IF ((EE.GE.SN+ECN).AND.(PROBN.GT.1.D-3)) GO TO 30              
C        IF ((EE.GE.EF).AND.(PROBF.GT.1.D-3)) GO TO 30                  
C        SORTIE = 1                                                     
C        GO TO 10                                                       
C      END IF                                                           
C                                                                       

      if (itest.EQ.1) then
        write(6,*)'SN,SBP,SBA,EF',SN,SBP,SBA,EF
        write(6,*)'PROBN,PROBP,PROBA,PROBF,PTOTL'
     &            ,PROBN,PROBP,PROBA,PROBF,PTOTL 
      endif

C CALCULATION OF THE DAUGHTER NUCLEUS                                   
C                                                                       
      AF = AF-AMOINS                                                    
      ZF = ZF-ZMOINS                                                
      EE = EE-EPSILN                                                    
      IF (EE.LE.0.01D0) EE = 0.01E0                                     
      MTOTA = MTOTA+MALPHA                                              
C                                                                       
C La methode ci dessous est fausse --> non spherique!                                                                       
C      IRNDM = IRNDM+1                                                   
C      RNDX = DBLE(RNDM(IRNDM)) - 0.5D0                                  
C      IRNDM = IRNDM+1                                                   
C      RNDY = DBLE(RNDM(IRNDM)) - 0.5D0                                  
C      IRNDM = IRNDM+1                                                   
C      RNDZ = DBLE(RNDM(IRNDM)) - 0.5D0                                  
C      RNDN = DSQRT(RNDX**2+RNDY**2+RNDZ**2)                             
C----> (ancien)      PTEVA = PTEVA + PC * RNDX/RNDN

	IF(FF.EQ.0) THEN
           CALL RIBM(RND,IY(9))
           CTET1 = 2.*RND-1.
           CALL RIBM(RND,IY(5))
           PHI1 = RND*2.*3.141592654
	   STET1 = DSQRT(1.-CTET1**2)
      xcv(iv) =  STET1*DCOS(PHI1)
      ycv(iv) =  STET1*DSIN(PHI1)
      zcv(iv) =  CTET1                                 
      PXEVA = PXEVA - PC * xcv(iv)
      PYEVA = PYEVA - PC * ycv(iv)
      PLEVA = PLEVA - PC * CTET1                                    
	END IF
C                                                                       
C     WRITE(6,*)'-- PC=',PC                                             
C     WRITE(6,*)'   PLEVA=',PLEVA,' PTEVA=',PTEVA                       
C                                                                       
C                                                                       
C    CONDITION FOR END OF EVAPORATION                                   
C                                                                       
      IF ((AF.LT.2.5D0).OR.(FF.EQ.1)) GO TO 100                         
      GO TO 10                                                          
 100  CONTINUE                                                          
C                                                                       
      RETURN                                                            
      END                                                               
C                                                                       
C********************************************************************** 
C                                                                       
C=====CORRECTED BY KUDYAEV==============================================
      SUBROUTINE DIRECT(ZPRF,A,EE,JPRF,PROBP,PROBN,PROBA,PROBF,PTOTL,   
     +		SN,SBP,SBA,ECN,ECP,ECA,BP,BA,INTTYPE,INUM,itest)        
C=======================================================================
C                                                                       
C     ________________________________________________________________  
C     /                                                               / 
C     /  CALCULATION OF PARTICLE-EMISSION PROBABILITIES & FISSION     / 
C     /  BASED ON THE SIMPLIFIED FORMULAS FOR THE DECAY WIDTH BY      / 
C     /  MORETTO, ROCHESTER MEETING TO AVOID COMPUTING TIME           / 
C     /  INTENSIVE INTEGRATION OF THE LEVEL DENSITIES                 / 
C     /  USES EFFECTIVE COULOMB BARRIERS AND AN AVERAGE KINETIC ENERGY/ 
C     /  OF THE EVAPORATED PARTICLES                                  / 
C     /  COLLECTIVE ENHANCMENT OF THE LEVEL DENSITY IS INCLUDED       / 
C     /  DYNAMICAL HINDRANCE OF FISSION IS INCLUDED BY A STEP FUNCTION/ 
C     /  APPROXIMATION. SEE A.R. JUNGHANS DIPLOMA THESIS              / 
C     /  SHELL AND PAIRING STRUCTURES IN THE LEVEL DENSITY IS INCLUDED/ 
C     /_______________________________________________________________/ 
C                                                                       
C     INPUT:                                                            
C            ZPRF,A,EE  CHARGE, MASS, EXCITATION ENERGY OF COMPOUND     
C                       NUCLEUS                                         
C            JPRF       ROOT-MEAN-SQUARED ANGULAR MOMENTUM                           
C                                                                       
C     DEFORMATIONS AND G.S. SHELL EFFECTS                               
C     COMMON /ECLD/   ECGNZ,ECFNZ,VGSLD,ALPHA                           
C                                                                       
C     ECGNZ - GROUND STATE SHELL CORR. FRLDM FOR A SPHERICAL G.S.       
C     ECFNZ - SHELL CORRECTION FOR THE SADDLE POINT (NOW: == 0)         
C     VGSLD - DIFFERENCE BETWEEN DEFORMED G.S. AND LDM VALUE            
C     ALPHA - ALPHA GROUND STATE DEFORMATION (THIS IS NOT BETA2!)       
C             BETA2 = SQRT((4PI)/5) * ALPHA                             
C                                                                       
C     OPTIONS AND PARAMETERS FOR FISSION CHANNEL                        
C     COMMON /FISS/    AKAP,BET,HOMEGA,KOEFF,IFIS,                       
C                            OPTSHP,OPTXFIS,OPTLES,OPTCOL               
C---------------------------------------------------------------------  
C                                                                       
C     AKAP   - HBAR**2/(2* MN * R_0**2) = 10 MEV, R_0 = 1.4 FM          
C     BET    - REDUCED NUCLEAR FRICTION COEFFICIENT IN (10**21 S**-1)   
C     HOMEGA - CURVATURE OF THE FISSION BARRIER = 1 MEV                 
C     KOEFF  - COEFFICIENT FOR THE LD FISSION BARRIER == 1.0            
C     IFIS   - 0/1 FISSION CHANNEL OFF/ON                               
C     OPTSHP - INTEGER SWITCH FOR SHELL CORRECTION IN MASSES/ENERGY     
C            = 0 NO MICROSCOPIC CORRECTIONS IN MASSES AND ENERGY        
C            = 1 SHELL ,  NO PAIRING                                    
C            = 2 PAIRING, NO SHELL                                      
C            = 3 SHELL AND PAIRING                                      
C     OPTCOL - 0/1 COLLECTIVE ENHANCEMENT SWITCHED ON/OFF               
C     OPTXFIS- 0,1,2 FOR MYERS & SWIATECKI, DAHLINGER, ANDREYEV         
C              FISSILITY PARAMETER.                                     
C     OPTLES - CONSTANT TEMPERATURE LEVEL DENSITY FOR A,Z > TH-224      
C     OPTCOL - 0/1 COLLECTIVE ENHANCEMENT OFF/ON                        
C---------------------------------------------------------------------  
C     LEVEL DENSITY PARAMETERS                                          
C     COMMON /ALD/    AV,AS,AK,OPTAFAN                                  
C     AV,AS,AK - VOLUME,SURFACE,CURVATURE DEPENDENCE OF THE             
C                LEVEL DENSITY PARAMETER                                
C     OPTAFAN - 0/1  AF/AN >=1 OR AF/AN ==1                             
C               RECOMMENDED IS OPTAFAN = 0                              
C---------------------------------------------------------------------  
C     FISSION BARRIERS                                                  
C     COMMON /FB/     EFA                                               
C     EFA    - ARRAY OF FISSION BARRIERS                                
C---------------------------------------------------------------------  
C                                                                       
C                                                                       
C     OUTPUT: PROBN,PROBP,PROBA,PROBF,PTOTL:                            
C           - EMISSION PROBABILITIES FOR N EUTRON, P ROTON,  A LPHA     
C             PARTICLES, F ISSION AND NORMALISATION                     
C             SN,SBP,SBA: SEPARATION ENERGIES N P A                     
C             INCLUDING EFFECTIVE BARRIERS                              
C             ECN,ECP,ECA,BP,BA                                         
C           - AVERAGE KINETIC ENERGIES (2*T) AND EFFECTIVE BARRIERS     
C                                                                       
C                                                                       
      SAVE                                                              
      INTEGER*4  REFMOD,K,J,IFIS,OPTSHP,AFP,IZ,IN,OPTXFIS,OPTLES,       
     &      INTTYPE,OPTCOL,OPTAFAN,INUM,ILAST,imaxwell,itest
      REAL*4 RAT,RPT,RNT                           
      REAL*8 A,ZPRF,EE,MAZ,MA1Z,MA1Z1,MA4Z2,SN,SP,SA,BP,SBP,SBA         
      REAL*8 DENSN,DENSP,DENSA,RN,RA,RP,BA,PROBP,PROBN,PROBA,PTOTL      
      REAL*8 DENOMI,ECP,ECN,ECA,GNGF,TCONST,DCONST,NPRF,PARC            
      REAL*8 EF,SBF,DENSF,RF,GF,PROBF,EFA,KOEFF,TAUC,EER                
      REAL*8 AKAP,CRAM,TAU,BET,HOMEGA,PI,DENSG,GP,GN,GA,TS1,TSUM        
      REAL*8 CF,WF,TEMP,AT,PT,NT,FT,BSHELL,WFEX,GSUM                
      REAL*8 XX,Y,BS,BK,BIPOL,AV,AS,AK,FISSILITY,EDYN,HBAR              
      REAL*8 ECGNZ(0:153,0:98),ECFNZ(0:153,0:98),VGSLD(0:153,0:98)      
      REAL*8 ALPHA(0:153,0:98), DEFBET,JPRF,SPDEF                       
      COMMON /ECLD/ECGNZ,ECFNZ,VGSLD,ALPHA                              
   
      COMMON /FISS/AKAP,BET,HOMEGA,KOEFF,IFIS,OPTSHP,                    
     &        OPTXFIS,OPTLES,OPTCOL
                                          
      COMMON /ALD/AV,AS,AK,OPTAFAN                                      
      COMMON /FB/EFA(0:100,0:160)                                       
      PARAMETER(PI=3.1415926535D0,HBAR=6.582122D-22)                    

      REAL*4 RND
      DIMENSION IY(19)
C ial generateur pour le cascade (et les IY pour eviter les correlations)
       common/hazard/ial,IY

C
C     Switch to calculate Maxwellian distribution of kinetic energies (1=Max)                                                  
      imaxwell = 1                           
                   
C
C LIMITING OF EXCITATION ENERGY WHERE FISSION OCCURS                    
C !!! THIS IS NOT THE DYNAMICAL HINDRANCE (SEE END OF ROUTINE) !!!      
C                                                                       
      EDYN = 1000.D0                                                    
C                                                                       
C NO LIMIT IF STATISTICAL MODEL IS CALCULATED.                         
C                                                                       
      IF (BET.LE.1.D-16) EDYN = 10000.D0                                
C                                                                       
C JUST A CHANGE OF NAME UNTIL THE END OF THIS SUBROUTINE                
C                                                                       
      EER = EE                                                          
      IF (INUM.EQ.1) ILAST = 1                                          
C                                                                       
C                                                                       
C-------CALCULATION OF MASSES                                           
C                                                                       
C     REFMOD = 1 ==> MYERS,SWIATECKI MODEL                              
C     REFMOD = 0 ==> WEIZSAECKER MODEL                                  
C                                                                       
      REFMOD = 1                                                        
C                                                                       
      IF (REFMOD.EQ.1) THEN                                             
        CALL MGLMS(A,ZPRF,OPTSHP,MAZ)                                   
        CALL MGLMS(A-1.0,ZPRF,OPTSHP,MA1Z)                              
        CALL MGLMS(A-1.0,ZPRF-1.0,OPTSHP,MA1Z1)                         
        CALL MGLMS(A-4.0,ZPRF-2.0,OPTSHP,MA4Z2)                         
      ELSE                                                              
        CALL MGLW(A,ZPRF,OPTSHP,MAZ)                                    
        CALL MGLW(A-1.0,ZPRF,OPTSHP,MA1Z)                               
        CALL MGLW(A-1.0,ZPRF-1.0,OPTSHP,MA1Z1)                          
        CALL MGLW(A-4.0,ZPRF-2.0,OPTSHP,MA4Z2)                          
      END IF                                                            
C                                                                       
C--------SEPARATION ENERGIES AND EFFECTIVE BARRIERS                     
C                                                                       
      SN = MA1Z-MAZ                                                     
      SP = MA1Z1-MAZ                                                    
      SA = MA4Z2-MAZ-28.29688D0             
      IF (ZPRF.LT.1.D0) THEN                                            
         SBP = 1.D75                                                    
         GO TO 30                                                       
      END IF                                                            
C     Parameterisation GAIMARD:
C     BP = 1.44D0*(ZPRF-1.D0)/(1.22D0*(A-1.D0)**(1.D0/3.D0)+5.6D0)     
C     Parameterisation KHS (12-99)
      BP = 1.44D0*(ZPRF-1.D0)/(2.1D0*(A-1.D0)**(1.D0/3.D0)+0.0D0)
C                                                                       
C coulomb barrier manipulation  ONLY FOR TESTING                        
C                                                                       
C      IF (BP.GT.1.D0) THEN                                             
C         BP = BP+2.D0                                                  
C      ELSE                                                             
C         BP = 0.0001D0                                                 
C      END IF                                                           
C                                                                       
      SBP = SP+BP                                                       
      IF (A-4.0.LE.0.0) THEN                                            
        SBA = 1.D+75                                                    
        GO TO 30                                                        
      END IF                                                   
c*** New effective barrier for alpha evaporation d=6.1: KHS          
C     BA = 2.88D0*(ZPRF-2.D0)/(1.22D0*(A-4.D0)**(1.D0/3.D0)+6.1D0)
C     Parametrisation KHS (12-99)
C     BA =2.88D0*(ZPRF-2.D0)/(2.3D0*(A-4.D0)**(1.D0/3.D0)+0.0D0)
      BA =2.88D0*(ZPRF-2.D0)/(2.2D0*(A-4.D0)**(1.D0/3.D0)+0.0D0)
C                                                                       
C coulomb barrier manipulation ONLY FOR TESTING                         
C                                                                       
C      IF (BA.GT.1.D0) THEN                                             
C         BA = BA-1.D0                                                  
C      ELSE                                                             
C         BA = 0.0001D0                                                 
C      END IF                                                           
C                                                                       
      SBA = SA+BA                                                       
 30   CONTINUE                                                          
C                                                                       
C---------CALCULATION OF SURFACE AND CURVATURE INTEGRALS NEEDED TO      
C         TO CALCULATE THE LEVEL DENSITY PARAMETER (IN DENSNIV)         
C                                                                       
      IF (IFIS.GT.0) THEN                                               
        K  = IDNINT(ZPRF)                                               
        J  = IDNINT(A - ZPRF)                                           
C
C     Now EF is calculated from EFA that depends on the subroutine
C     BARFIT which takes into account the modification on the ang. mom.
C     JB MVR 6-aug-1999
C     NOTE *** shell correction! (ECGNZ)  JB MVR 20-7-1999
C                                                                       
        IL  = IDNINT(JPRF)                                     
        CALL BARFIT(K,K+J,IL,SBFIS,SEGS,SELMAX)
          IF ((OPTSHP.EQ.1).OR.(OPTSHP.EQ.3)) THEN                     
            EFA(K,J) = DBLE(sbfis) -  ECGNZ(J,K)
C           efa(k,j) = DBLE(bfkhs(real(k),real(k+j)))-ecgnz(j,k)
          ELSE                
            EFA(K,J) = DBLE(sbfis)
C           efa(k,j) = DBLE(bfkhs(real(k),real(k+j)))
          END IF 
        EF = EFA(K,J)                                                   
C
C TO AVOID NEGATIVE VALUES FOR IMPOSSIBLE NUCLEI                        
C THE FISSION BARRIER IS SET TO ZERO IF SMALLER THAN ZERO.              
C                                                                       
          IF (EFA(K,J).LT.0.D0) EFA(K,J)=0.D0                          
C                                                                       
C  FACTOR WITH JPRF SHOULD BE 0.0025D0 - 0.01D0 FOR                     
C  APPROXIMATE INFLUENCE OF ANG. MOMENTUM ON BFIS  A.J. 22.07.96        
C  0.D0 MEANS NO ANGULAR MOMENTUM                                       
C                                                                       
C        EF = (EF - 0.0D0*JPRF)                                          
C                                                                       
        IF (EF.LT.0.D0) EF = 0.D0                                         
        XX = FISSILITY((K+J),K,OPTXFIS)                                 
        Y=1.0D0-XX                                                      
        IF (Y.LT.0.D0) Y = 0.D0                                         
        IF (Y.GT.1.D0) Y = 1.D0                                         
        BS = BIPOL(1,Y)                                                 
        BK = BIPOL(2,Y)                                                 
      ELSE                                                              
        EF=1.0D40                                                       
        BS = 1.D0                                                       
        BK = 1.D0                                                       
      END IF                                                            
      SBF=EE-EF                                                         
C     WRITE(6,*)' EF= ',EF,' Z= ',ZZPRF,' N=',ZNPRF, ' EE=',EE          
      AFP = IDNINT(A)                                                   
      IZ = IDNINT(ZPRF)                                                 
      IN = AFP - IZ                                                     
      BSHELL = ECFNZ(IN,IZ)                                             
C                                                                       
C CONST. SADDLE POINT DEF. USED BY MYERS&SWIATECKI                      
C                                                                       
C      DEFBET    = 1.5D0 * 0.52D0                                       
C                                                                       
C  LD SADDLE POINT DEFORMATION                                          
C                                                                       
C HERE: BETA2 = SQRT(5/(4PI)) * ALPHA2                                  
C                                                                       
C FOR THE GROUND STATE DEF. 1.5D0 SHOULD BE USED                        
C BECAUSE THIS WAS JUST THE FACTOR TO PRODUCE THE                       
C ALPHA-DEFORMATION TABLE 1.5D0 SHOULD BE USED                          
C A.R.J. 6.8.97                                                         
C                                                                       
      DEFBET = 1.58533D0 * SPDEF(IDNINT(A),IDNINT(ZPRF),OPTXFIS)        
C                                                                       
C LEVEL DENSITY AND TEMPERATURE AT THE SADDLE POINT                     
C                                                                       
      CALL DENSNIV(A,ZPRF,EE,EF,DENSF,BSHELL,BS,BK,TEMP,                
     &             OPTSHP,OPTCOL,DEFBET)                                
      FT=TEMP                                                           
      IF (IZ.GE.2) THEN                                                 
       BSHELL = ECGNZ(IN,IZ-1) - VGSLD(IN,IZ-1)                         
       DEFBET    = 1.5D0 *  ALPHA(IN,IZ-1)                              
C                                                                       
C LEVEL DENSITY AND TEMPERATURE IN THE PROTON DAUGHTER                  
C                                                                       
       CALL DENSNIV(A-1.0,ZPRF-1.0D0,EE,SBP,DENSP,                      
     &	            BSHELL,1.D0,1.D0,TEMP,OPTSHP,OPTCOL,DEFBET)          
C       
       PT = TEMP 
       if (imaxwell.eq.1) then                                                       
c*** Valentina - random kinetic energy in a Maxwelliam distribution
C Modif Juin/2002 A.B. C.V. for light targets; limit on the energy
C   from the Maxwell distribution.
        RPT = PT
	ECP=2.D0 * PT 
	IF(RPT.LE.1.E-3) GOTO 2914
	    IFLAG=0
1914        ECP = FMAXHAZ(KKK,RPT)
	    IFLAG=IFLAG+1
	    IF(IFLAG.GE.10) THEN
		call ribm(RND,IY(6))
		ECP=SQRT(RND)*(EER-SBP)
	    	GOTO 2914
	    ENDIF
	    IF((ECP+SBP).GT.EER) GO TO 1914
       else
        ECP = 2.D0 * PT 
       endif                                                 
2914   CONTINUE
      ELSE                                                              
       DENSP = 0.D0                                                     
       ECP = 0.D0                                                       
       PT = 0.D0                                                        
      END IF 
      
C      IF(ECP.GE.30.)
C     s WRITE(6,*) 'SP,SBP,ECP,E*,A',SP,SBP,ECP,EER,A
                                                                 
C                                                                       
      IF (IN.GE.2) THEN                                                 
       BSHELL = ECGNZ(IN-1,IZ) - VGSLD(IN-1,IZ)                         
       DEFBET = 1.5D0 * ALPHA(IN-1,IZ)                                  
C                                                                       
C LEVEL DENSITY AND TEMPERATURE IN THE NEUTRON DAUGHTER                 
C                                                                       
       CALL DENSNIV(A-1.0,ZPRF,EE,SN,DENSN,BSHELL,                      
     &		    1.D0,1.D0,TEMP,OPTSHP,OPTCOL,DEFBET)                        
       NT = TEMP                                                        
       if (imaxwell.eq.1) then
c*** Valentina - random kinetic energy in a Maxwelliam distribution
C Modif Juin/2002 A.B. C.V. for light targets; limit on the energy
C   from the Maxwell distribution.
        RNT = NT
	ECN=2.D0 * NT
	IF(RNT.LE.1.E-3) GOTO 2915
	    IFLAG=0
1915        ECN = FMAXHAZ(KKK,RNT)
	    IFLAG=IFLAG+1
	    IF(IFLAG.GE.10) THEN
		call ribm(RND,IY(7))
		ECN=SQRT(RND)*(EER-SN)
	    	GOTO 2915
	    ENDIF            
	    IF((ECN+SN).GT.EER) GO TO 1915
       else                                            
        ECN = 2.D0 * NT
       endif                                                  
2915   CONTINUE
      ELSE                                                              
       DENSN = 0.D0                                                     
       ECN = 0.D0                                                       
       NT = 0.D0                                                        
      END IF                                                            
C                                                                       
      IF ((IN.GE.3).AND.(IZ.GE.3)) THEN                                 
       BSHELL = ECGNZ(IN-2,IZ-2) - VGSLD(IN-2,IZ-2)                     
       DEFBET = 1.5D0 * ALPHA(IN-2,IZ-2)                                
C                                                                       
C LEVEL DENSITY AND TEMPERATURE IN THE ALPHA DAUGHTER                   
C                                                                       
       CALL  DENSNIV(A-4.0,ZPRF-2.0D0,EE,SBA,DENSA,                     
     &		    BSHELL,1.D0,1.D0,TEMP,OPTSHP,OPTCOL,DEFBET)                 
c*** Valentina - random kinetic energy in a Maxwelliam distribution
       AT = TEMP   
       if (imaxwell.eq.1) then                                                     
C Modif Juin/2002 A.B. C.V. for light targets; limit on the energy
C   from the Maxwell distribution.
        RAT = AT
	ECA= 2.D0 * AT
	IF(RAT.LE.1.E-3) GOTO 2916 
	    IFLAG=0
1916        ECA = FMAXHAZ(KKK,RAT)
	    IFLAG=IFLAG+1
	    IF(IFLAG.GE.10) THEN
		call ribm(RND,IY(8))
		ECA=SQRT(RND)*(EER-SBA)
	    	GOTO 2916
	    ENDIF            
	    IF((ECA+SBA).GT.EER) GO TO 1916
       else                                                 
        ECA = 2.D0 * AT
       endif         
2916   CONTINUE
      ELSE                                                              
       DENSA = 0.D0                                                     
       ECA = 0.D0                                                       
       AT = 0.D0                                                        
      END IF                                                            
C                                                                       
C SPECIAL TREATMENT FOR UNBOUND NUCLEI                                                
C                                                                       
      IF (SN.LT.0.D0) THEN                                              
         PROBN = 1.D0                                                   
         PROBP = 0.D0                                                   
         PROBA = 0.D0                                                   
         PROBF = 0.D0                                                   
         GO TO 70                                                       
      END IF                                                            
      IF (SBP.LT.0.D0) THEN                                             
         PROBP = 1.D0                                                   
         PROBN = 0.D0                                                   
         PROBA = 0.D0                                                   
         PROBF = 0.D0                                                   
         GO TO 70                                                       
      END IF                                                            
C                                                                       
C  NO FISSION IF E*> EDYN                                               
C  OR MASS < 50                                                         
C                                                                       
      IF ((A.LT.50.D0).OR.(EE.GT.EDYN)) DENSF = 0.D0                    
C                                                                       
C                                                                       
      BSHELL = ECGNZ(IN,IZ) - VGSLD(IN,IZ)                              
      DEFBET = 1.5D0 * ALPHA(IN,IZ)                                     
C                                                                       
C COMPOUND NUCLEUS LEVEL DENSITY                                        
C                                                                       
      CALL DENSNIV(A,ZPRF,EE,0.0D0,DENSG,BSHELL,                        
     &	     1.D0,1.D0,TEMP,OPTSHP,OPTCOL,DEFBET)                        
C                                                                       
      IF ( DENSG.GT.0.D0) THEN                                          
C                                                                       
C CALCULATION OF THE PARTIAL DECAY WIDTH                                
C USED FOR BOTH THE TIME SCALE AND THE EVAPORATION DECAY WIDTH
C                                                                       
       GP = (A**(2./3.)/AKAP)*DENSP/DENSG/PI*PT**2                      
       GN = (A**(2./3.)/AKAP)*DENSN/DENSG/PI*NT**2                      
       GA = (A**(2./3.)/AKAP)*DENSA/DENSG/PI*2.0D0*AT**2                
       GF= DENSF/DENSG/PI/2.0D0*FT                                      

      if (itest.eq.1) write(6,*) 'GN,GP,GA,GF',GN,GP,GA,GF
C                                                                       
      ELSE                                                              
       WRITE(6,*) 'DIRECT: DENSG <= 0.D0',A,ZPRF,EE                     
      END IF                                                            
C                                                                       

      GSUM = GA + GP + GN                                           
      IF (GSUM.GT.0.D0) THEN                                        
        TS1  = HBAR / GSUM                                          
      ELSE                                                          
          TS1  = 1.D99                                                
      END IF                                                        

C                                                                       
C NEW EVENT MEANS RESET THE TIME SCALE           
C                                                                       
      IF (INUM.GT.ILAST) TSUM = 0                                         

C       WRITE(6,'(I6,5E13.5)') INUM,ZPRF,A,EE,TS1,TSUM                  
C                                                                       
C CALCULATE THE RELATIVE PROBABILITIES FOR ALL DECAY CHANNELS        
C                                                                       
      IF (DENSF.EQ.0.0) THEN                                            
       IF (DENSP.EQ.0.0) THEN                                           
         IF (DENSN.EQ.0.0) THEN                                         
      	   IF (DENSA.EQ.0.0) THEN                                        
C                                                                       
C NO REACTION IS POSSIBLE                                               
C                                                                       
      	     PROBF = 0.E0                                                
      	     PROBP = 0.E0                                                
      	     PROBN = 0.E0                                                
      	     PROBA = 0.E0                                                
      	     GO TO 70                                                    
      	   END IF                                                        
C                                                                       
C ALPHA EVAPORATION IS THE ONLY OPEN CHANNEL                            
C                                                                       
      	   RF = 0.0                                                      
      	   RP = 0.0                                                      
      	   RN = 0.0                                                      
      	   RA = 1.0                                                      
      	   GO TO 50                                                      
      	 END IF                                                          
C                                                                       
C ALPHA EMISSION AND NEUTRON EMISSION                                   
C                                                                       
      	 RF = 0.0                                                        
       	 RP = 0.0                                                       
       	 RN = 1.0                                                       
       	 RA = DENSA*2.0/DENSN*(AT/NT)**2                              
       	 GO TO 50                                                       
       END IF                                                           
C                                                                       
C ALPHA, PROTON AND NEUTRON EMISSION                                    
C                                                                       
       RF = 0.0                                                         
       RP = 1.0                                                         
       RN = DENSN/DENSP*(NT/PT)**2                                      
       RA = DENSA*2.0/DENSP*(AT/PT)**2                                
       GO TO 50                                                         
      END IF                                                            
C                                                                       
C HERE FISSION HAS TAKEN PLACE                                          
C                                                                       
      RF=1.D0                                                           
C                                                                       
C CRAMERS AND WEIDENMUELLER FACTORS FOR THE DYNAMICAL HINDRANCES OF     
C FISSION                                                               
C                                                                       
      IF (BET.LE.1.0D-16) THEN                                          
        CF=1.0D0                                                        
        WF=1.0D0                                                        
      ELSE IF (SBF.GT.0.0D0) THEN                                       
        CF=CRAM(BET,HOMEGA)                                             
C                                                                       
C     IF FISSION BARRIER EF=0.D0 THEN FISSION IS THE ONLY POSSIBLE      
C     CHANNEL. TO AVOID LOG(0) IN FUNCTION TAU                          
C     A.J. 7/28/93                                                      
C     
        IF (EF.LE.0.D0) THEN                                            
      	  RP = 0.D0                                                      
      	  RN = 0.D0                                                      
      	  RA = 0.D0                                                      
      	  GO TO 50                                                       
        ELSE                                                            
C                                                                       
C WARNING: SHELL CORR. REMOVED FROM BARRIER A.J. 7.7.96                 
C                                                                       
C	  TAUC=TAU(BET,HOMEGA,EF+ECGNZ(IN,IZ),FT)                             
C                                                                       
C TRANSIENT TIME TAU()                                                  
C                                                                       
	  TAUC=TAU(BET,HOMEGA,EF,FT)                                           
C                                                                       
        END IF       
        WFEX = (TAUC-TSUM)/TS1
        IF (WFEX.LT.0.D0) THEN
          WF = 1.D0
        ELSE
          WF = DEXP(-WFEX)
        ENDIF

c        IF (TAUC.GT.TSUM) THEN                                               
c           WFEX = (TAUC-TSUM)/ TS1 
c        ELSE                                 
c          WFEX = 1.D0        !!! KHS 5.1.2000
C                             The value 1 is wrong, corrected to 0.
c           WFEX = 0.D0
c        END IF                     
C      WRITE(6,*)'WFEX,A,ZPRF ',WFEX,A,ZPRF                              
c        IF (WFEX.LT.1.D-98) THEN                                        
c      	    WF = 1.D0                                                    
c        ELSEIF (WFEX.LT.300.D0) THEN                                    
c            WF = DEXP(-WFEX)                                            
c        ELSE                                                            
c            WF = DEXP(-300.D0)                                          
c        END IF                                                          
      ELSE                                                              
          CF=1.0D0                                                      
          WF=1.0D0                                                      
      END IF
                                                            
      IF (itest.eq.1) write(6,*) 'TSUM,WF,CF',TSUM,WF,CF

C     WRITE(6,*)'A,ZPRF,EE,WFEX,WF, ',A,ZPRF,EE,WFEX,WF                              

      TSUM = TSUM + TS1                                                                                                                     
C     
C     FOR EXCITATION ENERGIES BELOW 17 MEV GAMMAN/GAMMAF IS             
C     CALCULATED FROM AN  EMPIRICAL PARAMETRISATION                     
C     C.F. A.J. PROT. BUCH S. 168 OR DIPLOMA THESIS                     
C                                                                       
C     Z**2/A > 25 REQUIRED                                              
C                                                                       
C     IF ((EE.LE.17.D0).AND.(OPTLES.EQ.1).AND.                         
C    &   ((DFLOAT(IZ**2)/DFLOAT(AFP)).GT.25.D0)) THEN                  
C       GNGF = DEXP(DLOG(10.D0) * ( 28.8009D0 -                        
C    &	    0.793121D0*DFLOAT(IZ)**2/DFLOAT(AFP)))                      
C                                                                       
C  CHANGE BY G.K. AND A.H. 5.9.95                                       
C                                                                       
        TCONST=.7D0                                                     
	DCONST=12.D0/DSQRT(A)                                                  
	NPRF=A-ZPRF                                                            
        IF (OPTSHP.GE.2) THEN                                           
	   CALL PARITE (NPRF,PARC)                                             
 	   DCONST=DCONST*PARC                                                 
        ELSE                                                            
           DCONST= 0.D0                                                 
        END IF                                                          
C	IF (INTTYPE.EQ.1.AND.OPTLES.EQ.1) THEN                                
C          GNGF=A**.666666*TCONST/10.D0*DEXP((EF-SN+DCONST)/TCONST)     
C      IF ((EE.LE.17.D0).AND.(OPTLES.EQ.1).AND.                         
C     &   ((DFLOAT(IZ**2)/DFLOAT(AFP)).GT.25.D0)) THEN                  
      IF ((EE.LE.17.D0).AND.(OPTLES.EQ.1).AND.                          
     &    (IZ.GE.90).AND.(IN.GE.134)) THEN                              
C                                                                      
C CONSTANT CHANGED TO 5.0 ACCORD TO MORETTO & VANDENBOSCH A.J. 19.3.96  
C                                                                       
          GNGF=A**(2.D0/3.D0)*TCONST/10.D0*DEXP((EF-SN+DCONST)/TCONST)    
                                                                          
C                                                                       
C  IF THE EXCITATION ENERGY IS SO LOW THAT DENSN=0 ==> GN = 0           
C  FISSION REMAINS THE ONLY CHANNEL.                                    
C  A. J. 10.1.94                                                        
C                                                                       
        IF (GN.EQ.0.D0) THEN                                            
           RN = 0.D0                                                    
           RP = 0.D0                                                    
           RA = 0.D0                                                    
        ELSE                                                            
           RN=GNGF                                                      
           RP=GNGF*GP/GN                                                
           RA=GNGF*GA/GN                                                
        END IF                                                          
      ELSE  
       RN = GN/(GF*CF)                                                  
       RP = GP/(GF*CF)                                                  
       RA = GA/(GF*CF)                                                  
      END IF                                                            
 50   CONTINUE                                                          
C                                                                       
C RELATIVE DECAY PROBABILITIES                                          
C          

      DENOMI = RP+RN+RA+RF                                              

C Decay probabilities after transient time
      PROBF = RF/DENOMI                                                  
      PROBP = RP/DENOMI                                                 
      PROBN = RN/DENOMI                                                 
      PROBA = RA/DENOMI                                                 

C     New treatment of Grange-Weidenmueller Factor, 5.1.2000, KHS !!!

C Decay probabilites with transient time included
      PROBF = PROBF * WF
      PROBP = PROBP * (WF + (1.D0-WF)/(1.D0-PROBF)) 
      PROBN = PROBN * (WF + (1.D0-WF)/(1.D0-PROBF))
      PROBA = PROBA * (WF + (1.D0-WF)/(1.D0-PROBF))

  70  PTOTL = PROBP+PROBN+PROBA+PROBF   

      EE = EER                                                          
      ILAST = INUM                                                      
C      IF (INTTYPE.EQ.1) WRITE(1,*) IZ,IN,EE,PROBF                      
c      WRITE(1,*) IZ,IN,EE,PROBF
      RETURN                                                            
      END                                                               
C                                                                       
C********************************************************************** 
C                                                                       
      SUBROUTINE DENSNIV(A,Z,EE,ESOUS,DENS,BSHELL,BS,BK,TEMP,           
     &                   OPTSHP,OPTCOL,DEFBET)                          
C                                                                       
C     INPUT:                                                            
C             A,EE,ESOUS,OPTSHP,BS,BK,BSHELL,DEFBET                     
C                                                                       
C     LEVEL DENSITY PARAMETERS                                          
C     COMMON /ALD/    AV,AS,AK,OPTAFAN                                  
C     AV,AS,AK - VOLUME,SURFACE,CURVATURE DEPENDENCE OF THE             
C                LEVEL DENSITY PARAMETER                                
C     OPTAFAN - 0/1  AF/AN >=1 OR AF/AN ==1                             
C               RECOMMENDED IS OPTAFAN = 0                              
C---------------------------------------------------------------------  
C     OUTPUT: DENS,TEMP                                                 
C                                                                       
C   ____________________________________________________________________
C  /                                                                    
C  /  PROCEDURE FOR CALCULATING THE STATE DENSITY OF A COMPOUND NUCLEUS 
C  /____________________________________________________________________
C                                                                       
      INTEGER AFP,IZ,OPTSHP,OPTCOL,J,OPTAFAN                            
      REAL*8 A,EE,ESOUS,DENS,E,Y0,Y1,Y2,Y01,Y11,Y21,PA,BS,BK,TEMP       
C=====INSERTED BY KUDYAEV===============================================
      COMMON /ALD/ AV,AS,AK,OPTAFAN                                     
      REAL*8 ECR,ER,DELTAU,Z,DELTPP,PARA,PARZ,FE,HE,ECOR,ECOR1,Pi6      
      REAL*8 BSHELL,DELTA0,AV,AK,AS,PONNIV,PONFE,DEFBET,QR,SIG,FP       
C=======================================================================
C                                                                       
C                                                                       
C-----------------------------------------------------------------------
C     A                 MASS NUMBER OF THE DAUGHTER NUCLEUS             
C     EE                EXCITATION ENERGY OF THE MOTHER NUCLEUS         
C     ESOUS             SEPARATION ENERGY PLUS EFFECTIVE COULOMB BARRIER
C     DENS              STATE DENSITY OF DAUGHTER NUCLEUS AT EE-ESOUS-EC
C     BSHELL            SHELL CORRECTION                                
C     TEMP              NUCLEAR TEMPERATURE                             
C     E        LOCAL    EXCITATION ENERGY OF THE DAUGHTER NUCLEUS       
C     E1       LOCAL    HELP VARIABLE                                   
C     Y0,Y1,Y2,Y01,Y11,Y21                                              
C              LOCAL    HELP VARIABLES                                  
C     PA       LOCAL    STATE-DENSITY PARAMETER                         
C     EC                KINETIC ENERGY OF EMITTED PARTICLE WITHOUT      
C                        COULOMB REPULSION                              
C     IDEN              FAKTOR FOR SUBSTRACTING KINETIC ENERGY IDEN*TEMP
C     DELTA0            PAIRING GAP 12 FOR GROUND STATE                 
C                       14 FOR SADDLE POINT                             
C     EITERA            HELP VARIABLE FOR TEMPERATURE ITERATION         
C-----------------------------------------------------------------------
C                                                                       
C                                                                       
      PI6 = 3.1415926535D0**2 / 6.D0                                    
      ECR=10.0D0                                                        
      ER=28.0D0                                                         
      AFP=IDNINT(A)                                                     
      IZ=IDNINT(Z)                                                      
C                                                                       
C LEVEL DENSITY PARAMETER                                               
C                                                                       
      IF(OPTAFAN.EQ.1) THEN                                             
        PA = AV*A + AS*A**(2.D0/3.D0) + AK*A**(1.D0/3.D0)               
      ELSE                                                              
        PA = AV*A + AS*BS*A**(2.D0/3.D0) + AK*BK*A**(1.D0/3.D0)         
      END IF                                                            
C                                                                       
C PERPENDICULAR MOMENT OF INERTIA                                       
C                                                                       
C  FP = 2/5*M0*R0**2/HBAR**2 * A**(5/3) * (1 + DEFBET/3)                
C                                                                       
C  M0 NUCLEON MASS, R0 =1.2 FM                                          
C                                                                       
      FP = 0.01377937231D0 * A**(5.D0/3.D0) * (1.D0 + DEFBET/3.D0)      
C                                                                       
C PAIRING CORRECTIONS                                                   
C                                                                       
      IF (BS.GT.1.D0) THEN                                              
         DELTA0 = 14                                                    
      ELSE                                                              
         DELTA0 = 12                                                    
      END IF                                                            
      IF (ESOUS.GT.1.0D30) THEN                                         
        DENS = 0.D0                                                     
        TEMP=0.D0                                                       
        GO TO 100                                                       
      END IF                                                            
      E = EE-ESOUS                                                      
C     WRITE(6,*)'E ',E                                                  
      IF (E.LT.0.D0) THEN                                               
        DENS = 0.D0                                                     
        TEMP=0.D0                                                       
        GO TO 100                                                       
      END IF                                                            
C                                                                       
C SHELL CORRECTIONS                                                     
C                                                                       
      IF (OPTSHP.GT.0) THEN                                             
        DELTAU = BSHELL                                                 
        IF (OPTSHP.EQ.2) DELTAU = 0.D0                                  
        IF (OPTSHP.GE.2) THEN                                           
C                                                                       
C  pairing energy shift with condensation energy A.R.J. 10.03.97        
C                                                                       
          DELTPP = -0.25D0* (DELTA0/DSQRT(A))**2 * PA /PI6              
     &             + 2.D0*DELTA0/DSQRT(A)                               
C                                                                       
          CALL PARITE(A,PARA)                                           
          IF (PARA.LT.0.D0) THEN                                        
            E = E-DELTA0/DSQRT(A)                                       
          ELSE                                                          
            CALL PARITE(Z,PARZ)                                         
            IF (PARZ.GT.0.D0) THEN                                      
                E = E-2.D0*DELTA0/DSQRT(A)                              
            ELSE                                                        
                E = E                                                   
            END IF                                                      
          END IF                                                        
        ELSE                                                            
          DELTPP = 0.D0                                                 
        END IF                                                          
      ELSE                                                              
        DELTAU = 0.D0                                                   
        DELTPP = 0.D0                                                   
      END IF                                                            
      IF (E.LT.0.D0) THEN                                               
        E = 0.D0                                                        
        TEMP=0.D0                                                       
      END IF                                                            
C                                                                       
C     PONFE = -2.5D0*PA*E*A**(-4.D0/3.D0)                               
C  WASHING OUT IS MADE STRONGER ! G.K. 3.7.96                           
C                                                                       
      PONFE = -2.5D0*PA*E*A**(-4.D0/3.D0)                               
C                                                                       
      IF (PONFE.LT.-700.D0)  PONFE=-700.D0                              
      FE = 1.D0-DEXP(PONFE)                                             
      IF (E.LT.ECR) THEN                                                
C                                                                       
C     priv. comm. K.-H. SCHMIDT                                         
C                                                                       
         HE = 1.D0 - (1.D0 - E/ECR)**2                                  
C                                                                       
      ELSE                                                              
        HE = 1.D0                                                       
      END IF                                                            
C                                                                       
C EXCITATION ENERGY CORRECTED FOR PAIRING AND SHELL EFFECTS             
C WASHING OUT WITH EXCITATION ENERGY IS INCLUDED                        
C                                                                       
      ECOR = E+DELTAU*FE+DELTPP*HE                                      
C                                                                       
      IF (ECOR.LE.0.1D0) ECOR = 0.1D0                                   
C                                                                       
C statt 170.d0 a.r.j. 8.11.97                                           
C                                                                       
C                                                                       
C ITERATIVE PROCEDURE ACCORDING TO GROSSJEAN AND FELDMEIER              
C TO AVOID THE SINGULARITY E = 0                                        
C                                                                       
      IF (EE.LT.5.0D0) THEN                                             
        Y1 = DSQRT(PA*ECOR)                                             
        DO 20  J = 1,5,1                                                
          Y2 = PA*ECOR*(1.D0-DEXP(-Y1))                                 
          Y1 = DSQRT(Y2)                                                
 20     CONTINUE                                                        
        Y0 = PA/Y1                                                      
        TEMP=1.0D0/Y0                                                   
        DENS = DEXP(Y0*ECOR)/                                           
     &  ((ECOR**3*Y0)**0.5D0*(1.D0-0.5D0*Y0*ECOR*DEXP(-Y1))**0.5D0)*    
     &  DEXP(Y1)*(1.D0-DEXP(-Y1))*0.1477045D0                           
        IF (ECOR.LT.1.D0) THEN                                          
          ECOR1=1.0D0                                                   
          Y11 = DSQRT(PA*ECOR1)                                         
 25       CONTINUE                                                      
          DO 30  J = 1,7,1                                              
            Y21 = PA*ECOR1*(1.D0-DEXP(-Y11))                            
            Y11 = DSQRT(Y21)                                            
 30       CONTINUE                                                      
          Y01 = PA/Y11                                                  
          DENS = DENS*(Y01/Y0)**1.5D0                                   
          TEMP = TEMP*(Y01/Y0)**1.5D0                                   
        END IF                                                          
      ELSE                                                              
C                                                                       
C       WRITE(6,*)'PA,ECOR ',PA,ECOR                                    
C       ECOR = DMIN1(3000.D0,ECOR)                                      
C                                                                       
        PONNIV = 2.D0*DSQRT(PA*ECOR)                                    
        IF (PONNIV.GT.700.D0) PONNIV = 700.D0                           
C                                                                       
C FERMI GAS STATE DENSITY                                               
C                                                                       
        DENS = PA**(-0.25D0)*ECOR**(-1.25D0)*DEXP(PONNIV)               
     +  * 0.1477045D0                                                   
        TEMP = DSQRT(ECOR/PA)                                           
      END IF                                                            
 100  CONTINUE                                                          
C                                                                       
C SPIN CUTOFF PARAMETER                                                 
C                                                                       
      SIG = FP * TEMP                                                   
C                                                                       
C COLLECTIVE ENHANCEMENT                                                
C                                                                       
      IF (OPTCOL.EQ.1) THEN                                             
         CALL QROT(Z,A,DEFBET,SIG,ECOR,QR)                              
      ELSE                                                              
         QR   = 1.D0                                                    
      END IF                                                            
      DENS = DENS * QR                                                  
C                                                                       
C     WRITE(6,*)'AFP, IZ, ECOR, ECOR1',AFP,IZ,ECOR,ECOR1                
C                                                                       
      RETURN                                                            
      END                                                               
C                                                                       
C******************************************************************     
C                                                                       
      SUBROUTINE QROT(Z,A,BET,SIG,U,QR)                                 
C                                                                       
C QROT INCLUDING DAMPING                                                
C                                                                       
C INPUT: Z,A,BET,SIG,U                                                  
C                                                                       
C OUTPUT: QR - COLLECTIVE ENHANCEMENT FACTOR                            
C                                                                       
C SEE  JUNGHANS ET AL., NUCL. PHYS. A 629 (1998) 635                    
C                                                                       
C                                                                       
C   FR(U) EXPONENTIAL FUNCTION TO DEFINE DAMPING                        
C   UCR   CRITICAL ENERGY FOR DAMPING                                   
C   DCR   WIDTH OF DAMPING                                              
C   BET   BETA-DEFORMATION !                                            
C   SIG   PERPENDICULAR SPIN CUTOFF FACTOR                              
C     U   ENERGY                                                        
C    QR   COEFFICIENT OF COLLECTIVE ENHANCEMENT                         
C     A   MASS NUMBER                                                   
C     Z   CHARGE NUMBER                                                 
C                                                                       
      REAL*8 UCR,DCR,BET,SIG,U,QR,PONQ,A,Z,DN,N,DZ                      
C      FR(U) = 1.D0/(1.D0+DEXP((U-UCR)/DCR))                            
C      UCR = 120.D0 *BET*BET*A**(1.D0/3.D0)                             
C      IF (BET.GE.0.001D0) THEN                                         
C        DCR = 1400.D0*BET*BET/A**(2.D0/3.D0)                           
C      ELSE                                                             
C         DCR = 1.D-16                                                  
C      END IF                                                           
C                                                                       
       DCR   = 10.D0                                                    
C                                                                       
C      IF (DABS(BET).GT.(0.50D0)) THEN                                  
         UCR = 40.D0                                                    
C      ELSE                                                             
C         UCR = 40.D0                                                   
C      END IF                                                           
C                                                                       
      IF (DABS(BET)-0.15D0) 10,10,11                                    
C                                                                       
C 10   QR = 1.D0                                                        
C      GO TO 12                                                         
  10  N = A - Z                                                         
      DZ = DABS(Z-82.D0)                                                
      IF (N.GT.104) THEN                                                
       DN = DABS(N-126.D0)                                              
      ELSE                                                              
       DN = DABS(N-82.D0)                                               
      END IF                                                            
C      BET = 0.02D0 + 0.004D0 *(DN + DABS(Z-82.D0))                     
C      BET = 0.01D0 + 0.003D0*DN + 0.005D0*DZ                           
      BET = 0.022D0 + 0.003D0*DN + 0.005D0*DZ                           
C                                                                       
      SIG = 25.D0*BET**2 * SIG                                          
C                                                                       
C  NO VIBRATIONAL ENHANCEMENT                                           
C                                                                       
C      SIG = 1.D0                                                       
C                                                                       
 11   PONQ = (U-UCR)/DCR                                                
      IF (PONQ.GT.700.D0) PONQ = 700.D0                                 
      IF (SIG.LT.1.D0) SIG = 1.D0                                       
      QR = 1.D0/(1.D0+DEXP(PONQ)) * (SIG-1.D0) + 1.D0                   
C      QR = 1.D0/(1.D0+DEXP(PONQ)) * SIG                                
      IF (QR.LT.1.D0) QR = 1.D0                                         
 12   RETURN                                                            
      END                                                               
C                                                                       
C********************************************************************** 
C                                                                       
      SUBROUTINE MGLW(A,Z,REFOPT,EL)                                    
C                                                                       
C      _______________________________________________________          
C     /                                                      /          
C     /  MODEL DE LA GOUTTE LIQUIDE DE C. F. WEIZSACKER.     /          
C     /  USUALLY AN OBSOLETE OPTION                          /          
C     /______________________________________________________/          
C                                                                       
      INTEGER A1,Z1,REFOPT                                              
      REAL*8 A,Z,XV,XS,XC,XA,EL,ECNZ                                    
C     REAL*8 ECNZ                                                       
C     DIMENSION ECNZ(0:153,0:98)                                        
      COMMON /EC2SUB/ ECNZ(0:153,0:98)                                  
C                                                                       
      A1 = IDNINT(A)                                                    
      Z1 = IDNINT(Z)                                                    
      IF ((A.LE.0.01).OR.(Z.LT.0.01)) THEN                              
        EL = 1.0D38                                                    
        GO TO 50                                                        
      END IF                                                            
      XV = -15.56*A                                                     
      XS = 17.23*A**(2.0/3.0)                                           
      IF (A.GT.1.0) THEN                                                
        XC = 0.7*Z*(Z-1.0)*(A-1.0)**(-1.D0/3.D0)                        
        ELSE                                                            
          XC = 0.0                                                      
      END IF                                                            
      XA = 23.6*((A-2.0*Z)**2/A)                                        
      EL = XV+XS+XC+XA                                                  
 50   CONTINUE                                                          
C                                                                       
      RETURN                                                            
      END                                                               
C                                                                       
C********************************************************************** 
C                                                                       
C ATTENTION: HERE FINITE RANGE LIQUID DROP MODEL IS USED !!!            
C                                                                       
C                                                                       
      SUBROUTINE MGLMS(A,Z,REFOPT4,EL)                                  
C                                                                       
C     ON INPUT:  A,Z,REFOPT4                                            
C     ON OUTPUT: EL                                                     
C   ____________________________________________________________________
C  /  USING FUNCTION EFLMAC(IA,IZ,0)                                    
C  /                                                                    
C  /  REFOPT4 = 0 : WITHOUT MICROSCOPIC CORRECTIONS                     
C  /  REFOPT4 = 1 : WITH SHELL CORRECTION                               
C  /  REFOPT4 = 2 : WITH PAIRING CORRECTION                             
C  /  REFOPT4 = 3 : WITH SHELL- AND PAIRING CORRECTION                  
C  /____________________________________________________________________
C                                                                       
      INTEGER A1,Z1,REFOPT4                                             
      REAL*8 A,Z,EL,ECNZ,EFLMAC                                         
C     DIMENSION ECNZ(0:153,0:98)                                        
      COMMON /EC2SUB/ ECNZ(0:153,0:98)                                  
C                                                                       
C-----------------------------------------------------------------------
C     A1       LOCAL    MASS NUMBER (INTEGER VARIABLE OF A)             
C     Z1       LOCAL    NUCLEAR CHARGE (INTEGER VARIABLE OF Z)          
C     REFOPT4           OPTION, SPECIFYING THE MASS FORMULA (SEE ABOVE) 
C     A                 MASS NUMBER                                     
C     Z                 NUCLEAR CHARGE                                  
C     DEL               PAIRING CORRECTION                              
C     EL                BINDING ENERGY                                  
C     ECNZ( , )         TABLE OF SHELL CORRECTIONS                      
C-----------------------------------------------------------------------
C                                                                       
      A1 = IDNINT(A)                                                    
      Z1 = IDNINT(Z)                                                    
C                                                                       
C      IF (((A.LE.0.0).OR.(Z.LE.0.0)).OR.(A-Z.LE.0.0)) THEN             
C                                                                       
      IF ( (A1.LE.0).OR.(Z1.LE.0).OR.((A1-Z1).LE.0) )  THEN 
C Modif pour récupérer une masse p et n correcte:
         EL=0.                  
C        EL = 1.0E38                                                     
        GO TO 50                                                        
      END IF                                                            
C      WRITE(6,*)'A1,Z1 ',A1,Z1                                         
C                                                                       
C BINDING ENERGY INCL. PAIRING CONTR. IS CALCULATED FROM                
C FUNCTION EFLMAC                                                       
C                                                                       
      EL = EFLMAC(A1,Z1,0,REFOPT4)                                      
C                                                                       
      IF (REFOPT4.GT.0) THEN                                            
        IF (REFOPT4.NE.2) EL = EL+ECNZ(A1-Z1,Z1)                        
C                                                                       
C PAIRING IS DONE INSIDE OF FUNCTION EFLMAC A.J. 17.6.96                
C                                                                       
C        IF (REFOPT4.GT.1) THEN                                         
C          CALL APPARIEM(A,Z,DEL)                                       
C          EL = EL+DEL                                                  
C        END IF                                                         
      END IF                                                            
 50   CONTINUE                                                          
C                                                                       
      RETURN                                                            
      END                                                               

C                                                                       
C CHANGED TO CALCULATE TOTAL BINDING ENERGY INSTEAD OF MASS EXCESS.     
C SWITCH FOR PAIRING INCLUDED AS WELL.                                  
C BINDING = EFLMAC(IA,IZ,0,OPTSHP)                                      
C FORTRAN TRANSCRIPT OF /U/GREWE/LANG/EEX/FRLDM.C                       
C A.J. 15.07.96                                                         
C                                                                       
                                                                        
C/**********************************************************************
C *                                                                     
C * this function will calculate the liquid-drop nuclear mass for spheri
C * configuration according to the preprint NUCLEAR GROUND-STATE        
C * MASSES and DEFORMATIONS by P. M"oller et al. from August 16, 1993 p.
C * All constants are taken from this publication for consistency.      
C *                                                                     
C * Parameters:                                                         
C *   a:    nuclear mass number                                         
C *   z:    nuclear charge                                              
C *   flag:     0       - return mass excess                            
C *         otherwise   - return pairing (= -1/2 dpn + 1/2 (Dp + Dn))   
C **********************************************************************
C                                                                       
      DOUBLE PRECISION FUNCTION EFLMAC(IA,IZ,FLAG,OPTSHP)               
C                                                                       
      INTEGER*4 IA,IZ,FLAG,IN,OPTSHP                                    
      REAL*8 Z,N,A,AV,AS,A0,C1,C4,B1,B3,F,CA,W,DP,DN,DPN,EFL,PI         
      REAL*8 RMAC,BS,H,R0,KF,KS,KV,RP,AY,ADEN,X0,Y0,MH,MN,ESQ,AEL,I     
      PI = 3.141592653589793238D0                                       
C                                                                       
C double eflmac(int a, int z, int flag) {                               
C  double Z,N,A,aV,aS,a0,c1,c4,B1,B3,f,ca,W,Dp,Dn,dpn,Efl;              
C  double rmac,BS,h,r0,kF,kS,kV,rp,aY,aden,x0,y0,MH,Mn,esq,ael,I;       
C  int n;                                                               
                                                                        
C /* fundamental constants */                                           
C /* hydrogen-atom mass excess */                                       
      MH  = 7.289034D0                                                  
C /* neutron mass excess       */                                       
      MN  = 8.071431D0                                                  
C /* electronic charge squared */                                       
      ESQ = 1.4399764D0                                                 
C                                                                       
C /* constants from considerations other than nucl. masses */           
C /* electronic binding        */                                       
      AEL = 1.433D-5                                                    
C /* proton rms radius         */                                       
      RP  = 0.8D0                                                       
C /* nuclear radius constant   */                                       
      R0  = 1.16D0                                                      
C /* range of Yukawa-plus-expon. potential */                           
      AY  = 0.68D0                                                      
C /* range of Yukawa function used to generate                          
C     nuclear charge distribution */                                    
      ADEN= 0.70D0                                                      
C /* constants from considering odd-even mass differences */            
C /* average pairing gap       */                                       
      RMAC= 4.80D0                                                      
C /* neutron-proton interaction*/                                       
      H   = 6.6D0                                                       
C /* Wigner constant           */                                       
      W   = 30.D0                                                       
C                                                                       
C /* adjusted parameters */                                             
C /* Volume energy              */                                      
      AV  = 16.00126D0                                                  
C /* volume asymmetry           */                                      
      KV  =  1.92240D0                                                  
C /* Surface energy             */                                      
      AS  = 21.18466D0                                                  
C /* surface asymmetry          */                                      
      KS  =  2.345D0                                                    
C /* A^0 constant               */                                      
      A0  =  2.615D0                                                    
C /* charge asymmetry           */                                      
      CA  =  0.10289D0                                                  
C                                                                       
C /* we will account for deformation by using the microscopic           
C   corrections tabulated from p. 68ff */                               
      BS=1.D0                                                           
C                                                                       
      Z   = DFLOAT(IZ)                                                  
      A   = DFLOAT(IA)                                                  
      IN  = IA-IZ                                                       
      N   = DFLOAT(IN)                                                  
      DN  = RMAC*BS/N**(1.D0/3.D0)                                      
      DP  = RMAC*BS/Z**(1.D0/3.D0)                                      
      DPN = H/BS/A**(2.D0/3.D0)                                         
C                                                                       
      C1  = 3.D0/5.D0*ESQ/R0                                            
      C4  = 5.D0/4.D0*(3.D0/(2.D0*PI))**(2.D0/3.D0) * C1                
C                                                                       
C                                                                       
      KF  = (9.D0*PI*Z/(4.D0*A))**(1.D0/3.D0)/R0                        
      F   = -1.D0/8.D0*RP*RP*ESQ/R0**3 *                                
     &      (145.D0/48.D0 - 327.D0/2880.D0*KF**2 * RP**2 +              
     &      1527.D0/1209600.D0*KF**4 * RP**4)                           
C                                                                       
C      WRITE(6,*)'C1 ',C1                                               
C      WRITE(6,*)'C4 ',C4                                               
C      WRITE(6,*)'KF ',KF                                               
C      WRITE(6,*)'F  ',F                                                
C                                                                       
      I   = (N-Z)/A                                                     
C                                                                       
      X0  = R0 * A**(1.D0/3.D0) / AY                                    
      Y0  = R0 * A**(1.D0/3.D0) / ADEN                                  
C                                                                       
      B1  = 1.D0 - 3.D0/(X0**2) +                                       
     &      (1.D0 + X0) * (2.D0 + 3.D0/X0 + 3.D0/X0**2)                 
     &       * DEXP(-2.D0*X0)                                           
C                                                                       
      B3  = 1.D0 - 5.D0/Y0**2 * (1.D0 - 15.D0/(8.D0*Y0) +               
     &      21.D0/(8.D0 * y0**3) -                                      
     &      3.D0/4.D0 * (1.D0 + 9.D0/(2.D0*Y0) + 7.D0/Y0**2 +           
     &      7.D0/(2.D0 * Y0**3)) * DEXP(-2.D0*Y0))                      
C                                                                       
C NOW CALULATION OF TOTAL BINDING ENERGY A.J. 16.7.96                   
C                                                                       
C      EFL = MH*Z + MN*N - AV*(1.D0-KV*I*I)*A                           
      EFL = - AV*(1.D0-KV*I*I)*A                                        
     &	    + AS*(1.D0-KS*I*I)*B1 * A**(2.D0/3.D0) + A0                  
     &	    + C1*Z*Z*B3/A**(1.D0/3.D0)                                   
     &	    - C4*Z**(4.D0/3.D0)/A**(1.D0/3.D0)                           
     &	    + F*Z**2/A -CA*(N-Z) - AEL * Z**(2.39D0)                     
C                                                                       
      IF ((IN.EQ.IZ).AND.(MOD(IN,2).EQ.1).AND.(MOD(IZ,2).EQ.1)) THEN                  
C #ifdef DEBUG                                                          
C  printf("z and n are odd and equal !\n");                             
C        WRITE(6,*)'z and n are odd and equal !'                        
C #endif                                                                
C       /* N and Z odd and equal */                                     
        EFL = EFL + W*(DABS(I)+1.D0/A)                                  
      ELSE                                                              
        EFL= EFL + W* DABS(I)                                           
      END IF                                                            
C                                                                       
C PAIRING IS MADE OPTIONAL                                              
C                                                                       
      IF (OPTSHP.GE.2) THEN                                             
C                                                                       
C /* average pairing */                                                 
        IF ((MOD(IN,2).EQ.1).AND.(MOD(IZ,2).EQ.1)) THEN                               
C #ifdef DEBUG                                                          
C     printf("z and n are odd !\n");                                    
C        WRITE(6,*)'z and n are odd !'                                  
C #endif                                                                
          EFL = EFL - DPN                                               
        END IF                                                          
        IF (MOD(IN,2).EQ.1) THEN                                             
          EFL = EFL + DN                                                
C#ifdef DEBUG                                                           
C    printf("n is odd !\n");                                            
C         WRITE(6,*)'n is odd !'                                        
C#endif                                                                 
        END IF                                                          
C                                                                       
        IF (MOD(IZ,2).EQ.1) THEN                                             
C #ifdef DEBUG                                                          
C     printf("z is odd !\n");                                           
C           WRITE(6,*)'z is odd !'                                      
C #endif                                                                
           EFL    = EFL + DP                                            
        END IF                                                          
C                                                                       
C END IF FOR PAIRING TERM                                               
C                                                                       
      END IF                                                            
C                                                                       
      IF (FLAG.NE.0) THEN                                                    
         EFLMAC =  (0.5D0*(DN + DP)-0.5D0*DPN)                          
      ELSE                                                              
         EFLMAC = EFL                                                   
      END IF                                                            
      RETURN                                                            
      END                                                               
C                                                                       
C*******************************************************************    
C                                                                       
      SUBROUTINE APPARIEM(A,Z,DEL)                                      
C                                                                       
C     ON INPUT:  A,Z                                                    
C     ON OUTPUT: DEL                                                    
C   ____________________________________________________________________
C  /                                                                    
C  /  CALCUL DE LA CORRECTION, DUE A L'APPARIEMENT, DE L'ENERGIE DE     
C  /  LIAISON D'UN NOYAU                                                
C  /                                                                    
C  /  PROCEDURE FOR CALCULATING THE PAIRING CORRECTION TO THE BINDING   
C  /  ENERGY OF A SPECIFIC NUCLEUS                                      
C  /____________________________________________________________________
C                                                                       
      REAL*8 A,Z,PARA,PARZ,DEL                                          
C                                                                       
C-----------------------------------------------------------------------
C     A                 MASS NUMBER                                     
C     Z                 NUCLEAR CHARGE                                  
C     PARA              HELP VARIABLE FOR PARITY OF A                   
C     PARZ              HELP VARIABLE FOR PARITY OF Z                   
C     DEL               PAIRING CORRECTION                              
C-----------------------------------------------------------------------
C                                                                       
      CALL PARITE(A,PARA)                                               
      IF (PARA.LT.0.D0) THEN                                            
        DEL = 0.D0                                                      
        ELSE                                                            
          CALL PARITE(Z,PARZ)                                           
          IF (PARZ.GT.0.D0) THEN                                        
            DEL = -12.D0/DSQRT(A)                                       
            ELSE                                                        
              DEL = 12.D0/DSQRT(A)                                      
          END IF                                                        
      END IF                                                            
C                                                                       
      RETURN                                                            
      END                                                               
C                                                                       
C*******************************************************************    
C                                                                       
      SUBROUTINE PARITE(N,PAR)                                          
C                                                                       
C     ON INPUT:  N                                                      
C     ON OUTPUT: PAR                                                    
C   ____________________________________________________________________
C  /                                                                    
C  /  CALCUL DE LA PARITE DU NOMBRE N                                   
C  /                                                                    
C  /  PROCEDURE FOR CALCULATING THE PARITY OF THE NUMBER N.             
C  /  RETURNS -1 IF N IS ODD AND +1 IF N IS EVEN                        
C  /____________________________________________________________________
C                                                                       
      REAL*8 N,PAR,N1,N2,N3                                             
C                                                                       
C-----------------------------------------------------------------------
C     N                 NUMBER TO BE TESTED                             
C     N1,N2             HELP VARIABLES                                  
C     PAR               HELP VARIABLE FOR PARITY OF N                   
C-----------------------------------------------------------------------
C                                                                       
      N3 = DFLOAT(IDNINT(N))                                            
      N1 = N3/2.0                                                       
      N2 = N1-DINT(N1)                                                  
      IF (N2.GT.0.D0) THEN                                              
        PAR = -1.0                                                      
        ELSE                                                            
          PAR = 1.0                                                     
      END IF                                                            
C                                                                       
      RETURN                                                            
      END                                                               
C                                                                       
C******************************************************************     
C                                                                       
      DOUBLE PRECISION FUNCTION TAU(BET,HOMEGA,EF,T)                    
C                                                                       
C   INPUT : BET, HOMEGA, EF, T                                          
C   OUTPUT: TAU - RISE TIME IN WHICH THE FISSION WIDTH HAS REACHED      
C           90 PERCENT OF ITS FINAL VALUE                               
C                                                                       
C   BETA   - NUCLEAR VISCOSITY                                          
C   HOMEGA - CURVATURE OF POTENTIAL                                     
C   EF     - FISSION BARRIER                                            
C   T      - NUCLEAR TEMPERATURE                                        
C                                                                       
      REAL*8 BET,HOMEGA,EF,T,TLIM,tau1,tau2                                       
      TLIM = 8.D0 * EF                                                  
      IF (T.GT.TLIM) T = TLIM   

C      modified BJ and KHS 6.1.2000:
C     IF (BET/(20.D0*(HOMEGA/6.582122D0)).LE.1.0D0) THEN                
      IF (BET/(SQRT(2.D0)*10.D0*(HOMEGA/6.582122D0)).LE.1.0D0) THEN                
       TAU = DLOG(10.0D0*EF/T)/(BET*1.0D+21)                          
      ELSE                                                              
       TAU = DLOG(10.0D0*EF/T)/                                       
     & (2.D0*(10.D0*HOMEGA/6.582122D0)**2)*(BET*1.0D-21)
      END IF                                                            
      RETURN                                                            
      END                                                               
C                                                                       
C******************************************************************     
C                                                                       
      DOUBLE PRECISION FUNCTION CRAM(BET,HOMEGA)                        
C                                                                       
C  INPUT : BET, HOMEGA  NUCLEAR VISCOSITY + CURVATURE OF POTENTIAL      
C  OUTPUT: KRAMERS FAKTOR  - REDUCTION OF THE FISSION PROBABILITY       
C          INDEPENDENT OF EXCITATION ENERGY                             
C                                                                       
      REAL*8 BET,HOMEGA,REL                                             
      REL= BET/(20.D0*HOMEGA/6.582122D0)                                
      CRAM=DSQRT(1.0D0+REL**2)-REL 
C     limitation introduced   6.1.2000  by  KHS
      if (cram.gt.1.D0) cram = 1.D0                                     
      RETURN                                                            
      END                                                               
C                                                                       
C********************************************************************** 
C                                                                       
      DOUBLE PRECISION FUNCTION BIPOL(IFLAG,Y)                          
C                                                                       
C     CALCULATION OF THE SURFACE BS OR CURVATURE BK OF A NUCLEUS        
C     RELATIVE TO THE SPHERICAL CONFIGURATION                           
C     BASED ON  MYERS, DROPLET MODEL FOR ARBITRARY SHAPES               
C                                                                       
C     INPUT: IFLAG - 0/1 BK/BS CALCULATION                              
C             Y    - (1 - X) COMPLEMENT OF THE FISSILITY                
C                                                                       
C     LINEAR INTERPOLATION OF BS BK TABLE                               
C                                                                       
      INTEGER IFLAG,I                                                   
      REAL*8 Y,BS(1:52),BK(1:52)                                        
      DATA                                                              
     & BK/1.00000,1.00087,1.00352,1.00799,1.01433,1.02265,1.03306,      
     & 1.04576,1.06099,1.07910,1.10056,1.12603,1.15651,1.19348,         
     & 1.23915,1.29590,1.35951,1.41013,1.44103,1.46026,1.47339,         
     & 1.48308,1.49068,1.49692,1.50226,1.50694,1.51114,1.51502,         
     & 1.51864,1.52208,1.52539,1.52861,1.53177,1.53490,1.53803,         
     & 1.54117,1.54473,1.54762,1.55096,1.55440,1.55798,1.56173,         
     & 1.56567,1.56980,1.57413,1.57860,1.58301,1.58688,1.58688,         
     & 1.58688,1.58740,1.58740/                                         
      DATA                                                              
     & BS/1.00000,1.00086,1.00338,1.00750,1.01319,                      
     &    1.02044,1.02927,1.03974,                                      
     & 1.05195,1.06604,1.08224,1.10085,1.12229,1.14717,1.17623,1.20963, 
     & 1.24296,1.26532,1.27619,1.28126,1.28362,1.28458,1.28477,1.28450, 
     & 1.28394,1.28320,1.28235,1.28141,1.28042,1.27941,1.27837,1.27732, 
     & 1.27627,1.27522,1.27418,1.27314,1.27210,1.27108,1.27006,1.26906, 
     & 1.26806,1.26707,1.26610,1.26514,1.26418,1.26325,1.26233,1.26147, 
     & 1.26147,1.26147,1.25992,1.25992/                                 
      I = IDINT(Y/2.D-02) + 1                                           
      IF (IFLAG.EQ.1) THEN                                              
        BIPOL = BS(I) + (BS(I+1)-BS(I))/2.D-02 * (Y - 2.D-02*(I-1))     
      ELSE                                                              
        BIPOL = BK(I) + (BK(I+1)-BK(I))/2.D-02 * (Y - 2.D-02*(I-1))     
      END IF                                                            
      RETURN                                                            
      END                                                               

C                                                                       
C******************************************************************     
C                                                                       
      SUBROUTINE BARFIT(IZ,IA,IL,SBFIS,SEGS,SELMAX)                     
C     VERSION FOR 32BIT COMPUTER                                        
C     THIS SUBROUTINE RETURNS THE BARRIER HEIGHT BFIS, THE              
C     GROUND-STATE ENERGY SEGS, IN MEV, AND THE ANGULAR MOMENTUM        
C     AT WHICH THE FISSION BARRIER DISAPPEARS, LMAX, IN UNITS OF        
C     H-BAR, WHEN CALLED WITH INTEGER AGUMENTS IZ, THE ATOMIC           
C     NUMBER, IA, THE ATOMIC MASS NUMBER, AND IL, THE ANGULAR           
C     MOMENTUM IN UNITS OF H-BAR. (PLANCK'S CONSTANT DIVIDED BY         
C     2*PI).                                                            
C                                                                       
C        THE FISSION BARRIER FO IL = 0 IS CALCULATED FROM A 7TH         
C     ORDER FIT IN TWO VARIABLES TO 638 CALCULATED FISSION              
C     BARRIERS FOR Z VALUES FROM 20 TO 110. THESE 638 BARRIERS ARE      
C     FIT WITH AN RMS DEVIATION OF 0.10 MEV BY THIS 49-PARAMETER        
C     FUNCTION.                                                         
C     IF BARFIT IS CALLED WITH (IZ,IA) VALUES OUTSIDE THE RANGE OF      
C     THE BARRIER HEIGHT IS SET TO 0.0, AND A MESSAGE IS PRINTED        
C     ON THE DEFAULT OUTPUT FILE.                                       
C                                                                       
C        FOR IL VALUES NOT EQUAL TO ZERO, THE VALUES OF L AT WHICH      
C     THE BARRIER IS 80% AND 20% OF THE L=0 VALUE ARE RESPECTIVELY      
C     FIT TO 20-PARAMETER FUNCTIONS OF Z AND A, OVER A MORE             
C     RESTRICTED RANGE OF A VALUES, THAN IS THE CASE FOR L = 0.         
C     THE VALUE OF L WHERE THE BARRIER DISAPPEARS, LMAX IS FIT TO       
C     A 24-PARAMETER FUNCTION OF Z AND A, WITH THE SAME RANGE OF        
C     Z AND A VALUES AS L-80 AND L-20.                                  
C        ONCE AGAIN, IF AN (IZ,IA) PAIR IS OUTSIDE OF THE RANGE OF      
C     VALIDITY OF THE FIT, THE BARRIER VALUE IS SET TO 0.0 AND A        
C     MESSAGE IS PRINTED. THESE THREE VALUES (BFIS(L=0),L-80, AND       
C     L-20) AND THE CONSTRINTS OF BFIS = 0 AND D(BFIS)/DL = 0 AT        
C     L = LMAX AND L=0 LEAD TO A FIFTH-ORDER FIT TO BFIS(L) FOR         
C     L>L-20. THE FIRST THREE CONSTRAINTS LEAD TO A THIRD-ORDER FIT     
C     FOR THE REGION L < L-20.                                          
C                                                                       
C        THE GROUND STATE ENERGIES ARE CALCULATED FROM A                
C     120-PARAMETER FIT IN Z, A, AND L TO 214 GROUND-STATE ENERGIES     
C     FOR 36 DIFFERENT Z AND A VALUES.                                  
C     (THE RANGE OF Z AND A IS THE SAME AS FOR L-80, L-20, AND          
C     L-MAX)                                                            
C                                                                       
C        THE CALCULATED BARRIERS FROM WHICH THE FITS WERE MADE WERE     
C     CALCULATED IN 1983-1984 BY A. J. SIERK OF LOS ALAMOS              
C     NATIONAL LABORATORY GROUP T-9, USING YUKAWA-PLUS-EXPONENTIAL      
C     DOUBLE FOLDED NUCLEAR ENERGY, EXACT COULOMB DIFFUSENESS           
C     CORRECTIONS, AND DIFFUSE-MATTER MOMENTS OF INERTIA.               
C     THE PARAMETERS OF THE MODEL R-0 = 1.16 FM, AS 21.13 MEV,          
C     KAPPA-S = 2.3, A = 0.68 FM.                                       
C     THE DIFFUSENESS OF THE MATTER AND CHARGE DISTRIBUTIONS USED       
C     CORRESPONDS TO A SURFACE DIFFUSENESS PARAMETER (DEFINED BY        
C     MYERS) OF 0.99 FM. THE CALCULATED BARRIERS FOR L = 0 ARE          
C     ACCURATE TO A LITTLE LESS THAN 0.1 MEV; THE OUTPUT FROM           
C     THIS SUBROUTINE IS A LITTLE LESS ACCURATE. WORST ERRORS MAY BE    
C     AS LARGE AS 0.5 MEV; CHARACTERISTIC UNCERTAINY IS IN THE RANGE    
C     OF 0.1-0.2 MEV. THE RMS DEVIATION OF THE GROUND-STATE FIT         
C     FROM THE 214 INPUT VALUES IS 0.20 MEV. THE MAXIMUM ERROR          
C     OCCURS FOR LIGHT NUCLEI IN THE REGION WHERE THE GROUND STATE      
C     IS PROLATE, AND MAY BE GREATER THAN 1.0 MEV FOR VERY NEUTRON      
C     DEFICIENT NUCLEI, WITH L NEAR LMAX. FOR MOST NUCLEI LIKELY TO     
C     BE ENCOUNTERED IN REAL EXPERIMENTS, THE MAXIMUM ERROR IS          
C     CLOSER TO 0.5 MEV, AGAIN FOR LIGHT NUCLEI AND L NEAR LMAX.        
C                                                                       
C     WRITTEN BY A. J. SIERK, LANL T-9                                  
C     VERSION 1.0 FEBRUARY, 1984                                        
C                                                                       
C     THE FOLLOWING IS NECESSARY FOR 32-BIT MACHINES LIKE DEC VAX,      
C     IBM, ETC                                                          
C                                                                       
C                                                                       
C      IMPLICIT REAL*8 (A-H,O-Z)                                        
      REAL*4 SBFIS,SEGS,SELMAX                                          
      REAL*8 ELZCOF,ELMCOF,EMNCOF,PA,PZ,PL,EMXCOF,EGSCOF,EGS1,EGS2      
      REAL*8 EGS3,EGS4,A,Z,AMIN,AMAX,AMIN2,AMAX2,AA,ZZ,BFIS           
      REAL*8 BFIS0,ELL,EL,EGS,EL80,EL20,ELMAX,SEL80,SEL20,X,Y,Q,QA,QB   
      REAL*8 AJ,AK,A1,A2                                                
      INTEGER*4 IA,IZ,IL,I,J,K,L,M                                        
      DIMENSION  ELZCOF(7,7),ELMCOF(5,4),EMNCOF(5,4),PA(7),PZ(7),PL(10) 
      DIMENSION  EMXCOF(6,4),EGSCOF(5,6,4),EGS1(5,6),EGS2(5,6),         
     1           EGS3(5,6),EGS4(5,6)                                    
      EQUIVALENCE (EGS1(1,1),EGSCOF(1,1,1)), (EGS2(1,1),EGSCOF(1,1,2)), 
     1            (EGS3(1,1),EGSCOF(1,1,3)), (EGS4(1,1),EGSCOF(1,1,4))  
      DATA EMNCOF                                                       
     1/-9.01100D+2,-1.40818D+3, 2.77000D+3,-7.06695D+2, 8.89867D+2,     
     2  1.35355D+4,-2.03847D+4, 1.09384D+4,-4.86297D+3,-6.18603D+2,     
     3 -3.26367D+3, 1.62447D+3, 1.36856D+3, 1.31731D+3, 1.53372D+2,     
     4  7.48863D+3,-1.21581D+4, 5.50281D+3,-1.33630D+3, 5.05367D-2/     
      DATA ELMCOF                                                       
     1 /1.84542D+3,-5.64002D+3, 5.66730D+3,-3.15150D+3, 9.54160D+2,     
     2 -2.24577D+3, 8.56133D+3,-9.67348D+3, 5.81744D+3,-1.86997D+3,     
     3  2.79772D+3,-8.73073D+3, 9.19706D+3,-4.91900D+3, 1.37283D+3,     
     4 -3.01866D+1, 1.41161D+3,-2.85919D+3, 2.13016D+3,-6.49072D+2/     
      DATA EMXCOF /                                                     
     19.43596D4,-2.241997D5,2.223237D5,-1.324408D5,4.68922D4,-8.83568D3,
     2-1.655827D5,4.062365D5,-4.236128D5,2.66837D5,-9.93242D4,1.90644D4,
     3 1.705447D5,-4.032D5,3.970312D5,-2.313704D5,7.81147D4,-1.322775D4,
     4-9.274555D4,2.278093D5,-2.422225D5,1.55431D5,-5.78742D4,9.97505D3/
      DATA ELZCOF                                                       
     1 /5.11819909D+5,-1.30303186D+6, 1.90119870D+6,-1.20628242D+6,     
     2  5.68208488D+5, 5.48346483D+4,-2.45883052D+4,                    
     3 -1.13269453D+6, 2.97764590D+6,-4.54326326D+6, 3.00464870D+6,     
     4 -1.44989274D+6,-1.02026610D+5, 6.27959815D+4,                    
     5  1.37543304D+6,-3.65808988D+6, 5.47798999D+6,-3.78109283D+6,     
     6  1.84131765D+6, 1.53669695D+4,-6.96817834D+4,                    
     7 -8.56559835D+5, 2.48872266D+6,-4.07349128D+6, 3.12835899D+6,     
     8 -1.62394090D+6, 1.19797378D+5, 4.25737058D+4,                    
     9  3.28723311D+5,-1.09892175D+6, 2.03997269D+6,-1.77185718D+6,     
     A  9.96051545D+5,-1.53305699D+5,-1.12982954D+4,                    
     B  4.15850238D+4, 7.29653408D+4,-4.93776346D+5, 6.01254680D+5,     
     C -4.01308292D+5, 9.65968391D+4,-3.49596027D+3,                    
     D -1.82751044D+5, 3.91386300D+5,-3.03639248D+5, 1.15782417D+5,     
     E -4.24399280D+3,-6.11477247D+3, 3.66982647D+2/                    
      DATA EGS1 /                                                       
     2 1.927813D5, 7.666859D5, 6.628436D5, 1.586504D5,-7.786476D3,      
     3-4.499687D5,-1.784644D6,-1.546968D6,-4.020658D5,-3.929522D3,      
     4 4.667741D5, 1.849838D6, 1.641313D6, 5.229787D5, 5.928137D4,      
     5-3.017927D5,-1.206483D6,-1.124685D6,-4.478641D5,-8.682323D4,      
     6 1.226517D5, 5.015667D5, 5.032605D5, 2.404477D5, 5.603301D4,      
     7-1.752824D4,-7.411621D4,-7.989019D4,-4.175486D4,-1.024194D4/      
      DATA EGS2 /                                                       
     1-6.459162D5,-2.903581D6,-3.048551D6,-1.004411D6,-6.558220D4,      
     2 1.469853D6, 6.564615D6, 6.843078D6, 2.280839D6, 1.802023D5,      
     3-1.435116D6,-6.322470D6,-6.531834D6,-2.298744D6,-2.639612D5,      
     4 8.665296D5, 3.769159D6, 3.899685D6, 1.520520D6, 2.498728D5,      
     5-3.302885D5,-1.429313D6,-1.512075D6,-6.744828D5,-1.398771D5,      
     6 4.958167D4, 2.178202D5, 2.400617D5, 1.167815D5, 2.663901D4/      
      DATA EGS3 /                                                       
     1 3.117030D5, 1.195474D6, 9.036289D5, 6.876190D4,-6.814556D4,      
     2-7.394913D5,-2.826468D6,-2.152757D6,-2.459553D5, 1.101414D5,      
     3 7.918994D5, 3.030439D6, 2.412611D6, 5.228065D5, 8.542465D3,      
     4-5.421004D5,-2.102672D6,-1.813959D6,-6.251700D5,-1.184348D5,      
     5 2.370771D5, 9.459043D5, 9.026235D5, 4.116799D5, 1.001348D5,      
     6-4.227664D4,-1.738756D5,-1.795906D5,-9.292141D4,-2.397528D4/      
      DATA EGS4 /                                                       
     1-1.072763D5,-5.973532D5,-6.151814D5, 7.371898D4, 1.255490D5,      
     2 2.298769D5, 1.265001D6, 1.252798D6,-2.306276D5,-2.845824D5,      
     3-2.093664D5,-1.100874D6,-1.009313D6, 2.705945D5, 2.506562D5,      
     4 1.274613D5, 6.190307D5, 5.262822D5,-1.336039D5,-1.115865D5,      
     5-5.715764D4,-2.560989D5,-2.228781D5,-3.222789D3, 1.575670D4,      
     6 1.189447D4, 5.161815D4, 4.870290D4, 1.266808D4, 2.069603D3/      
C-----------------------------------------------------------------------
C     THE PROGRAM STARTS HERE                                           
C                                                                       
      IF (IZ.LT.19 .OR. IZ.GT.111) GOTO 900                             
      IF(IZ.GT.102  .AND. IL.GT.0) GOTO 902                             
      Z=FLOAT(IZ)                                                       
      A=FLOAT(IA)                                                       
      EL=FLOAT(IL)                                                      
      AMIN= 1.2D0*Z + 0.01D0*Z*Z                                        
      AMAX= 5.8D0*Z - 0.024D0*Z*Z                                       
      IF(A .LT. AMIN .OR. A .GT. AMAX) GOTO 910                         
C-----------------------------------------------------------------------
C                                ANGUL.MOM.ZERO BARRIER                 
      AA=2.5D-3*A                                                       
      ZZ=1.D-2*Z                                                        
      ELL=1.D-2*EL                                                      
      BFIS0=0.D0                                                        
      CALL LPOLY(ZZ,7,PZ)                                               
      CALL LPOLY(AA,7,PA)                                               
      DO 10 I=1,7                                                       
      DO 10 J=1,7                                                       
      BFIS0=BFIS0+ELZCOF(J,I)*PZ(J)*PA(I)                               
 10   CONTINUE                                                          
      BFIS=BFIS0                                                        
      SBFIS=BFIS                                                        
      EGS=0.D0                                                          
      SEGS=EGS                                                          
C-----------------------------------------------------------------------
C                               VALUES OF L AT WHICH THE BARRIER        
C                            IS 20%(EL20) AND 80%(EL80) OF L=0 VALUE    
      AMIN2=1.4D0*Z + 0.009D0*Z*Z                                       
      AMAX2=20.D0 + 3.0D0*Z                                             
      IF((A.LT.AMIN2-5.D0 .OR. A.GT.AMAX2+10.D0).AND. IL.GT.0) GOTO 920 
      CALL LPOLY(ZZ,5,PZ)                                               
      CALL LPOLY(AA,4,PA)                                               
      EL80=0.D0                                                         
      EL20=0.D0                                                         
      ELMAX=0.D0                                                        
      DO 20 I=1,4                                                       
      DO 20 J=1,5                                                       
      EL80=EL80+ELMCOF(J,I)*PZ(J)*PA(I)                                 
      EL20=EL20+EMNCOF(J,I)*PZ(J)*PA(I)                                 
 20   CONTINUE                                                          
      SEL80=EL80                                                        
      SEL20=EL20                                                        
C     WRITE(6,25) SEL80,SEL20                                           
 25   FORMAT(' EL80=',F8.1,4X,'EL20=',F8.1)                             
C-----------------------------------------------------------------------
C                               VALUE OF L (ELMAX) WHERE BARRIER DISAPP.
      CALL LPOLY(ZZ,6,PZ)                                               
      CALL LPOLY(ELL,9,PL)                                              
      DO 30 I= 1,4                                                      
      DO 30 J=1,6                                                       
      ELMAX=ELMAX+EMXCOF(J,I)*PZ(J)*PA(I)                               
 30   CONTINUE                                                          
      SELMAX=ELMAX                                                      
C---------------------------------------------------------------------- 
C                              VALUE OF BARRIER AT ANG.MOM.  L          
      IF(IL.LT.1) RETURN                                                
      X=SEL20/SELMAX                                                    
      Y=SEL80/SELMAX                                                    
      IF(EL.GT.SEL20) GOTO 40                                           
C-------------------                                 LOW L              
C                                                                       
      Q=0.2D0/(SEL20**2*SEL80**2*(SEL20-SEL80))                         
      QA=Q*(4.D0*SEL80**3 - SEL20**3)                                   
      QB=-Q*(4.D0*SEL80**2 - SEL20**2)                                  
      BFIS=BFIS*(1.D0 + QA*EL**2 + QB*EL**3)                            
      GOTO 50                                                           
C-------------------                                 HIGH L             
C                                                                       
 40   AJ=(-20.D0*X**5 + 25.D0*X**4 - 4.D0)*(Y-1.D0)**2*Y*Y              
      AK=(-20.D0*Y**5 + 25.D0*Y**4 - 1.D0) * (X-1.D0)**2*X*X            
      Q= 0.2D0/((Y-X)*((1.D0-X)*(1.D0-Y)*X*Y)**2)                       
      QA=Q*(AJ*Y - AK*X)                                                
      QB=-Q*(AJ*(2.D0*Y+1.D0) - AK*(2.D0*X+1.D0))                       
      Z=EL/SELMAX                                                       
      A1=4.D0*Z**5 - 5.D0*Z**4 + 1.D0                                   
      A2=QA*(2.D0*Z+1.D0)                                               
      BFIS=BFIS*(A1 + (Z-1.D0)*(A2 + QB*Z)*Z*Z*(Z-1.D0))                
 50   IF(BFIS.LE.0.0D0) BFIS=0.0D0                                      
      IF(EL.GT.SELMAX) BFIS=0.0D0                                       
      SBFIS=BFIS 
C---------------------------------------------------------------        
C     NOW CALCULATE ROTATING GROUND STATE ENERGY                        
C                                                                       
      IF(EL.GT.SELMAX) RETURN                                           
      DO 70 K=1,4                                                       
         DO 70 L=1,6                                                    
            DO 70 M=1,5                                                 
            EGS=EGS+EGSCOF(M,L,K)*PZ(L)*PA(K)*PL(2*M-1)                 
 70   CONTINUE                                                          
      SEGS=EGS                                                          
      IF(SEGS.LT.0.D0) SEGS=0.0D0                                       
      RETURN                                                            
 900  CONTINUE                                                          
C      WRITE(6,1000)                                                    
      SBFIS=0.0                                                         
C                                                                       
C FOR Z<19 SBFIS SET TO 1.E3                                            
C                                                                       
      IF (IZ.LT.19)  SBFIS = 1.E3                                       
      SEGS=0.0                                                          
      SELMAX=0.0                                                        
      RETURN                                                            
 902  CONTINUE                                                          
C      WRITE(6,1002)                                                    
      SBFIS=0.0                                                         
      SEGS=0.0                                                          
      SELMAX=0.0                                                        
      RETURN                                                            
 910  CONTINUE                                                          
C      WRITE(6,1010) IA,IZ                                              
      SBFIS=0.0                                                         
      SEGS=0.0                                                          
      SELMAX=0.0                                                        
      RETURN                                                            
 920  CONTINUE                                                          
C      WRITE(6,1020) IA,IL                                              
      SBFIS=0.0                                                         
      SEGS=0.0                                                          
      SELMAX=0.0                                                        
      RETURN                                                            
C                                                                       
 1000 FORMAT('**** BARFIT CALLED WITH Z LESS THAN 19 OR',               
     1        ' GREATER THAN 111; BFIS IS SET 0.0 ****'/)               
 1002 FORMAT('**** BARFIT CALLED WITH Z GREATER THAN 102',              
     1        ' AND L NOT EQUAL TO ZERO; BFIS IS SET 0.0 ****'/)        
 1010 FORMAT('**** BARFIT CALLED WITH A =',I3,' OUTSIDE ',              
     1        ' THE ALLOWED VALUES FOR Z = ',I3,' ****'/)               
 1020 FORMAT('**** BARFIT CALLED WITH A =',I3,' OUTSIDE ' ,             
     1        'THE ALLOWED VALUES FOR Z = ',I3/26X,'FOR NONZERO L =',I3 
     2         ,' ****'/)                                               
      END                                                               
C                                                                       
C---------------------------------------------------------------------- 
C                                                                       
      REAL FUNCTION BFMS67(Zms,Ams)
C                                                                       
C     This subroutine calculates the fission barriers                                                                  
C     of the liquid-drop model of Myers and Swiatecki (1967).                                                                 
C     Analytic parameterization of Dahlinger 1982 
C     replaces tables. Barrier heights from Myers and Swiatecki !!!                                                                 
C                                                                       
      REAL*4  Zms,Ams,Nms,Ims,ksims,xms 

      Nms = Ams - Zms
      Ims = (Nms-Zms)/Ams
      ksims= 50.15E0 * (1.- 1.78 * Ims**2)
      xms = Zms**2 / (Ams * ksims)
      ums = 0.368E0-5.057E0*xms+8.93E0*xms**2-8.71*xms**3
      Bfms67 = 0.7322E0*Zms**2/Ams**(0.333333E0)*10.E0**ums

      RETURN
      End
C                                                                       
C---------------------------------------------------------------------- 
C                                                                       
      SUBROUTINE LPOLY(X,N,PL)                                          
C                                                                       
C     THIS SUBROUTINE CALCULATES THE ORDINARY LEGENDRE POLYNOMIALS OF   
C     ORDER 0 TO N-1 OF ARGUMENT X AND STORES THEM IN THE VECTOR PL.    
C     THEY ARE CALCULATED BY RECURSION RELATION FROM THE FIRST TWO      
C     POLYNOMIALS.                                                      
C     WRITTEN BY A.J.SIERK  LANL  T-9  FEBRUARY, 1984                   
C                                                                       
C     NOTE: PL AND X MUST BE DOUBLE PRECISION ON 32-BIT COMPUTERS!      
C                                                                       
      DOUBLE PRECISION PL,X                                             
C                                                                       
       INTEGER*4 I,N                                                    
       DIMENSION PL(20)                                                 
       PL(1)=1.0                                                        
       PL(2)=X                                                          
       DO 10 I=3,N                                                      
       PL(I) = ((2*I-3)*X*PL(I-1)-(I-2)*PL(I-2))/(I-1)                  
 10    CONTINUE                                                         
       RETURN                                                           
       END                                                              
                                                                        
                                                                        
C=======================================================================
        FUNCTION EXPOHAZ(K,T)
C*** TIRAGE ALEATOIRE DANS UNE EXPONENTIELLLE : Y=EXP(-X/T)
         EXPOHAZ=-T*LOG(HAZ(K))
         RETURN
        END
C**********************************************************************
C------------------  DISTRIBUTION DE MAXWELL  -------------------------
        FUNCTION FD(E)
           FD=E*EXP(-E)
        END
C**********************************************************************
C------------------ FONCTION INTEGRALE DE FD(E) -----------------------
        FUNCTION F(E)
           F=1-(E+1)*EXP(-E)
        END
C**********************************************************************
C------------  TIRAGE ALEATOIRE DANS UNE MAXWELLIENNE  ----------------  
        FUNCTION FMAXHAZ(K,T) 
C       T : TEMPERATURE
C------ DECLARATION DES VARIABLES

        SAVE
	
        DIMENSION P(100)
       DIMENSION IY(19)
C ial generateur pour le cascade (et les IY pour eviter les correlations)
       common/hazard/ial,IY

        DATA ITEST/0/
C------ PROGRAMME PRINCIPAL
        IF (ITEST.EQ.1) GO TO 120
C       CALCUL DES P(I) PAR APPROXIMATION DE NEWTON
        P(100)=8.
        X=0.1
        DO I=1,99
20         X1=X-(F(X)-I/100.)/FD(X)
           X=X1
           IF (ABS(F(X)-I/100.).LT.1E-5) GO TO 100
           GO TO 20               
100        P(I)=X
        END DO
        ITEST=1
C------ TIRAGE ALEATOIRE ET CALCUL DU X CORRESPONDANT 
C       PAR REGRESSION LINEAIRE
C120     Y=HAZ(K)
120     CALL RIBM(Y,IY(18))
        I=NINT(Y*100)
        
C Ici on evite froidement les depassements de tableaux....(A.B. 3/9/99)        
        IF(I.EQ.0) GOTO 120

        IF (I.EQ.1) THEN
            X=P(I)*Y*100
        ELSE
            X=(P(I)-P(I-1))*(Y*100-I)+P(I)
        END IF
        FMAXHAZ=X*T
        END
  
C****************************************************************************
C===============================================================================

      	subroutine inipace(RACINE)

      	DIMENSION XM(15)
      	common/pace/DM(500,0:500)
      	character*2 a4

      CHARACTER*80 FILENAME,RACINE,FILEDAT              
 
c        open(unit=22,status='OLD',readonly,shared,
c     &    file='/import/projet/spall/Cascabla/PACE2.DATA')
c      	open(unit=22,status='OLD',readonly,shared,file=
c     1  '[BENLLIURE.SIMULATION.EVAPORA]pace2.DATA')
C    	open(22,file='/TABMASS DATA A',form='formatted')

      FILEDAT='pace2.data'
      long=0
      i=0
      DO WHILE(long.EQ.0)
      	i=i+1
	IF(RACINE(i:i).EQ.' ') long=i-1
      END DO
      FILENAME=RACINE(1:long)//FILEDAT                                                         
        open(22,file=FILENAME,form='formatted')
        
        rewind 22
      	DO I=1,500
      	 DO J=0,500
      	  DM(I,J)=0.
      	 END DO
        ENDDO
C   	do i=1,261
      	do i=1,262
      	 read(22,1) a4,a,zi,zf
      	 JF=INT(ZF-ZI)+1
      	 READ(22,*)(XM(J),J=1,JF)
      	 KI=INT(ZI)
      	 KF=INT(ZF)
      	 DO J=KI,KF
      	  JJ=J-KI+1
      	  if(jj.gt.15) then
      	   print *,'pb de dimension'
       	   print *,jj,ki,kf,xm(jj)
          end if
          DM(I,J)=XM(JJ)
         END DO
        end do
      	close(unit=22)
1     	format(A2,f3.0,4x,f3.0,4x,f3.0)
2     	format(11f10.0)

      	return
      	end

C===============================================================================
C=	"PACE2"                                                                =
C=  	Cette fonction retourne le defaut de masse du noyau A,Z en MeV
C              Révisée pour a, z flottants 25/4/2002	                       =
C===============================================================================

      	function pace2(a,z)

        REAL*8 A,Z

      	common/pace/DM(500,0:500)
      	parameter(u=931500)

c        write(6,*)' Pace2: A,Z',a,z

      	ii=IDINT(a+0.5)
      	jj=IDINT(z+0.5)
      	if(ii.le.0.or.jj.lt.0) then
         pace2=0.
         return
      	end if
      	if(jj.gt.300) then
         pace2=0.
      	else
         pace2=dm(ii,jj)
      	end if
      	pace2=pace2/1000.
      	if(dm(ii,jj).eq.0.) then
	if(ii.lt.12) then
	pace2=-500.
	else
         call guet(a,z,pace2)
         pace2=pace2-ii*931.5
         pace2=pace2/1000.
     	end if
        end if
      	return
      	end
  
C===============================================================================
C=	"GUET"                                                                 =
C=	TABLE DE MASSES ET FORMULE DE MASSE TIRE DU PAPIER DE BRACK-GUET       =
C=   	Gives the theoritical value for mass excess...                         =
C              Révisée pour x, z flottants 25/4/2002	                       =
C===============================================================================
 
      	subroutine guet(X,Z,FIND)

        REAL*8 X,Z
        DIMENSION Q(0:50,0:70)

	IX=INT(X+0.5)
        IZ=INT(Z+0.5)
        ZZ=IZ
        XX=IX
        find=0.
        AVOL=15.776
        ASUR=-17.22
        AC=-10.24
        AZER=8.
        XJJ=-30.03
        QQ=-35.4
        C1=-.737
        C2=1.28

        IF(IX.GT.7) GO TO 70
        Q(0,1)=939.50
        Q(1,1)=938.21
        Q(1,2)=1876.1
        Q(1,3)=2809.39
        Q(2,4)=3728.34
        Q(2,3)=2809.4
        Q(2,5)=4668.8
        Q(2,6)=5606.5
        Q(3,5)=4669.1
        Q(3,6)=5602.9
        Q(3,7)=6535.27
        Q(4,6)=5607.3
        Q(4,7)=6536.1
        Q(5,7)=6548.3
        FIND=Q(IZ,IX)
        RETURN

70      CONTINUE
        XNEU=XX-ZZ
        SI=(XNEU-ZZ)/XX
        X13=XX**.333
        EE1=C1*ZZ*ZZ/X13
        EE2=C2*ZZ*ZZ/XX
        AUX=1.+(9.*XJJ/4./QQ/X13)
        EE3=XJJ*XX*SI*SI/AUX
        EE4=AVOL*XX+ASUR*(XX**.666)+AC*X13+AZER
        TOTA=EE1+EE2+EE3+EE4
        FIND=939.55*XNEU+938.77*ZZ-TOTA

        RETURN 
        END
C=======================================================================
                                                                        
