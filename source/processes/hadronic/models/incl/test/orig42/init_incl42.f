
      SUBROUTINE INIT_INCL(INIT_GRAINE)
      
C**************************************************************************
C Subroutine for initialisation of intranuclear cascade INCL
C
C     This will  read some specific parameters for INCL,
C                prepare the Saxon-Wood density for each nucleus
C                compute the deuteron momentum space density from Paris pot.
C                print some global informations 
C**************************************************************************
C input: should contain Z and A for NBMAT nucleus considered in this problem 
      INTEGER ZMAT,AMAT     
      COMMON/MAT/ZMAT(500),AMAT(500),BMAX_GEO(6,500),NBMAT
C input: should contain a seed (ial, odd and of 5 digits) to start the work.     
      COMMON/hazard/ial,igraine(19)

C Dialogue with INCL for nucleus density and parameters.
      COMMON/WS/R0,ADIF,RMAXWS,DRWS,NOSURF,XFOISA,NPAULSTR,BMAX
C Dialogue with INCL: function R(q/pf) for each nucleus
      COMMON/SAXW/ XX(30,500),YY(30,500),SS(30,500),NBPINTER,JMAT

C Input for INCL, to be known at the CALL PNUC (SUPPRIMER icoup!!)
      COMMON/CALINCL/F(15),icoup

C For the 19 secondary seeds of /hazard/:
      DIMENSION nbtirhaz(19)
      DATA nbtirhaz/38,82,76,18,39,31,41,59,26,54,
     s              14,84,13,15,91,89,10,6,52/
     
C FILE des donnees ecrites en dur pour WORSHOP-DEBUG
C	OPEN(5,file='cu4KHS3.in',status='old')
C	OPEN(6,file='cu4KHS3.out',status='old')	
cjcd
      common /inout/ in, io, itty, iscrt                                
cjcd
      write(io,*)' '	
      write(io,*)'*************** VERSION INCL 4.2*********************'
      write(io,*)'* stopping time and potential can be changed     ****'
      write(io,*)'* input of first random numbers                  ****'
      write(io,*)'* bimpact is output instead of sepa as in INCL 3.0 **'
      write(io,*)'* implementation of surface W.S.          4/2000  ***'
      write(io,*)'* interaction only with "participants"    4/2000  ***'
      write(io,*)'* CDPP:Coherent Dynamical Pauli Principle 5/2001  ***'
      write(io,*)'* Paris momentum density for the deuteron 4/2001  ***'
      write(io,*)'* ND-NN (*3) and phase space Delta width  4/2001  ***'
      write(io,*)'************* INCL 4.0 -> INCL4.1 *******************'
      write(io,*)'* No lower cut on pi-N interaction	   2/2002  ****'       
      write(io,*)'* Init of first avatars for participants  2/2002  ***'     
      write(io,*)'* Output of excit energy for absorption   2/2002  ***'
      write(io,*)'************* INCL 4.1 -> INCL4.2 *******************'
      write(io,*)'* Increased absorption below 100 MeV 	   5/2002  ****'       
      write(io,*)'* Coulomb transmis. on projo from LAHET  5/2002  ****'       
      write(io,*)'* light targets(gaussian A<5 MHO 4<A<19) 6/2002  ****'       
      write(io,*)'*****************************************************'
      write(io,*)'* Corrections 11/2003 numerical from J Hendricks    *'
      write(io,*)'*           L of the remnant (from actual c.m.)     *'
      write(io,*)'*****************************************************'
      write(io,*)''
      
C specific parameters for INCL:	 
c espace de phases test (R et P) pour PAULI: 
C Valeur recommandee par J.C. V-test=0.592 h**3:
      rbl=2.
      pbl=200.
C Valeur pour avoir V-test=2 h**3 (avec pbl=200)
      rbl=3.1848
C                                                                       P-N02780
CCC   CONSTANTS AND DERIVED DATA                                        P-N02790
C                                                                       P-N02800
      HC=197.328
      FMP=938.2796
      PF=1.37*HC                                                        P-N02820
      TF=SQRT(PF*PF+FMP*FMP)-FMP
      BINDING=f(5)-TF
      write(io,*)'V0 nuclear potential (A>5) ',f(5),'(MeV) and scaling',
     s  ' stopping time factor: ',f(6)
      write(io,*) 'Fermi momentum: ',PF,' Binding energy: ',BINDING
      write(io,*) 'R*P cell for pauli stat: ',RBL,PBL
      write(io,*) 'Pauli strict (1) or statistic (0): YOUR CHOICE: '
     s                                                 ,NPAULSTR
      write(io,*) ' '
      write(io,*) ' Your choice, NOSURF=  ',NOSURF,'  means'
      write(io,*) 'NOSURF=-2, with W.S. density and INCL4 stopping time'
      write(io,*) 'NOSURF=1, sharp surface'
      write(io,*) 'NOSURF=0, with W.S. density, time without B dep'
      write(io,*) 'NOSURF=-1, with W.S. density, time with B dep.'
      write(io,*) 'RMAXWS=RO+XFOISA*A, XFOISA= ',XFOISA,' above A=19'
      
      K1=5                                                              P-N02530
      K2=0                                                              P-N02540
      k3=0
      k4=2
      k5=1
      k6=0
      write(io,*)' '
      write(io,*) 'K1,K2,K3,K4,K5,K6 ',K1,K2,K3,K4,K5,K6
      write(io,*) '   Meaning: K1=5, REFERENCE FRAME=LAB SYSTEM'
      write(io,*) '            K2=0, RELATIVISTIC KINEMATICS'
      write(io,*) '            K3=0, DELTAS ARE PRODUCED'
      write(io,*) '            K4=2, DELTA HAS A EXPONENTIALLY'
     s                                     ,' RANDOM LIFETIME'
      write(io,*) '            K5=1, DELTA-NUCLEON=DELTA-DELTA'
     s ,'=NUCLEON-NUCLEON ELASTIC X-SECTION'
C      write(io,*) '            K6=0, NO ANGULAR MOMENTUM CONSERVATION'
      write(io,*)' '
        
c-----------------------------LD's MODIF, 7/8/2001-------------------
c      ZMAT(1)=f(2)+0.5
c      AMAT(1)=f(1)+0.5
c      NBMAT=1      
c--------------------------------------------------------------------
C*************************************************************************
C preparation of 19 other seeds (can also be initialized from outside):
      IF(INIT_GRAINE.EQ.0) GO TO 496
      ialdep=ial
      DO i=1,19
         DO j=1,nbtirhaz(i)
           CALL RIBM(xrand,ial)
         ENDDO
5        CALL RIBM(xrand,ial)
	 IF(xrand.EQ.0.) GO TO 5
	 xrand = xrand*100000
10       CONTINUE
	 IF (xrand.lt.10000) THEN
	      xrand=xrand*10
	      GO TO 10
	 ELSE
	      igraine(i)=xrand
	      IF(igraine(i).EQ.(igraine(i)/2)*2) igraine(i)=igraine(i)+1
         ENDIF
      ENDDO
      ial=ialdep
C      write(6,1000) ialdep
1000  FORMAT('  Primary:', i8, ' and subsidiary seeds for INCL')
C	 write(6,1001) (igraine(i),i=1,6) 
C	 write(6,1001) (igraine(i),i=7,12) 
C	 write(6,1001) (igraine(i),i=13,19)
1001  FORMAT(7(i8,2x)) 
496   CONTINUE
C*************************************************************************
C Calculation with realistic nuclear density (Saxon-Wood)
	IF (NOSURF.LE.0)  THEN
	 
C prepare nucleus density for NBMAT nucleus defined in common/MAT/
      IF(NBMAT.GT.500) THEN
            write(6,1002) NBMAT	 
1002  FORMAT(' You need ',i5,' nucleus in your problem. This overpass',/,
     s ' the dimension of 20 in commons /MAT/ and /SAXW/',/,' You have', 
     s ' to enlarge dimensions before making this calculation.')
            STOP
      ENDIF
      
      DO i=1,NBMAT
         imat=i
	 IZMAT=ZMAT(i)
	 IAMAT=AMAT(i)
         CALL INIT_MAT(IZMAT,IAMAT,imat)
	 
      ENDDO
        ENDIF
	
C**************************************************************************	
C Deuteron density in momentum space:
      CALL DENS_DEUT

    
      RETURN
      END
      
C**************************************************************************
C Computes the R(q/pf) function for a given nuclei Z,A and put it
C   in the common/SAXW/
C**************************************************************************

      SUBROUTINE INIT_MAT(IZMAT,IAMAT,IMAT)
      
      INTEGER ZMAT,AMAT,IZMAT,IAMAT     
      COMMON/MAT/ZMAT(500),AMAT(500),BMAX_GEO(6,500),NBMAT      
C Dialogue with INCL for nucleus density and parameters.
      COMMON/WS/R0,ADIF,RMAXWS,DRWS,NOSURF,XFOISA,NPAULSTR,BMAX
C Dialogue with INCL: function R(q/pf) for each nucleus
      COMMON/SAXW/ XX(30,500),YY(30,500),SS(30,500),NBPINTER,JMAT
      
      
C RMS espace R, espace P, Fermi momentum and energy for light gauss nuc.      
      COMMON/light_gaus_nuc/rms1t(9),pf1t(9),pfln(9),tfln(9),vnuc(9)
      DATA rms1t,pf1t/0.,0.,0.,0.,0.,2.10,1.80,1.80,1.63,
     -0.,0.,0.,0.,0.,77.,110.,110.,153./


      
      
C Fermi 2 param from A=19 to 28, modified harm oscil A=6 to 18
C (H. De Vries et al. At. Data and Nuc. Data Tab. 36 (1987) 495)
      COMMON/light_nuc/R_light_nuc(30),a_light_nuc(30)
      DATA R_light_nuc/0.,0.,0.,0.,0.,0.334,0.327,0.479,0.631,0.838,
     s 0.811,1.07,1.403,1.335,1.25,1.544,1.498,1.513,2.58,2.77,
     s 2.775,2.78,2.88,2.98,3.22,3.03,2.84,3.14,0.,0./
      DATA a_light_nuc/0.,0.,0.,0.,0.,1.78,1.77,1.77,1.77,1.71,
     s 1.69,1.69,1.635,1.730,1.81,1.833,1.798,1.841,0.567,0.571,
     s 0.560,0.549,0.550,0.551,0.580,0.575,0.569,0.537,0.,0./


      
cjcd
      common /inout/ in, io, itty, iscrt                                
cjcd

      EXTERNAL WSAX,DERIVWSAX,DMHO,DERIVMHO,DERIVGAUS

C print of the function:
      write(io,*) '***************************************************'
      write(io,1789) IZMAT,IAMAT,IMAT
1789  FORMAT(' Nuclear density for nucleus Z,A: ',2i5,'   imat=',i3)


      FMP=938.2796	!From INCL data

C parametres moyens de densite de la cible (fermi 2 parametres)
      IF (IAMAT.GE.28) THEN
      	R0 = (2.745E-4*IAMAT+1.063)*IAMAT**0.33333333
      	ADIF = 1.63E-4*IAMAT+0.510
        RMAXWS = R0 + XFOISA*ADIF
      ELSE IF(IAMAT.GE.19) THEN
      	R0 = R_light_nuc(IAMAT)
      	ADIF = a_light_nuc(IAMAT)
        RMAXWS = R0 + XFOISA*ADIF
      ELSE IF(IAMAT.GE.6) THEN
      	R0 = R_light_nuc(IAMAT)
      	ADIF = a_light_nuc(IAMAT)
	RMAXWS = 5.5 + 0.3*(IAMAT-6.)/12.
      ELSE IF(IAMAT.GE.2) THEN
      	IF(IAMAT.EQ.2) THEN
		R0=rms1t(6)
		PFLN(6)=pf1t(6)*1.291   !SQRT(5/3)=1.291
	        TFLN(6)=SQRT(PFLN(6)**2+FMP*FMP)-FMP
		VNUC(6)=TFLN(6)+2.22
		WRITE(io,1806) VNUC(6),PFLN(6),TFLN(6)
	ENDIF
      	IF(IAMAT.EQ.3.AND.IZMAT.EQ.1) THEN
		R0=rms1t(7)
		PFLN(7)=pf1t(7)*1.291   !SQRT(5/3)=1.291
	        TFLN(7)=SQRT(PFLN(7)**2+FMP*FMP)-FMP
		VNUC(7)=TFLN(7)+4.24
		WRITE(io,1806) VNUC(7),PFLN(7),TFLN(7)
	ENDIF
      	IF(IAMAT.EQ.3.AND.IZMAT.EQ.2) THEN
		R0=rms1t(8)
		PFLN(8)=pf1t(8)*1.291   !SQRT(5/3)=1.291
	        TFLN(8)=SQRT(PFLN(8)**2+FMP*FMP)-FMP
		VNUC(8)=TFLN(8)+3.86
		WRITE(io,1806) VNUC(8),PFLN(8),TFLN(8)
	ENDIF
      	IF(IAMAT.EQ.4) THEN
		R0=rms1t(9)
		PFLN(9)=pf1t(9)*1.291   !SQRT(5/3)=1.291
	        TFLN(9)=SQRT(PFLN(9)**2+FMP*FMP)-FMP
		VNUC(9)=TFLN(9)+9.43
		WRITE(io,1806) VNUC(9),PFLN(9),TFLN(9)
	ENDIF
	ADIF=0.57735*R0
	RMAXWS=R0+2.5
      END IF
      DRWS = RMAXWS/29.
1806  FORMAT(' nuclear pot:',F8.2,' Mev, fermi momentum and energy',
     s F8.2,' MeV/c',F8.2,' MeV')
C Bmax for sigma geom and various projectiles (p,n,pion/d/t/He3/He4/)
      DO i=1,6
      	 j=i
      	 IF(i.GT.2) j=i+3
         BMAX_GEO(i,IMAT)=RMAXWS+rms1t(j)
      END DO
C preparation de la distribution W.S.:
      IF (IAMAT.GE.19) THEN
      step=0.2
        CALL INTEG(0.,13.5,step,DERIVWSAX,RES_DWS)
      ELSE 
C preparation de la distribution M.H.O.:
      	IF(IAMAT.GE.6) THEN
      step=0.1
        CALL INTEG(0.,10.,step,DERIVMHO,RES_DWS)
      	ELSE
C preparation de la distribution Gaussienne:
        CTE=ADIF**3*SQRT(2.*3.141592654)        
	RES_DWS=3.*CTE/2.
      	END IF
      END IF
      FNOR=RES_DWS
C calcul de q/PF=F(R)      
      nbr=RMAXWS/DRWS + 1.5
      rcour=-DRWS
      j=0
      DO i=1,nbr
      	rcour=rcour+DRWS
      	IF(i.EQ.1) THEN
      		j=j+1
      		F_R=0.
      		XX(j,IMAT)=F_R
      		YY(j,IMAT)=0.		!On impose x(1)=0., y(1)=0.
      		RES_DWS=0.
      	ELSE
      		step=rcour/20.
      		IF(step.GT.0.05) step=0.05
	      IF (IAMAT.GE.19) THEN
      		CALL INTEG(0.,rcour,step,DERIVWSAX,RES_DWS)
		F_R=RES_DWS/FNOR
	      ELSE 
	      	IF(IAMAT.GE.6) THEN
		CALL INTEG(0.,rcour,step,DERIVMHO,RES_DWS)
		F_R=RES_DWS/FNOR
	      	ELSE 
		CALL INTEG(0.,rcour,step,DERIVGAUS,RES_DWS)
		F_R=RES_DWS/FNOR
	      	END IF
	      END IF    		
C Modif le 20/10/2003; eviter les valeurs negatives AVANT **1/3 !
C      		F_R=F_R**(1./3.)
      	ENDIF
      		IF(F_R.GT.0.0) THEN
      		F_R=F_R**(1./3.)
      			j=j+1
      			XX(j,IMAT)=F_R
      			YY(j,IMAT)=rcour
      		ENDIF
      END DO
      nbpinter=j
      XX(j,IMAT)=1.			!On impose x(nbpinter)=1. (y()=rmax)

C interpolation de F_inv(R) (fonction inverse de F(R))       
      CALL FLIN2(IMAT)
      
      
      IF(IAMAT.GE.19) write(io,1870) R0,ADIF
1870  FORMAT(' Wood-Saxon density, R0= ',F7.3,5x,'A= ',F7.3)
      IF(IAMAT.GE.6.AND.IAMAT.LT.19) write(io,1871) R0,ADIF     
1871  FORMAT(' Modif Harm. Oscil. density, alpha= ',F7.3,5x,'A= ',F7.3)
      IF(IAMAT.GE.2.AND.IAMAT.LT.6) write(io,1872) R0,ADIF     
1872  FORMAT(' Gaussian density, R.M.S.= ',F7.3,5x,'Sigma= ',F7.3)

      GEOM=31.41592653*RMAXws**2
      write(io,1873) RMAXws,GEOM
1873  FORMAT('FOR INCIDENT Nucleons OR Pions, RMAX=',F8.3,/,
     s ' and geometrical (PI*RMAXws**2) reaction cross',
     s ' section (mb) is: ',F10.3) 
      write(io,1792) (BMAX_GEO(k,IMAT),k=3,6)    
1792  FORMAT(' RMAXws for d/t/3He/4He: ',4F10.3)
      write(io,*) ' '
c      write(io,1790)
1790  FORMAT(1x,'Exact calculation of the R(q) function for the ',
     s'target nucleus density',/,'    q/PF        R(q/PF)')

C      DO i=1,nbpinter
C      		write(io,1791) XX(i,IMAT),YY(i,IMAT)
C      ENDDO            
1791  FORMAT(F10.6,F10.2)
C      write(io,*) 'lin. interpolation prepared for',nbpinter,' points'
      write(io,*) '***************************************************'
      
      RETURN
      END
c-------------------------------------------------------------------------------


      FUNCTION WSAX(R)
      COMMON/WS/R0,ADIF,RMAXWS,DRWS,NOSURF,XFOISA,NPAULSTR,BMAX
      WSAX = R**2/(1.+EXP((R-R0)/ADIF)) 
      RETURN
      END
           
      FUNCTION DERIVWSAX(R)
      COMMON/WS/R0,ADIF,RMAXWS,DRWS,NOSURF,XFOISA,NPAULSTR,BMAX
      DERIVWSAX=R**3*EXP((R-R0)/ADIF)/(1.+EXP((R-R0)/ADIF))**2
      DERIVWSAX=DERIVWSAX/ADIF
      RETURN
      END
      
      FUNCTION DMHO(R)
      COMMON/WS/R0,ADIF,RMAXWS,DRWS,NOSURF,XFOISA,NPAULSTR,BMAX
      ARG=(R/ADIF)**2
      DMHO=R*R*(1.+R0*ARG)*exp(-ARG)
      RETURN
      END
      
      FUNCTION DERIVMHO(R)
      COMMON/WS/R0,ADIF,RMAXWS,DRWS,NOSURF,XFOISA,NPAULSTR,BMAX
      ARG=(R/ADIF)**2
      DERIVMHO=-2.*R*R*ARG*(R0 -1.-R0*ARG)*exp(-ARG)
      RETURN
      END

      FUNCTION DERIVGAUS(R)
      COMMON/WS/R0,ADIF,RMAXWS,DRWS,NOSURF,XFOISA,NPAULSTR,BMAX
      ARG=(R/ADIF)**2
      DERIVGAUS=R*R*ARG*exp(-ARG/2.)      
      RETURN
      END
     
      SUBROUTINE INTEG(AMI,AMA,DR,FONC,RES)                             INT00010
C  SOUS-PROGRAMME D'INTEGRATION PAR LA METHODE D'ALKHAZOV               INT00020
      DIMENSION X1(5)                                                   INT00030
      PAS=DR                                                            INT00040
      RI=AMI                                                            INT00050
      RA=AMA                                                            INT00060
      ACONT=1.                                                          INT00070
      IF(AMA.GT.AMI)GO TO 1                                             INT00080
      ACONT=-1.                                                         INT00090
      RI=AMA                                                            INT00100
      RA=AMI                                                            INT00110
1     CONTINUE                                                          INT00120
      X1(1)=95./288.                                                    INT00130
      X1(2)=317./240.                                                   INT00140
      X1(3)=23./30.                                                     INT00150
      X1(4)=793./720.                                                   INT00160
      X1(5)=157./160.                                                   INT00170
      NB=(RA-RI)/DR+1.0000000001                                        INT00180
      DR=(RA-RI)/(NB-1)                                                 INT00190
      RES=0.                                                            INT00200
      IF (NB.LT.10) GO TO 100                                           INT00210
      DO 10 I=1,5                                                       INT00220
      RES=RES+(FONC(RI)+FONC(RA))*X1(I)                                 INT00230
      RI=RI+DR                                                          INT00240
      RA=RA-DR                                                          INT00250
10    CONTINUE                                                          INT00260
      NB=NB-10                                                          INT00270
      IF (NB.EQ.0) GO TO 20                                             INT00280
      DO 11 I=1,NB                                                      INT00290
      RES=RES+FONC(RI)                                                  INT00300
      RI=RI+DR                                                          INT00310
11    CONTINUE                                                          INT00320
20    CONTINUE                                                          INT00330
      RES=RES*DR*ACONT                                                  INT00340
      DR=PAS                                                            INT00350
      RETURN                                                            INT00360
100   WRITE (6,101)                                                     INT00370
101   FORMAT(1H ,33HPAS ASSEZ DE POINTS D|INTEGRATION)                  INT00380
      DR=PAS                                                            INT00390
      RETURN                                                            INT00400
      END                                                               INT00410

      SUBROUTINE FLIN2(K)                                   
      COMMON/SAXW/ X(30,500),Y(30,500),S(30,500),N,JMAT
      COMMON/WS/R0,ADIF,RMAXWS,DRWS,NOSURF,XFOISA,NPAULSTR,BMAX
C PREMIER SP A APPELLER POUR METHODE D'INTERPOLATION LINEAIRE
C X(1:N) CONTIENT LES ABSISSES , Y LES ORDONNEES                       
C FLIN2 CALCULE LES DERIVEES premieres QU'IL MET DANS S               
      DO i=1,N-1
      		S(i,K)=(Y(i+1,K)-Y(i,K))/(X(i+1,K)-X(i,K))
      ENDDO
      S(N,K)=S(N-1,K)
      RETURN
      END
      FUNCTION FLIN(XV)
      COMMON/SAXW/ X(30,500),Y(30,500),S(30,500),N,K
      COMMON/WS/R0,ADIF,RMAXWS,DRWS,NOSURF,XFOISA,NPAULSTR,BMAX
C FONCTION D'INTERPOLATION AU POINT XV ( MEME HORS BORNES )             
C DE LA FN X->Y DONT LES DERIVEES premieres (S) ONT ETE                 
C EVALUEES PAR L'APPEL PREALABLE DE FLIN2                              
C LES INDICES VONT DE 1 A N                                             
      TZ=XV-X(1,K)                                                        
      IF(TZ)1,2,3                                                       
1     FLIN=Y(1,K)+S(1,K)*TZ
      call hfill(201, FLIN, 0.0, 1.0);
      call hfill(200, FLIN, XV, 1.0)
      RETURN                                                            
2     FLIN=Y(1,K)                                                       
      call hfill(201, FLIN, 0.0, 1.0);
      call hfill(200, FLIN, XV, 1.0)
      RETURN                                                           
3     CONTINUE
                                                       
      DO 10 i=2,N
                j=i
      		TZ=XV-X(i,K)                                  
      		IF(TZ)8,9,10                                                     
10    CONTINUE
      GO TO 8 
                                                        
9     FLIN=Y(j,K)                                                     
      call hfill(201, FLIN, 0.0, 1.0);
      call hfill(200, FLIN, XV, 1.0)
      RETURN
                                                                
8     j=j-1                                                             
      DGX=XV-X(j,K)
      FLIN=Y(j,K)+S(j,K)*DGX                   
      call hfill(201, FLIN, 0.0, 1.0);
      call hfill(200, FLIN, XV, 1.0)
      RETURN                                                            
      END                                                               
                                          
C ********************************************************************
c-------------------------------------------------------------------------------
      SUBROUTINE DENS_DEUT
C Ce subroutine appele sur le premier tir va calculer la densite du deuton
C   dans l'espace des impulsions et preparer l'interpolation permettant ensuite
C   le tir au hasard d'un module de l'impulsion (q).
C Ce subroutine remplit le common /SPL2/:
C          XSP(0:1), YSP integrale normalisee de la densite de 0 a q.
C          A(),B(),C() coefs des NSP points pour une interpolation du second degre.
C   q est en fm-1. 
      EXTERNAL DENS
      DIMENSION Q(100),F(100)
      COMMON/SPL2/ XSP(100),YSP(100),A(100),B(100),CC(100),NBP
      
      REAL*8 C(13),D(13),FN                                             DEU00350
      COMMON /DTON/C,D,FN                                               DEU00370
      DATA C/.88688076D+00,-.34717093D+00,-.30502380D+01,               DEU00380
     1 .56207766D+02,-.74957334D+03,.53365279D+04,-.22706863D+05,       DEU00390
     2 .60434469D+05,-.10292058D+06,.11223357D+06,-.75925226D+05,       DEU00400
     3 .29059715D+05,-.48157368D+04/                                    DEU00410
      DATA D/.23135193D-01,-.85604572D+00,.56068193D+01,                DEU00420
     1 -.69462922D+02,.41631118D+03,-.12546621D+04,.12387830D+04,       DEU00430
     2 .33739172D+04,-.13041151D+05,.19512524D+05,-.15634324D+05,       DEU00440
     3 .66231089D+04,-.11698185D+04/                                    DEU00450
      DATA FN/.28212D+00/                                               DEU00460
C AVEC FN=.28212 LES FO RADIALES SUIVANTES SONT NORMALISEES A:          DEU00470
C SOMME(0,INFINI)(DEUT0(Q)**2 + DEUT2(Q)**2))*Q*Q*DQ = 1./4*PI          DEU00480
C ET CECI DANS L'ESPACE R ET DANS L'ESPACE Q. PD=5.74%                  DEU00490
cjcd
      common /inout/ in, io, itty, iscrt                                
cjcd
      DQ=0.01
      Q(1)=0.
      DO i=2,50
      	  Q(i)=Q(i-1)+DQ
      ENDDO
      
      NBP=77   !nombre de points de calcul

      DQ=0.1
      DO i=51,NBP
      	  Q(i)=Q(i-1)+DQ
      ENDDO
      
      F(1)=0.
      SUMINT=0.
      DO i=2,NBP
        DQ=(Q(i)-Q(i-1))/10.
        CALL INTEG(Q(i-1),Q(i),DQ,DENS,RES)
	SUMINT=SUMINT+RES
        F(i)=SUMINT
C	write(6,*) i,Q(i),SUMINT,DENS(Q(i))
      ENDDO

      DO i=1,NBP
        XSP(i)=F(i)/F(NBP)
	YSP(i)=Q(i)
      ENDDO
      CALL SPL2AB

      write(io,100) NBP,Q(NBP)
100   FORMAT('  Deuteron density in q space from PARIS potential',/,
     s '   ',I4,' exact values from 0 to ',F6.2,' fm-1')
           
      RETURN
      END
C-----------------------------------------------------------------------------------
      FUNCTION DEUTV(L,Q)                                               DEU00510
C FONCTION D'ONDE DU DEUTON ESPACE Q DE VINH MAU P.L.101B,139.          DEU00520
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*4 Q,DEUTV
      DIMENSION C(13),D(13)                                             DEU00540
      COMMON /DTON/C,D,FN                                               DEU00550
      PI = 3.141592654                                                  DEU00560
      Q2 = Q*Q                                                          DEU00570
      RES = 0.                                                          DEU00580
      IF (L.NE.0) GO TO 51                                              DEU00590
      DO 52 I=1,13                                                      DEU00600
52    RES = RES + C(I)/(Q2+FM2(I))                                      DEU00610
      GO TO 53                                                          DEU00620
51    CONTINUE                                                          DEU00630
      DO 54 I=1,13                                                      DEU00640
54    RES = RES + D(I)/(Q2+FM2(I))                                      DEU00650
53    DEUTV = RES*DSQRT(2./PI)*FN                                       DEU00660
      RETURN                                                            DEU00670
      END                                                               DEU00680

      FUNCTION FM2(J)                                                   DEU00910
      IMPLICIT REAL*8(A-H,O-Z)                                          DEU00920
      DATA A/0.23162461/                                                DEU00930
      FM2 = (A + (J - 1))*(A + (J - 1))                                 DEU00940
      RETURN                                                            DEU00950
      END
      
      FUNCTION DENS(Q)
C      REAL*8 DEUTV,DENS,Q
      DENS=Q*Q*(DEUTV(0,Q)**2+DEUTV(2,Q)**2)
      RETURN
      END      
C-------------------------------------------------------------------------------------      
      SUBROUTINE SPL2AB
C      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/SPL2/ X(100),Y(100),A(100),B(100),C(100),N    
      DO i=1,N-2
        j=i+1
	k=i+2
       
        C(i)=(Y(k)-Y(i))*(X(j)-X(i))-(X(k)-X(i))*(Y(j)-Y(i))
        C(i)=C(i)/((X(j)-X(i))*(X(k)-X(i))*(X(k)-X(j)))
	
	B(i)=(Y(j)-Y(i))/(X(j)-X(i))
	
	A(i)=Y(i)
      ENDDO
      DO i=N-1,N
        C(i)=C(N-2)
	B(i)=B(N-2)
	A(i)=A(N-2) 
      ENDDO     
      RETURN
      END
      FUNCTION SPLINEAB(XV)
C      IMPLICIT REAL*8(A-H,O-Z)
      REAL*4 XV
      COMMON/SPL2/ X(100),Y(100),A(100),B(100),C(100),N
      TZ=XV-X(1)                                                        
      IF(TZ)1,2,3                                                       
1     SPLINEAB=A(1)+B(1)*TZ+C(1)*TZ*(XV-X(2))
      RETURN                                                            
2     SPLINEAB=Y(1)                                                       
      RETURN                                                           
3     CONTINUE
      DO 10 i=2,N-1
                j=i
      		TZ=XV-X(i)                                  
      		IF(TZ)8,9,10                                                     
10    CONTINUE
      GO TO 8 
                                                        
9     SPLINEAB=Y(j)                                                     
      RETURN
                                                                
8     j=j-1 
      TZ=XV-X(j)                                                            
      SPLINEAB=A(j)+B(j)*TZ+C(j)*TZ*(XV-X(j+1))                                  
      RETURN
      END

