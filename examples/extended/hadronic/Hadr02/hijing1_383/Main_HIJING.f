      program test
      CHARACTER FRAME*8,PROJ*8,TARG*8
      COMMON/HIPARNT/HIPR1(100),IHPR2(50),HINT1(100),IHNT2(50)
      SAVE  /HIPARNT/
      COMMON/HIMAIN1/NATT,EATT,JATT,NT,NP,N0,N01,N10,N11
      SAVE  /HIMAIN1/
      COMMON/HIMAIN2/KATT(130000,4),PATT(130000,4)
      SAVE  /HIMAIN2/
      COMMON/RANSEED/NSEED
      SAVE  /RANSEED/
C
C       
C....initialize HIJING for Au+Au collisions at c.m. energy of 200 GeV:
c        FRAME="CMS"
c        PROJ="A"
c	TARG="A"
        IAP=1     
	IZP=1       
	IAT=197      
	IZT=79    
        EFRM=200.0   !GeV 
        N_EVENT=100
        IHPR2(10)=0
c       ---------------------------------------------------------
c        initialize HIJING for geant4
	IHNT2(1)=IAP
	IHNT2(2)=IZP
	IHNT2(3)=IAT
	IHNT2(4)=IZT
	IHNT2(5)=211   !pi+
	IHNT2(6)=0
cc     --------------Ap>1----------------------
         
        HINT1(8)=MAX(ULMASS(2112),ULMASS(2212))
        HINT1(9)=HINT1(8)
        write(*,*)'mass at the beginning',HINT1(8), HINT1(9)
cc---------------------------------------------
        if(IHNT2(5).ne.0) then
        HINT1(8)=ULMASS(IHNT2(5))
        endif
	
c        write(*,*)'modified mass',HINT1(8)
c
c
c         ranseed=1097569630
c         print *,"seed",ranseed
c         call sseed(ranseed)

c         do i=1, 10
c         write(*,*)'ran',NSEED, RAN(NSEED), RLU(NSEED)
c         end do
C
C
C*** initialize HIJING
c        CALL HIJSET(EFRM,FRAME,PROJ,TARG,IAP,IZP,IAT,IZT)
         CALL HIJSET(EFRM,IAP,IZP,IAT,IZT)
c         CALL HIJSET(EFRM)
        
          write(*,*)'Sjet=', HINT1(11),'mb','Stot=',HINT1(13),'mb'
C
C....set BMIN=0 and BMAX=0.0 for central interaction
      BMIN=0.0
      BMAX=0.0    !HIPR1(34)+HIPR1(35)
C....generating N_EVENT events of central AA interaction:
      DO 200 IE=1,N_EVENT
c         CALL HIJING(FRAME,BMIN,BMAX)
         CALL HIJING(BMIN,BMAX)
         WRITE(*,*) IE,NATT,EATT
 200  continue
      STOP
      END 

c	FUNCTION RAN(NSEED)
c	RAN=RLU(NSEED)
c	RETURN
c	END
