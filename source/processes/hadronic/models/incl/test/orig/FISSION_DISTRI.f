	SUBROUTINE FISSION_DISTRI(A,Z,E,A1,Z1,E1,A2,Z2,E2)                     

C                                                                       
C*********************************************************************  
C  FISSION MODEL:                                                       
C            BENLLIURE ET AL., NUCL. PHYS. A 628 (1998) 458             
C*********************************************************************  
C                                                                       
C                                                                       
C 4-2-98: Modifications including the charge polarization for symmetric and 
c         asymmetric fission in the right way
c         --> symmetric fission: E = M1 + M2 + Ec(1,2)
c         --> asymmetric fission: sig_z(A=cte) = 0.5 (Lang)
c                                 sig_n(Z=cte) = 2.58*sig_z(A=cte) = 1.29
c                                 sig_z(N=cte) = 1.63*sig_z(A=cte) = 0.8
c
c         New parametrization of the potential stiffness according to
c         Mulgin et al NPA 640 (1998) 375
c
C
C                                                                       
c This program calculates isotopic distributions of fission fragments */
c with a semiempirical model                                          */
c Copy from SIMFIS3, KHS, 8. February 1995                            */
c Modifications made by Jose Benlliure and KHS in August 1996         */
C                                                                       
C The excitation energy at the saddle point is calculated from the      
c lowest potencial well (LDM-shell(82)-shell(89))  J.Benlliure Sept 1996
C                                                                       
C                                                                       
	REAL*4    N,Nlight1,Nlight2,Aheavy1,Aheavy2,Alight1,Alight2            
	REAL*4    Eheavy1,Eheavy2,Elight1,Elight2                              
	REAL*4    Zheavy1,Zheavy2,Bsym,Basym_1,Basym_2	                        
                                                                        
	REAL*4    Zheavy1_shell,Zheavy2_shell                                  
	REAL*4    Zlight1,Zlight2                                              
	REAL*4    MassCurv                                                     
	REAL*4    Sasymm1,Sasymm2,Ssymm,Ysum,Yasymm                            
	REAL*4    Ssymm_mode1,Ssymm_mode2                                      
	REAL*4    CZ_asymm1_saddle,CZ_asymm2_saddle                            
	REAL*4    WZasymm1_saddle, WZasymm2_saddle, WZsymm_saddle              
	REAL*4    WZasymm1_scission, WZasymm2_scission, WZsymm_scission        
	REAL*4    WZasymm1,WZasymm2,WZsymm                                     
	REAL*4    Nlight1_eff, Nlight2_eff                                     
        INTEGER*4 IMODE                                                 
	REAL*4    RMODE                                                        
	REAL*4    Z1mean, Z1width                                              
	REAL*4    Z1,Z2,N1R,N2R,A1R,N1,N2,A1,A2                                
	REAL*4    NE_MIN,NE_M1,NE_M2                                           
	REAL*4    Zsymm,Nsymm,Asymm                                            
	REAL*4    N1mean, N1width                                              
	REAL*4    R_E_O_MAX                                                    
	REAL*4    E_PAIR,E1,E2                                                 
        REAL*4    HAZ,GAUSSHAZ                                          
        INTEGER*4 I_HELP,IZ,KKK                                         
        INTEGER*4 icz,iis                                         
        REAL*4    El,A2R,RNE_M1,RNE_M2,RN1,RN2,ED1_LOW,ED2_LOW
        REAL*4    ED1_HIGH,ED2_HIGH,ED1,ED2,ATOT                    
                                                                        
                                                                        
c Parameters of the semiempirical fission model */                      
                                                                        
	REAL*4    A_levdens                                                    
	REAL*4    A_levdens_light1,A_levdens_light2                            
	REAL*4    A_levdens_heavy1,A_levdens_heavy2                            
	REAL*4    r_null                                                       
	REAL*4    mass,massdef,rr                                                      
	REAL*4    epsilon_1_saddle,epsilon0_1_saddle                           
	REAL*4    epsilon_2_saddle,epsilon0_2_saddle,epsilon_symm_saddle       
	REAL*4    epsilon_1_scission,epsilon0_1_scission                       
	REAL*4    epsilon_2_scission,epsilon0_2_scission                       
        REAL*4    epsilon_symm_scission                                 
	REAL*4    E_eff1_saddle,E_eff2_saddle                                  
                                                                        
	PARAMETER (r_null=1.16)                                                
                                                                        
	PARAMETER (kkk=10)                                                     
                                                                        
c*** Parameters with I: are calculated internally.)                     
c*** Spectrum R_Mode gets only one data point for each start!)          
                                                                        
	REAL*4 Z         ! Nuclear charge number                               
	REAL*4 A         ! Nuclear mass number                                 
	REAL*4 E         ! Excitation energy above fission barrier             
	REAL*4 Nheavy1   ! position of heavy peak valley 1                     
	REAL*4 Nheavy2   ! position of heavy peak valley 2                     
                                                                        
c values calculated from proposition of Gverdovskii et al.,             
c      Sov. J. Nucl. Phys. 55 (1992) 9                                  
                                                                        
                                                                        
c   	REAL*4  Ysymm           !I: Yield of symmetric mode                
   	REAL*4  Yasymm1         !I: Yield of asymmetric mode 1              
   	REAL*4  Yasymm2         !I: Yield of asymmetric mode 2              
   	REAL*4  Nheavy1_eff     !I: Effective position of valley 1          
   	REAL*4  Nheavy2_eff     !I: Effective position of valley 2          
   	REAL*4  Eexc1_saddle    !I: Excitation energy above saddle 1        
   	REAL*4  Eexc2_saddle    !I: Excitation energy above saddle 2        
   	REAL*4  Eexc_max        !I: Excitation energy above lowest saddle   
   	REAL*4  Aheavy1_mean    !I: Effective mass mode 1                   
   	REAL*4  Aheavy2_mean    !I: Effective mass mode 2                   
   	REAL*4  WAsymm_saddle   !I: Width of symmetric mode                 
   	REAL*4  WAheavy1_saddle !I: Width of asymmetric mode 1              
  	REAL*4  WAheavy2_saddle !I: Width of asymmetric mode 2               
   	REAL*4  WAsymm          !I: Width of symmetric mode                 
   	REAL*4  WAheavy1        !I: Width of asymmetric mode 1              
        REAL*4  WAheavy2        !I: Width of asymmetric mode 2          
   	REAL*4  R_E_O           !I: Even-odd effect in Z                    
   	REAL*4  CZ_symm         !I: Curveture of symmetric valley           
   	REAL*4  CN              !I: Curvature of mass distribution for fixed
                                                                        
        REAL*4  N1ucd           ! Neutron number given by UCD                                                     
        REAL*4  N2ucd           ! Neutron number given by UCD                                                     
                                                                        
C	Z=92    ! Nuclear charge number                                       
C	A=236   ! Nuclear mass number                                         
C	E=-2     ! Excitation energy above fission barrier                    
                                                                        
                                                                        
	REAL*4	E_crit                                                          
	REAL*4	Delta_U1_shell                                                  
	REAL*4	Delta_U2_shell                                                  
	REAL*4  CZ_asymm1_shell                                                
	REAL*4  CZ_asymm2_shell                                                
	REAL*4  Fwidth_asymm1                                                  
	REAL*4  Fwidth_asymm2                                                  
	REAL*4  CZ_asymm2_scission                                             
	REAL*4  FGAMMA1                                                        
	REAL*4  FGAMMA2                                                        
	REAL*4  gamma                                                          
	REAL*4  E_zero_point                                                   
	REAL*4  Friction_factor                                                
                                                                        
	REAL*4  Delta_U1                                                       
	REAL*4  Delta_U2                                                       
	REAL*4  gamma_heavy1                                                   
	REAL*4  gamma_heavy2                                                   
	REAL*4  I_EVA                                                          
	REAL*4  E_saddle_scission                                              
	REAL*4  Ysymm                                                          
                                                                        
                                                                        
     	INTEGER*4  I_COUNT                                                
     	INTEGER*4  I_MODE,II                                                 

        COMMON/NEVENT/II
                                                                        
	Nheavy1 = 82   ! position of heavy peak valley 1                      
	Nheavy2 = 89    ! position of heavy peak valley 2                      
	E_crit = 5    !Critical pairing energy                                 
	Delta_U1_shell = -2.5  !Shell effect for valley 1                      
	Delta_U2_shell = -5.5  !Shell effect for valley 2                      
	CZ_asymm1_shell = 0.7  !Curvature of asymmetric valley 1               
	CZ_asymm2_shell = 0.08 !Curvature of asymmetric valley 2               
	Fwidth_asymm1 = 1.2      !Factor for width of distr. valley 1          
	Fwidth_asymm2 = 1.0      !Factor for width of distr. valley 2          
	CZ_asymm2_scission = 0.12!Curvature of asymmetric valley 2             
	FGAMMA1 = 1.0      !Factor to gamma_heavy1                             
	FGAMMA2 = 1.0      !Factor to gamma_heavy2                             
	gamma = 0        !fading of shells (general)                           
	El = 30        !fading of shells for energy deformation       
	E_zero_point = 0.5 !Zero-point energy at saddle                        
	Friction_factor = 1    !Friction factor                                
                                                                        
	Delta_U1 = 0           !I: used shell effect                           
	Delta_U2 = 0           !I: used shell effect                           
	gamma_heavy1 = 0 !I: fading of shell 1                                 
	gamma_heavy2 = 0 !I: fading of shell 2                                 
	I_EVA = 0        !Calculate A =1  or Aprime =0                         
	E_saddle_scission = 10. !I: friction from saddle to scision            
	Ysymm = 0 !I: Yield of symmetric mode                                  
                                                                        
                
     	I_COUNT = 0                                                       
     	I_MODE = 0                                                        
                                                                        
c average Z of asymmetric and symmetric components: */                  
     	N = A - Z  ! neutron number of the fissioning nucleus */          

        icz = 0                                            
        IF (Z**2/A.lt.25.or.N.lt.Nheavy2.or.E.gt.500) THEN
            icz = -1
            GOTO 1002
        END IF

     	Nlight1 = N - Nheavy1                                             
     	Nlight2 = N - Nheavy2                   

c UCD - 1/2 for A=cte ==> sig_z(N=cte) = 0.8!!!!!!                         
     	Zheavy1_shell = Nheavy1 / N * Z - 0.8E0 
     	Zheavy2_shell = Nheavy2 / N * Z - 0.8E0   
                                                                        
     	E_saddle_scission = (3.535E0 * Z**2/A - 121.1E0) * Friction_factor
c from Wilkins et al. Proc. Int. Symp. Nucl. Fission and                
C   Heavy Ion Induced Reactions, Schroeder W. ed.,Harwood, 1986, 101 */ 
                                                                        
        IF (E_saddle_scission.lt.0) THEN                                
     	   E_saddle_scission = 0                                          
        END IF                                                          
                                                                        
                                                                        
c Semiempirical fission model: */                                       
                                                                        
c Fit to experimental result on curvature of potential at saddle */     
C        reference:                                              */     
c     	IF (Z**2/A.LE.33.15E0) THEN                                       
c           MassCurv = 30.5438538E0 - 4.00212049E0 * Z**2/A              
c     #                          + 0.11983384E0 * Z**4 / (A**2)          
c     	ELSE                                                              
c           MassCurv = 10.E0 ** (7.16993332E0 - 0.26602401E0 * Z**2/A    
c     #                          + 0.00283802E0 * Z**4 / (A**2))         
c        END IF                                                          

c** New parametrization of T.Enqvist according to Mulgin et al NPA 640(1998)375
     	IF (Z**2/A.LT.33.9186E0) THEN                                       
           MassCurv = 10.E0 ** (-1.093364E0 + 0.082933E0 * Z**2/A    
     #                          - 0.0002602E0 * Z**4 / (A**2))         
     	ELSE                                                              
           MassCurv = 10.E0 ** (3.053536E0 - 0.056477E0 * Z**2/A    
     #                          + 0.0002454E0 * Z**4 / (A**2))         
        END IF                                                          
                    
     	CZ_symm = 8.E0 / Z**2 * MassCurv                                  
c        write(6,*)'A,Z,e,CZ',a,z,e,CZ_symm

c     	CZ_symm = -1.                                  
        icz = 0                                            
        IF (CZ_symm.lt.0) THEN
            icz = -1
            GOTO 1002
        END IF
                                            
     	Zsymm  = Z / 2.E0  ! proton number in symmetric fission (centre) *
     	Nsymm  = N / 2.E0                                                 
     	Asymm = Nsymm + Zsymm                                             
                                                                        
     	Zheavy1 = (CZ_symm*Zsymm + CZ_asymm1_shell*Zheavy1_shell)         
     #                      /(CZ_symm + CZ_asymm1_shell)                
     	Zheavy2 = (CZ_symm*Zsymm + CZ_asymm2_shell*Zheavy2_shell)         
     #                      /(CZ_symm + CZ_asymm2_shell)                
                                                                        
c position of valley due to influence of liquid-drop potential */       
     	Nheavy1_eff = (Zheavy1 + 0.8E0)*N/Z ! UCD - 1/2 unit assumed KHS *
     	Nheavy2_eff = (Zheavy2 + 0.8E0)*N/Z ! from Wagemanns p. 397     */
     	Nlight1_eff = N - Nheavy1_eff                                     
     	Nlight2_eff = N - Nheavy2_eff                                     
     	Zlight1 = Z - Zheavy1  ! proton number of light fragments (centre)
     	Zlight2 = Z - Zheavy2  ! proton number of light fragments (centre)
     	Aheavy1 = Nheavy1_eff + Zheavy1                                   
     	Aheavy2 = Nheavy2_eff + Zheavy2                                   
     	Aheavy1_mean = Aheavy1                                            
     	Aheavy2_mean = Aheavy2                                            
     	Alight1 = Nlight1_eff + Zlight1                                   
     	Alight2 = Nlight2_eff + Zlight2                                   
     	Eheavy1 = E * Aheavy1 / A                                         
     	Eheavy2 = E * Aheavy2 / A                                         
     	Elight1 = E * Alight1 / A                                         
     	Elight2 = E * Alight2 / A                                         
     	A_levdens = A / 8.E0                                              
     	A_levdens_heavy1 = Aheavy1 / 8.E0                                 
     	A_levdens_heavy2 = Aheavy2 / 8.E0                                 
     	A_levdens_light1 = Alight1 / 8.E0                                 
     	A_levdens_light2 = Alight2 / 8.E0                                 
     	gamma = A_levdens / (0.4E0 * A**1.3333E0)                         
     	gamma_heavy1 = A_levdens_heavy1 / (0.4E0 * Aheavy1**1.3333E0)     
     #                 *FGAMMA1                                         
     	gamma_heavy2 = A_levdens_heavy2 / (0.4E0 * Aheavy2**1.3333E0)     
     #                 *FGAMMA2                                         
                                                                        
     	CZ_asymm1_saddle = CZ_asymm1_shell + CZ_symm                      
     	CZ_asymm2_saddle = CZ_asymm2_shell + CZ_symm                      
                                                                        
     	CN =  MASS(Zsymm , Nsymm + 1.E0) + MASS(Zsymm, Nsymm - 1.E0)      
     #            + 1.44E0 * (Zsymm)**2 /                               
     #            (r_null**2 * ((Asymm+1)**(1./3.) +                    
     #            (Asymm-1)**(1./3.))**2 )                              
     #            - 2.E0 * MASS(Zsymm,Nsymm)                            
     #            - 1.44E0 * (Zsymm)**2 /                               
     #            (r_null * 2.E0 * (Asymm)**(1./3.))**2                 
c        WRITE(6,*)'MASS1',MASS(Zsymm , Nsymm + 1.E0)                   
c        WRITE(6,*)'MASS2',MASS(Zsymm , Nsymm - 1.E0)                   
c        WRITE(6,*)'MASS3',MASS(Zsymm , Nsymm)                          
                                                                        
                                                                        
                                                                        
c*** Minimum potencial value for the mass asymmetry                     
                                                                        
     	Delta_U1 = Delta_U1_shell + (Zheavy1_shell - Zheavy1)**2 *        
     #               CZ_asymm1_shell                                    
     	Delta_U2 = Delta_U2_shell + (Zheavy2_shell - Zheavy2)**2 *        
     #               CZ_asymm2_shell                                    
                                                                        
        Bsym = 0.                                                       
        Basym_1 = Bsym + (Zheavy1-Zsymm)**2*CZ_symm + Delta_U1          
        Basym_2 = Bsym + (Zheavy2-Zsymm)**2*CZ_symm + Delta_U2          
                                                                        
                                                                        
        IF (Bsym.lt.Basym_1.and.Bsym.lt.Basym_2) THEN                   
                                                                        
C***       EXCITATION ENERGIES AT THE SADDLE POINT      ****************
                                                                        
c***  excitation energy at saddle without and with shell effect         
     	   epsilon0_1_saddle = (E+E_zero_point - (Zheavy1-Zsymm)**2       
     #                           *CZ_symm)                              
     	   epsilon_1_saddle = epsilon0_1_saddle - Delta_U1                
                                                                        
c***  excitation energy at saddle without and with shell effect         
     	   epsilon0_2_saddle = (E+E_zero_point - (Zheavy2-Zsymm)**2       
     #                         *CZ_symm)                                
     	   epsilon_2_saddle = epsilon0_2_saddle - Delta_U2                
                                                                        
     	   epsilon_symm_saddle = E+E_zero_point                           
                                                                        
     	   EEXC1_saddle = epsilon_1_saddle                                
     	   EEXC2_saddle = epsilon_2_saddle                                
                                                                        
                                                                        
C***       EXCITATION ENERGIES AT THE SCISSION POINT      **************
                                                                        
c excitation energy of heavy fragment without and with shell effect     
     	epsilon0_1_scission = (E+E_saddle_scission - (Zheavy1-Zsymm)**2
     #                           * CZ_symm) * Aheavy1/A                 
     	epsilon_1_scission = epsilon0_1_scission - (Delta_U1)*Aheavy1/A
                                                                        
c excitation energy of heavy fragment without and with and shell effect 
     	epsilon0_2_scission = (E+E_saddle_scission - (Zheavy2-Zsymm)**2
     #                        * CZ_symm) * Aheavy2/A                    
c excitation energy of heavy fragment */                                
     	epsilon_2_scission = epsilon0_2_scission - (Delta_U2)*Aheavy2/A
                                                                        
     	   epsilon_symm_scission = E+E_saddle_scission                    
                                                                        
        ELSE IF (Basym_1.lt.Bsym.and.Basym_1.lt.Basym_2) THEN           
                                                                        
C***       EXCITATION ENERGIES AT THE SADDLE POINT      ****************
                                                                        
c***  excitation energy at saddle without and with shell effect         
     	   epsilon_symm_saddle = (E+E_zero_point + Delta_U1 +             
     #                           (Zheavy1-Zsymm)**2*CZ_symm)            
                                                                        
c***  excitation energy at saddle without and with shell effect         
     	epsilon0_2_saddle = (epsilon_symm_saddle - (Zheavy2-Zsymm)**2  
     #                         *CZ_symm)                                
     	   epsilon_2_saddle = epsilon0_2_saddle - Delta_U2                
                                                                        
     	   epsilon0_1_saddle = E+E_zero_point + Delta_U1                  
     	   epsilon_1_saddle = E+E_zero_point                              
                                                                        
                                                                        
     	   EEXC1_saddle = epsilon_1_saddle                                
     	   EEXC2_saddle = epsilon_2_saddle                                
                                                                        
                                                                        
C***       EXCITATION ENERGIES AT THE SCISSION POINT      **************
                                                                        
c excitation energy of heavy fragment without and with shell effect     
     	   epsilon_symm_scission = (E+E_saddle_scission +                 
     &     (Zheavy1-Zsymm)**2 * CZ_symm + Delta_U1)                     
                                                                        
c excitation energy of heavy fragment without and with and shell effect 
     	   epsilon0_2_scission = (epsilon_symm_scission -                 
     &     (Zheavy2-Zsymm)**2 * CZ_symm) * Aheavy2/A                    
     	   epsilon_2_scission = epsilon0_2_scission -                     
     &     (Delta_U2)*Aheavy2/A                                         
           epsilon0_1_scission=(E+E_saddle_scission + Delta_U1)         
     &     * Aheavy1/A                                                  
     	   epsilon_1_scission = (E+E_saddle_scission)*Aheavy1/A           
                                                                        
        ELSE IF (Basym_2.lt.Bsym.and.Basym_2.lt.Basym_1) THEN           
                                                                        
C***       EXCITATION ENERGIES AT THE SADDLE POINT      ****************
                                                                        
c***  excitation energy at saddle without and with shell effect         
     	   epsilon_symm_saddle = (E+E_zero_point + Delta_U2 +             
     #                           (Zheavy2-Zsymm)**2*CZ_symm)            
                                                                        
c***  excitation energy at saddle without and with shell effect         
     	   epsilon0_1_saddle = (epsilon_symm_saddle -                     
     &     (Zheavy1-Zsymm)**2 *CZ_symm)                                 
     	   epsilon_1_saddle = epsilon0_1_saddle - Delta_U1                
     	   epsilon0_2_saddle = E+E_zero_point + Delta_U2                  
     	   epsilon_2_saddle = E+E_zero_point                              
                                                                        
                                                                        
     	   EEXC1_saddle = epsilon_1_saddle                                
     	   EEXC2_saddle = epsilon_2_saddle                                
                                                                        
                                                                        
C***       EXCITATION ENERGIES AT THE SCISSION POINT      **************
                                                                        
c excitation energy of heavy fragment without and with shell effect     
     	   epsilon_symm_scission = (E+E_saddle_scission +                 
     &      (Zheavy2-Zsymm)**2 * CZ_symm + Delta_U2)                    
                                                                        
c excitation energy of heavy fragment without and with and shell effect 
     	   epsilon0_1_scission = (epsilon_symm_scission -                 
     &     (Zheavy1-Zsymm)**2 * CZ_symm) * Aheavy1/A                    
     	   epsilon_1_scission = epsilon0_1_scission - (Delta_U1)          
     &     * Aheavy1/A                                                  
           epsilon0_2_scission=(E+E_saddle_scission + Delta_U2)
     &     * Aheavy2/A
     	   epsilon_2_scission = (E+E_saddle_scission)*Aheavy2/A           
                                                                        
                                                                        
	END IF                                                                 
                                                                        
        IF (epsilon_1_saddle.lt.0)epsilon_1_saddle = 0                  
        IF (epsilon0_1_saddle.lt.0)epsilon0_1_saddle = 0                
        IF (epsilon_2_saddle.lt.0)epsilon_2_saddle = 0                  
        IF (epsilon0_2_saddle.lt.0)epsilon0_2_saddle = 0                
        IF (epsilon_symm_saddle.lt.0)epsilon_symm_saddle = 0            
                                                                        
        IF (epsilon_1_scission.lt.0)epsilon_1_scission = 0              
        IF (epsilon0_1_scission.lt.0)epsilon0_1_scission = 0            
        IF (epsilon_2_scission.lt.0)epsilon_2_scission = 0              
        IF (epsilon0_2_scission.lt.0)epsilon0_2_scission = 0            
        IF (epsilon_symm_scission.lt.0)epsilon_symm_scission = 0        
                                                                        
c        write(6,*)'E,E1,E2,Es',e,epsilon_1_saddle,epsilon_2_saddle,    
c     #             epsilon_symm_saddle                                 
                                                                        
                                                                        
c Calculate widhts at the saddle:                */                     
                                                                        
     	E_eff1_saddle = epsilon0_1_saddle -                               
     #              Delta_U1 * exp(-epsilon_1_saddle*gamma)             
     	IF (E_eff1_saddle .gt. 0) THEN                                    
           WZasymm1_saddle = sqrt(0.5E0 *                               
     #     sqrt(1.E0/A_levdens*E_eff1_saddle) /                          
     #     (CZ_asymm1_shell * exp(-epsilon_1_saddle*gamma) + CZ_symm)) 
     	ELSE                                                              
  	   WZasymm1_saddle = 1                                               
        END IF                                                          
                                                                        
     	E_eff2_saddle = epsilon0_2_saddle -                               
     #              Delta_U2 * exp(-epsilon_2_saddle*gamma)             
     	IF (E_eff2_saddle.GT.0) THEN  
           WZasymm2_saddle = sqrt(0.5E0 *                               
     #     sqrt(1.E0/A_levdens*E_eff2_saddle) /                          
     #    (CZ_asymm2_shell * exp(-epsilon_2_saddle*gamma) + CZ_symm))  
     	ELSE                                                              
	   WZasymm2_saddle = 1                                                 
        END IF                                                          
                                                                        
     	IF ((E+E_zero_point).GT.0) THEN                                   
     	   WZsymm_saddle = sqrt(0.5E0 *                                   
     #     sqrt(1.E0/A_levdens*(E+epsilon_symm_saddle)) / CZ_symm)      
     	ELSE                                                              
	   WZsymm_saddle = 1                                                   
        END IF                                                          
                                                                        
c Calculate widhts at the scission point: */                            
c fits of ref. Beizin 1991 (Plots brought to GSI by Sergei Zhdanov) */  
                                                                        
     	WZsymm_scission = WZsymm_saddle                                   
                                                                        
     	IF (E_saddle_scission.EQ.0) THEN                                  
           WZasymm1_scission = WZasymm1_saddle                          
           WZasymm2_scission = WZasymm2_saddle                          
     	ELSE                                                              
           IF (Nheavy1_eff.GT.75) THEN                                  
              WZasymm1_scission = sqrt(21.E0)*Z/A                       
              RR = (70.E0-28.E0)/3.E0 * (Z**2/A - 35.E0) + 28.E0        
              IF (RR.GT.0) THEN                                         
                 WZasymm2_scission = sqrt(RR)*Z/A                       
              ELSE                                                      
                 WZasymm2_scission = 0.                                 
              END IF                                                    
C              WZasymm2_scission = sqrt(MAX( (70.E0-28.E0)/3.E0 *       
C     #                         (Z**2/A - 35.E0) + 28.E0 , 0))*Z/A      
           ELSE                                                         
              WZasymm1_scission = WZasymm1_saddle                       
              WZasymm2_scission = WZasymm2_saddle                       
           END IF                                                       
        END IF                                                          
                                                                        
     	WZasymm1_scission = MAX(WZasymm1_scission,WZasymm1_saddle)        
     	WZasymm2_scission = MAX(WZasymm2_scission,WZasymm2_saddle)        
                                                                        
                                                                        
     	WZasymm1 = WZasymm1_scission * Fwidth_asymm1                      
     	WZasymm2 = WZasymm2_scission * Fwidth_asymm2                      
     	WZsymm = WZsymm_scission                                          
                                                                        
     	WAsymm = WZsymm * A/Z                                             
     	WAheavy1 = WZasymm1 * A/Z                                         
     	WAheavy2 = WZasymm2 * A/Z                                         
                                                                        
     	WAsymm_saddle = WZsymm_saddle * A/Z                               
     	WAheavy1_saddle = WZasymm1_saddle * A/Z                           
     	WAheavy2_saddle = WZasymm2_saddle * A/Z                           
                                                                        
                                                                        
c        write(6,*)'al,e,es,cn',A_levdens,E,E_saddle_scission,cn        

c*** sig_0 = quantum fluctuation = 0.45 z units for A=cte
c                                = 0.45*2.58 = 1.16 n units for Z=cte
c                       sig_0**2  = 1.16*2 = 1.35    n units for Z=cte

     	N1width = sqrt(0.5E0 *                                            
     #         sqrt(1.E0/A_levdens*(E+E_saddle_scission)) / CN + 1.35)         
                                                                        
     	IF (epsilon0_1_saddle - Delta_U1 *                                
     #               EXP(-epsilon_1_saddle*gamma_heavy1).LT.0)  THEN    
     	   Sasymm1 = -10                                                  
     	ELSE                                                              
     	   Sasymm1 = 2.E0 * SQRT(A_levdens*(epsilon0_1_saddle -           
     #	   Delta_U1 * EXP(-epsilon_1_saddle*gamma_heavy1)))              
	END IF                                                                 
                                                                        
     	IF (epsilon0_2_saddle - Delta_U2 *                                
     #               EXP(-epsilon_2_saddle*gamma_heavy2).LT.0)  THEN    
     	   Sasymm2 = -10                                                  
     	ELSE                                                              
     	Sasymm2 = 2.E0 * SQRT(A_levdens*(epsilon0_2_saddle - Delta_U2 *
     #               EXP(-epsilon_2_saddle*gamma_heavy2)))              
      	END IF                                                           
                                                                        
     	IF (epsilon_symm_saddle.GT.0) THEN                                
     	   Ssymm = 2.E0 * SQRT(A_levdens*(epsilon_symm_saddle))           
     	ELSE                                                              
	   Ssymm = -10                                                         
        END IF                                                          
                                                                        
     	IF (Ssymm.GT.-10) THEN                                            
       	   Ysymm = 1.E0                                                 
                                                                        
           IF (epsilon0_1_saddle.LT.0) THEN                             
       	      Yasymm1 = EXP(Sasymm1-Ssymm)                              
     #             * WZasymm1_saddle / WZsymm_saddle * 2.E0  ! low energ
c factor of 2 for symmetry classes */                                   
           ELSE                          ! high energy */               
              Ssymm_mode1 = 2.E0 * sqrt(A_levdens*(epsilon0_1_saddle))  
              Yasymm1 = (EXP(Sasymm1-Ssymm) - EXP(Ssymm_mode1 - Ssymm)) 
     #             * WZasymm1_saddle / WZsymm_saddle * 2.E0             
           END IF                                                       
                                                                        
           IF (epsilon0_2_saddle.LT.0) THEN                             
              Yasymm2 = EXP(Sasymm2-Ssymm)                              
     #             * WZasymm2_saddle / WZsymm_saddle * 2.E0  ! low energ
c factor of 2 for symmetry classes */                                   
           ELSE                         ! high energy */                
              Ssymm_mode2 = 2.E0 * sqrt(A_levdens*(epsilon0_2_saddle))  
              Yasymm2 = (EXP(Sasymm2-Ssymm) - EXP(Ssymm_mode2 - Ssymm)) 
     #             * WZasymm2_saddle / WZsymm_saddle * 2.E0             
           END IF                                                       
c difference in the exponent in order */                                
c to avoid numerical overflow         */                                
     	ELSE                                                              
           IF ((Sasymm1.GT.-10).AND.(Sasymm2.GT.-10)) THEN              
              Ysymm = 0                                                 
              Yasymm1 = EXP(Sasymm1) * WZasymm1_saddle * 2.E0           
              Yasymm2 = EXP(Sasymm2) * WZasymm2_saddle * 2.E0           
           END IF                                                       
        END IF                                                          
                                                                        
     	Ysum = Ysymm + Yasymm1 + Yasymm2  ! normalize */                  
     	IF (Ysum.GT.0) THEN                                               
           Ysymm = Ysymm / Ysum                                         
           Yasymm1 = Yasymm1 / Ysum                                     
           Yasymm2 = Yasymm2 / Ysum                                     
           Yasymm = Yasymm1 + Yasymm2                                   
     	ELSE                                                              
           Ysymm = 0                                                    
           Yasymm1 = 0                                                  
           Yasymm2 = 0                                                  
c search minimum threshold and attribute all events to this mode */     
           IF ((epsilon_symm_saddle.LT.epsilon_1_saddle).AND.           
     #        (epsilon_symm_saddle.LT.epsilon_2_saddle)) THEN           
              Ysymm = 1                                                 
           ELSE                                                         
              IF (epsilon_1_saddle.LT.epsilon_2_saddle) THEN            
                 Yasymm1 = 1                                            
              ELSE                                                      
                 Yasymm2 = 1                                            
              END IF                                                    
           END IF                                                       
        END IF                                                          
                                                                        
                                                                        
                                                                        
c even-odd effect */                                                    
c simple parametrization KHS */                                         
     	EEXC_MAX = MAX(Eexc1_saddle,Eexc2_saddle)                         
     	EEXC_MAX = MAX(EEXC_MAX,E)                                        
        IZ = NINT(Z)                                                    
c        write(6,*)'mod(z,2)',jmod(iz,2)                                
     	IF (MOD(IZ,2).EQ.0) THEN                                         
            R_E_O_MAX = 0.3E0 * (1.E0 -  0.2E0 * ( Z**2/A -             
     #              92.E0**2.E0/238.E0 ))                               
           E_PAIR = 2.E0 * 12.E0 / SQRT(A)                              
           IF (EEXC_MAX.GT.(E_CRIT + E_PAIR)) THEN                      
              R_E_O = 0                                                 
           ELSE                                                         
              IF (EEXC_MAX.LT.E_PAIR) THEN                              
                 R_E_O = R_E_O_MAX                                      
              ELSE                                                      
                 R_E_O =                                                
     #            ((EEXC_MAX-E_CRIT-E_PAIR)/(E_CRIT))**2 * R_E_O_MAX    
              END IF                                                    
           END IF                                                       
        ELSE                                                            
           R_E_O = 0                                                    
        END IF                                                          
c        write(6,*)'rmax',R_E_O_MAX                                     
                                                                        
c        write(6,*)'ex_max,ec,ep',eexc_max,e_crit,e_pair                
c        IF (r_e_o.gt.0) write(6,*)'e_crit,r_e_o',e_crit,r_e_o          
                                                                        
c  R_E_O = EXP( -( E + E_saddle_scission) / 4.E0) */                    
c  Nifenecker model */                                                  
                                                                        
                                                                        
                                                                        
                                                                        
c random decision: symmetric or asymmetric */                           
c   IMODE = 1 means symmetric fission,                                  
C   IMODE = 2 means asymmetric fission, mode 1,                         
C   IMODE = 3 means asymmetric fission, mode 2 */                       
     	RMODE = haz(kkk)                                                  
       	IF (RMODE.LT.Ysymm) THEN                                        
            IMODE = 1                                                   
       	ELSE IF (RMODE.LT.(Ysymm + Yasymm1)) THEN                       
            IMODE = 2                                                   
       	ELSE                                                            
            IMODE = 3                                                   
        END IF                                                          
                                                                        
c determine parameters of the Z distribution */                         
       	IF (IMODE.EQ.1) THEN                                            
           Z1mean = Zsymm                                               
           Z1width = WZsymm                                             
        ELSE IF (IMODE.EQ.2) THEN                                       
           Z1mean = Zheavy1                                             
           Z1width = WZasymm1                                           
        ELSE IF (IMODE.EQ.3) THEN                                       
           Z1mean = Zheavy2                                             
           Z1width = WZasymm2                                           
        END IF                                                          
                                                                        
c random decision: Z1 and Z2 at scission: */                            
        Z1 = 1.0
        Z2 = 1.0
        DO WHILE (Z1.LT.5.OR.Z2.LT.5)
           Z1 = GAUSSHAZ(KKK,Z1mean,Z1width)                               
     	   CALL EVEN_ODD(Z1,R_E_O,I_HELP)                                    
     	   Z1 = FLOAT(I_HELP)                                                
     	   Z2 = Z - Z1                                                       
        END DO
                                                     
c average N of both fragments: */                                       
       	IF (IMODE.EQ.1) THEN  
           N1ucd = Z1 * N/Z                               
           N2ucd = Z2 * N/Z
           re1 = massdef(Z1,N1ucd,0.6) + massdef(Z2,N2ucd,0.6) +
     1           ecoul(Z1,N1ucd,0.6,Z2,N2ucd,0.6,2.0)
           re2 = massdef(Z1,N1ucd+1.0,0.6) + massdef(Z2,N2ucd-1.0,0.6) +
     1           ecoul(Z1,N1ucd+1.0,0.6,Z2,N2ucd-1.0,0.6,2.0)
           re3 = massdef(Z1,N1ucd+2.0,0.6) + massdef(Z2,N2ucd-2.0,0.6) +
     1           ecoul(Z1,N1ucd+2.0,0.6,Z2,N2ucd-2.0,0.6,2.0)
           reps2 = (re1-2.0*re2+re3)/2.0
           reps1 = re2 - re1 -reps2
           rn1_pol = -reps1/(2.0*reps2)
           N1mean = N1ucd + rn1_pol                                            
       	ELSE                                                            
           N1mean = (Z1 + 0.5E0) * N/Z                                  
       	END IF                                                          
                                                                        
c N1mean = Nsymm + (Z1 - Zsymm) * 1.6E0  1.6 from 238U(nth,f) */        
c N1width = 0.9 + E * 0.002E0  KHS */                                   
                                                                        
c random decision: N1R and N2R at scission, before evaporation: */      
        N1R = 1.0
        N2R = 1.0
        DO WHILE (N1R.LT.5.OR.N2R.LT.5)
           N1R = GAUSSHAZ(KKK,N1mean,N1width)                              
C modification to have N1R as integer, and N=N1R+N2R rigorously A.B. 19/4/2001
	   I_INTER = N1R + 0.5
	   N1R = FLOAT(I_INTER)
	   
           N2R = N - N1R                                                     
        END DO                                                                
                                                                        
c neutron evaporation from fragments */                                 
     	IF (I_EVA.GT.0) THEN                                              
c Treatment SZ */                                                       
     	   NE_MIN = 0.095E0 * A - 20.4E0                                  
           IF (NE_MIN.LT.0) NE_MIN = 0.                                 
     	   NE_MIN = NE_MIN + E / 8.E0  ! 1 neutron per 8 MeV */           
     	   A1R = Z1 + N1R              ! Mass of first fragment */        
     	   NE_M1 = A1R / A * NE_MIN    ! Devide nbr. of neutrons acc. mass
     	   NE_M2 = NE_MIN - NE_M1      ! nmbr. of neutrons of 2. fragment 
     	   N1 = N1R - NE_M1            ! final neutron number 1. fragment 
     	   N2 = N2R - NE_M2            ! final neutron number 2. fragment 
     	ELSE                                                              
     	   N1 = N1R                                                       
     	   N2 = N2R                                                       
     	END IF  



c Excitation energy due to deformation                                 

     	A1 = Z1 + N1R              ! Mass of first fragment */        
     	A2 = Z2 + N2R              ! Mass of second fragment */        



        IF (A1.lt.80) THEN
     	   ED1_LOW = 0.
        ELSE IF (A1.ge.80.and.A1.lt.110) THEN
     	   ED1_LOW = (A1-80.)*20./30.
        ELSE IF (A1.ge.110.and.A1.lt.130) THEN
     	   ED1_LOW = -(A1-110.)*20./20. + 20.
        ELSE IF (A1.ge.130) THEN
     	   ED1_LOW = (A1-130.)*20./30.
        END IF

        IF (A2.lt.80) THEN
     	   ED2_LOW = 0.
        ELSE IF (A2.ge.80.and.A2.lt.110) THEN
     	   ED2_LOW = (A2-80.)*20./30.
        ELSE IF (A2.ge.110.and.A2.lt.130) THEN
     	   ED2_LOW = -(A2-110.)*20./20. + 20.
        ELSE IF (A2.ge.130) THEN
     	   ED2_LOW = (A2-130.)*20./30.
        END IF


        Ed1_high = 20.0*A1/(A1+A2)                                          
        Ed2_high = 20.0 - Ed1_high                                          

        Ed1 = Ed1_low*exp(-E/El) + Ed1_high*(1-exp(-E/El))
        Ed2 = Ed2_low*exp(-E/El) + Ed2_high*(1-exp(-E/El))


c        write(6,101)E,A1,A2,Ed1,Ed2,Ed1+Ed2
c        write(6,102)Ed1_low,Ed1_high,Ed2_low,Ed2_high

 
        E1 = E*A1/(A1+A2) + ED1
        E2 = E - E*A1/(A1+A2) + ED2

        Atot = A1+A2
        IF (Atot.gt.A+1) THEN
           write(6,*)'A,,A1,A2,Atot',A,A1,A2,Atot
C           write(6,*)'N,N1R,N2R',N,N1R,N2R
C           write(6,*)'Z,Z1,Z2',Z,Z1,Z2
        END IF


1002    CONTINUE
c**************************
c*** only symmetric fission
c**************************
        IF ((icz.eq.-1).or.(A1.lt.0.).or.(A2.lt.0.)) THEN
C           IF (z.eq.92) THEN
C              write(6,*)'symmetric fission'
C              write(6,*)'Z,A,E,A1,A2,icz,Atot',Z,A,E,A1,A2,icz,Atot
C           END IF


        N = A-Z
     	Zsymm  = Z / 2.E0  ! proton number in symmetric fission (centre) *
     	Nsymm  = N / 2.E0                                                 
     	Asymm = Nsymm + Zsymm                                             

     	A_levdens = A / 8.E0                                              

        MassCurv = 2
     	CZ_symm = 8.E0 / Z**2 * MassCurv                                  

     	WZsymm = sqrt(0.5E0 * sqrt(1.E0/A_levdens*E) / CZ_symm)      

        Z1mean = Zsymm                                              
        Z1width = WZsymm                     

c random decision: Z1 and Z2 at scission: */
        Z1 = 1.0                 
        Z2 = 1.0                 
        DO WHILE (Z1.lt.5.or.z2.lt.5)
           Z1 = GAUSSHAZ(KKK,Z1mean,Z1width)                               
     	   Z2 = Z - Z1
        END DO                                            
                                                      

c        write(6,*)'Zsymm,Nsymm,Asymm',Zsymm,Nsymm,Asymm

     	CN =  MASS(Zsymm , Nsymm + 1.E0) + MASS(Zsymm, Nsymm - 1.E0)      
     #            + 1.44E0 * (Zsymm)**2 /                               
     #            (r_null**2 * ((Asymm+1)**(1./3.) +                    
     #            (Asymm-1)**(1./3.))**2 )                              
     #            - 2.E0 * MASS(Zsymm,Nsymm)                            
     #            - 1.44E0 * (Zsymm)**2 /                               
     #            (r_null * 2.E0 * (Asymm)**(1./3.))**2                 


c*** sig_0 = quantum fluctuation = 0.45 z units for A=cte
c                                = 0.45*2.58 = 1.16 n units for Z=cte
c                       sig_0**2  = 1.16**2 = 1.35    n units for Z=cte

     	N1width = sqrt(0.5E0 *                                            
     #         sqrt(1.E0/A_levdens*E) / CN + 1.35)         


        N1ucd = Z1 * N/Z                               
        N2ucd = Z2 * N/Z
        re1 = massdef(Z1,N1ucd,0.6) + massdef(Z2,N2ucd,0.6) +
     1        ecoul(Z1,N1ucd,0.6,Z2,N2ucd,0.6,2.0)
        re2 = massdef(Z1,N1ucd+1.0,0.6) + massdef(Z2,N2ucd-1.0,0.6) +
     1        ecoul(Z1,N1ucd+1.0,0.6,Z2,N2ucd-1.0,0.6,2.0)
        re3 = massdef(Z1,N1ucd+2.0,0.6) + massdef(Z2,N2ucd-2.0,0.6) +
     1        ecoul(Z1,N1ucd+2.0,0.6,Z2,N2ucd-2.0,0.6,2.0)
        reps2 = (re1-2.0*re2+re3)/2.0
        reps1 = re2 - re1 -reps2
        rn1_pol = -reps1/(2.0*reps2)
        N1mean = N1ucd + rn1_pol                                            


c random decision: N1R and N2R at scission, before evaporation: */      
        N1R = NINT(GAUSSHAZ(KKK,N1mean,N1width))                              
     	N2R = N - N1R                                                     

     	A1 = Z1 + N1R              ! Mass of first fragment */        
     	A2 = Z2 + N2R              ! Mass of second fragment */        

        E1 = E*A1/(A1+A2)                                            
        E2 = E - E*A1/(A1+A2) 

        END IF
                                                         
c101     FORMAT(' E,A1,A2,Ed1,Ed2,S',6F8.2)
c102     FORMAT(' Ed1_l,Ed1_h,Ed2_l,Ed2_h',4F8.2)                            
                                                                       
        RETURN                                                          
	END        


                                                                        
C=======================================================================
                                                                        
   	SUBROUTINE EVEN_ODD(R_IN,R_EVEN_ODD,I_OUT)                          
c Procedure to calculate I_OUT from R_IN in a way that        */        
c on the average a flat distribution in R_IN results in a     */        
c flucuating distribution in I_OUT with an even-odd effect as */        
c given by R_EVEN_ODD                                         */        
c EXAMPLES :                                                  */        
c    If R_EVEN_ODD = 0 :                                      */        
c           CEIL(R_IN)  ----                                  */        
c              R_IN ->                                        */        
c                                                             */        
c           FLOOR(R_IN) ----       --> I_OUT                  */        
c    If R_EVEN_ODD > 0 :                                      */        
c      The interval for the above treatment is                */        
c         larger for FLOOR(R_IN) = even and                  */         
c         smaller for FLOOR(R_IN) = odd                        */       
c    For R_EVEN_ODD < 0 : just opposite treatment             */        
                                                                        
	REAL*4     R_IN,R_EVEN_ODD,R_REST,R_HELP                               
       	REAL*4     R_FLOOR                                              
       	REAL*4     R_MIDDLE                                             
       	INTEGER*4  I_OUT,I_IN,N_FLOOR                                   
                                                                        
                                                                        
        I_IN = INT(R_IN)                                                
	R_FLOOR = INT(R_IN)                                                    
                                                                        
c        write(6,*)'r_e_0',R_EVEN_ODD                                   
     	IF (R_EVEN_ODD.eq.0) THEN                                         
           I_OUT = INT(R_FLOOR)                                         
        ELSE                                                            
           R_REST = R_IN - R_FLOOR                                      
           R_MIDDLE = R_FLOOR + 0.5E0                                   
                                                                        
           N_FLOOR = INT(R_FLOOR)                                       
           IF (MOD(N_FLOOR,2).eq.0) THEN  ! even before modif. */      
              R_HELP = R_MIDDLE + (R_REST - 0.5E0) * (1.E0 - R_EVEN_ODD)
           ELSE   ! c odd before modification */                        
              R_HELP = R_MIDDLE + (R_REST - 0.5E0) * (1.E0 + R_EVEN_ODD)
           END IF                                                       
           I_OUT = INT(R_HELP)                                          
        END IF                                                          
                                                                        
                                                                        
        RETURN                                                          
        END	                                                            
                                                                        
C=======================================================================
                                                                        
   	REAL FUNCTION MASS(Z,N)                                                  
                                                                        
c liquid-drop mass, Myers & Swiatecki, Lysekil, 1967  */                
c pure liquid drop, without pairing and shell effects */                
                                                                        
     	REAL*4   Z,N,A                                                    
      	REAL*4   XCOM,XVS,XE,EL                                     
                                                                        
     	A = N + Z                                                         
     	XCOM = 1.E0 - 1.7826 * ((A - 2.E0*Z)/A)**2                        
     	XVS = - XCOM * (15.4941E0*A - 17.9439E0*A**(2.E0/3.E0))           
     	XE = Z**2 * (0.7053E0/A**(1.E0/3.E0) - 1.1529E0/A)                
     	EL = XVS + XE                                                     
     	MASS = XVS + XE                                                   
c        WRITE(6,*)'Z,N,XCOM,XVS,XE',Z,N,XCOM,XVS,XE                    
c        WRITE(6,*)'Z,N,MASS',Z,N,MASS                                  
                                                                        
     	RETURN                                                            
   	END                                                                 

C=======================================================================
                                                                        
   	FUNCTION MASSDEF(Z,N,beta)                                                  
                                                                        
c liquid-drop mass, Myers & Swiatecki, Lysekil, 1967  */                
c pure liquid drop, without pairing and shell effects */                
                                                                        
     	REAL*4   Z,N,A,beta,pi                                                    
      	REAL*4   XCOM,XVS,XE,EL,MASSDEF                                     
                              
c        write(6,*)'N,Z,beta',N,Z,beta
        pi = 3.1416                                 
     	A = N + Z                      
c        write(6,*)'A,Z,beta',A,Z,beta
        alpha = sqrt(5.0/(4.0*pi)) * beta                           
     	XCOM = 1.E0 - 1.7826 * ((A - 2.E0*Z)/A)**2  
c        write(6,*)'alpha,xcom',alpha,xcom
c** factor for asymmetry dependence of surface and volume term                      
     	XVS = - XCOM * (15.4941E0*A - 17.9439E0*A**(2.E0/3.E0) *
     1           ( 1.0 + 0.4*alpha**2))
c** sum of volume and suface energy           
     	XE = Z**2 * (0.7053E0/A**(1.E0/3.E0)*(1.0-0.2*alpha**2) - 
     1       1.1529E0/A)                

     	EL = XVS + XE                    
c        write(6,*)'xvs,xe,el',xvs,xe,el                         
     	MASSDEF = XVS + XE                                                   
c        WRITE(6,*)'Z,N,XCOM,XVS,XE',Z,N,XCOM,XVS,XE                    
c        WRITE(6,*)'Z,N,MASS',Z,N,MASSDEF                                  
                                                                        
     	RETURN                                                            
   	END                                                                 

C=======================================================================
                                                                        
   	FUNCTION ECOUL(Z1,N1,beta1,Z2,N2,beta2,d)                                                  
                                                                        
c Coulomb energy between two deformed nuclei (Wilkins 1976)                
                                                                        
     	REAL*4   Z1,N1,A1,beta1,Z2,N2,A2,beta2,pi                                                    
      	REAL*4   dtot,r0,ECOUL,d

c        write(6,*)'Z1,N1,beta1',Z1,N1,beta1
c        write(6,*)'Z2,N2,beta2,d',Z2,N2,beta2,d
        r0 = 1.16
        A1 = Z1+N1
        A2 = Z2+N2
c        write(6,*)'A1,A2',a1,a2
        dtot = r0*(A1**(1./3.) * (1.0+0.6666*beta1) +
     1             A2**(1./3.) * (1.0+0.6666*beta2) ) + d

        ECOUL = 1.44*Z1*Z2/dtot     
c        write(6,*)'dtot,ecoul',dtot,ecoul

     	RETURN                                                            
   	END                                                                 
                                                                        
C=======================================================================
                                                                        
       FUNCTION HAZ(K)                                                  
       REAL*4 P(110)                                                    
       INTEGER*4 IX,K,I                                                 
       REAL*4 X,Y,A,HAZ                                                 
C***      SI K=<-1 ON INITIALISE                                        
C***      SI K=-1 C'EST REPRODUCTIBLE                                   
C***      SI K<>-1 CE N'EST PAS REPRODUCTIBLE

       SAVE                           
                                                                        
          IF (K.LE.-1) THEN                                             
             IF(K.EQ.-1)THEN                                            
                              IX=0                                      
                        ELSE                                            
                              X=0.                                      
                              Y=SECNDS(X)                               
                              IX=Y*100+43543000                         
                              IF(MOD(IX,2).EQ.0)IX=IX+1                 
             END IF                                                     
             X=RANF(IX)                                                  
             DO I=1,110                                                 
                P(I)=RANF(IX)                                            
             END DO                                                     
             A=RANF(IX)                                                  
             K=0                                                        
          END IF                                                        
          I=NINT(100*A)+1                                               
          HAZ=P(I)                                                      
          A=RANF(IX)                                                     
          P(I)=A                                                        
          RETURN                                                        
        END                                                             
                                                                        
C********************************************************************** 
                                                                        
        FUNCTION GAUSSHAZ(K,XMOY,SIG)                                   
C*** TIRAGE ALEATOIRE DANS UNE GAUSSIENNE DE LARGEUR SIG ET MOYENNE XMOY
        INTEGER*4  ISET,K                                               
        REAL*4     V1,V2,HAZ,R,FAC,GSET,GAUSSHAZ,SIG,XMOY               
        
	SAVE
	                                                                
        DATA ISET/0/                                                    
        IF(ISET.EQ.0) THEN                                              
1          V1=2.*HAZ(K)-1.                                              
           V2=2.*HAZ(K)-1.                                              
           R=V1**2+V2**2                                                
           IF(R.GE.1)GO TO 1                                            
           FAC=SQRT(-2.*LOG(R)/R)                                       
           GSET=V1*FAC                                                  
           GAUSSHAZ=V2*FAC*SIG+XMOY                                     
           ISET=1                                                       
        ELSE                                                            
           GAUSSHAZ=GSET*SIG+XMOY                                       
           ISET=0                                                       
        END IF                                                          
        RETURN                                                          
        END                                                             
                                                                        
                                                                        
                 
