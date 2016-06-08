#. HCAL 96 Test beam description file for HCAL part
*DO HCal
#. Mother volume of HCal
#. Material     dy/2    dx/2	Xpos
   "Air"        375.0 1013.0	361.0
#. Boxes
#. Material	Number	dy/2	dx/2	Wall_T	Positions
   "Aluminium"	2	375.0	490.0	4.0	-523.0	509.0
#. Scintillator Layers
#. Type_Name		Number
#.	Type	Mother	Xpos	Type	Mother	Xpos
   "ScintillatorLayer"	28
	2	0	-427.0	1	0	-393.5
	0	0	-351.5	0	0	-314.5
	0	0	-277.5	0	0	-238.5
	0	0	-201.5	0	0	-164.5
	0	0	 -95.5	0	0	 -27.5
	0	0	  39.5	0	0	 107.5
	0	0	 175.5	0	0	 243.5
	0	0	 318.5	0	0	 388.5
	0	-1	 -14.5	0	1	-412.5
	0	1	-344.5	0	1	-274.5
	0	1	-205.5	0	1	-136.5
	2	1	 -47.0	1	1	  45.5
	2	1	 138.0	1	1	 232.5
	2	1	 328.0	0	-1	1009.5
#. Absorber Layers
#. Type_Name		Number
#.	Type	Mother	Xpos	Type	Mother	Xpos
   "AbsorberLayer"	28
	0	0	-410.0	1	0	-370.0
	1	0	-333.0	1	0	-296.0
	1	0	-257.0	1	0	-220.0
	1	0	-183.0	3	0	-129.0
	3	0	 -61.0	3	0	   6.0
	3	0	  74.0	3	0	 142.0
	3	0	 210.0	3	0	 285.0
	3	0	 355.0	3	0	 429.0
	3	1	-446.0	3	1	-378.0
	3	1	-308.0	3	1	-239.0
	3	1	-170.0	4	1	 -92.0
	4	1	  -1.0	4	1	  93.0
	4	1	 186.0	4	1	 283.0
	2	1	 358.0	2	1	 438.0
#. Absorber 
#. Material	Number	dY/2	Half thicknesses
   "Copper"	5	330.0	10.0	15.0	20.0	30.0	40.0
#. Scintillator
#. Material		Wrapper	 	Plastic		Number
#.	dY/2	dX/2(L)	dX/2(W)	dX/2(FP)	dX/2(BP)	dX/2(Sc)
   "Scintillator"	"Aluminium"	"Polystyrene"	3
	320.0	3.5	0.0	0.5		1.0		2.0
	320.0	6.5	0.0	0.5		1.0		5.0
	330.0	5.0	0.5	0.5		0.5		2.0
*ENDDO
