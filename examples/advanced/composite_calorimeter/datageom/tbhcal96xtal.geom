#. HCAL 96 Test beam description file for crystal matrix
*DO CrystalMatrix
#. Envelope
#. Material     Width	Length	Position in space and Rotation angles
   "Air"        320.0	682.0
	-1033.0	0.0	0.0	90.0	90.0	0.0	0.0	90.0	0.0
#.	-1033.0	0.0	0.0	90.0	270.0	0.0	0.0	90.0	180.0
#. Layer
#. Material	Number	Radius	Angle	FrontLength	Parameters
   "Air"	7	1500.0	0.01386	100.0
				10.5	11.9	73.5	85.0	115.5
#. Crystal
#. Material	Number	Length	Toler.	Parameters
   "E_PbWO4"	7	232.0	0.5	10.25	11.90	10.25	11.90	115.0
#. Support
#. Material	dX	dY	dZ	dist
   "Copper"	170.0	5.5	37.0	20.0	
*ENDDO
