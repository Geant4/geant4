
      PARAMETER (MaxAbs = 10)                 
      COMMON/CALOR/materAbs(MaxAbs),thickAbs(MaxAbs),NbAbsor,NbLayer,
     &             calorYZ,Field,thLayer,calorX,worldX,worldYZ
 
*                                                              
*  NbAbsor = Number of Absorbers                    (data card CALOR)
*  NbLayer = Number of Layers                       (data card CALOR) 
*  calorYZ = tranverse size of calor (cm)           (data card CALOR)
*  Field   = magnetic field (tesla)                 (data card CALOR)
*  materAbs(i) = Absorber material number.          (data card MATE )
*  thickAbs(i) = thickness of absorber   (cm)       (data card THICK)
*
*  thLayer = total thickness of a layer (cm)
*  calorX  = total  X-size of the calor (cm)
*  worldX  = total  X-size of the world (cm)
*  worldYZ = total YZ-size of the world (cm)
*                                                                               
