
      COMMON/CELOSS/ EdepAbs(MaxAbs),TrckAbs(MaxAbs),
     &               sumEdep(MaxAbs),sumTrck(MaxAbs),
     &               su2Edep(MaxAbs),su2Trck(MaxAbs)
     
      PARAMETER (MaxPlanes = MaxAbs*MaxLay + 1)
      COMMON/EFLOW/ EnerFlow(MaxPlanes),EleakLat(MaxPlanes)
