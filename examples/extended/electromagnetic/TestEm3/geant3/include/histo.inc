      CHARACTER*25 fileName
      COMMON/HISTO1/fileName
      
      COMMON/HISTO2/idhist,nbBins,valmin,valmax,valunit
      
      PARAMETER (MaxHist = 22)
      LOGICAL histo            
      COMMON/HISTO3/histo(MaxHist),histUnit(MaxHist),binWidth(MaxHist) 
