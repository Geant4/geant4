      REAL FUNCTION slope90(x)
      REAL x,y
      COMMON/PAWPAR/PAR(4)
      slope90=PAR(1)/(1.-PAR(2)*x)
      RETURN
      END
