      subroutine g4dpmjet_initialise_block_data ()
      
C
C ------------------------------------------------------------------------------
C
      EXTERNAL         RUNTT
      EXTERNAL         NONAME
      EXTERNAL         ZK
      EXTERNAL         BLKD43
      EXTERNAL         REACCH
      EXTERNAL         JTDATA
      EXTERNAL         BLKD41
      EXTERNAL         BOOKLE
      EXTERNAL         BLKD46
      EXTERNAL         BLKD47
      EXTERNAL         QPROP
      EXTERNAL         RADINI
      EXTERNAL         HADINI
      EXTERNAL         POMEN
      EXTERNAL         PYDATA
C
C
C The follow line is to force the mathlib least-squares polynomial fit to be
C loaded.  This is required as the -L and -lmathlib -lkernlib flags in the
C load process appear to come before the G4ParamType1GlauberDataSet.o object
C is loaded, resulting in dlsqpm not being found.
C
      external     dlsqpm
     
      return
      
      end
