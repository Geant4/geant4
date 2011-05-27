      subroutine ablainit()
      character*80 racine
      common/rseed/iseed1,iseed2
      racine = "./data/"
      iseed1=666
      iseed2=777
      call init_evapora(racine)
      call inipace(racine)

      return
      end
