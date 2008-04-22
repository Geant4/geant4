      subroutine ablainit()
      character*80 racine

      racine = "./data/"
      call init_evapora(racine)
      call inipace(racine)

      return
      end
