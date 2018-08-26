C La linia que ve a continuacio facilita la identificacio de la posicio
C de cada caracter en la versio impresa del codi font:
C23456789012345678901234567890123456789012345678901234567890123456789012

C     MOVIMENT BROWNIA AMB UNA PARET. Francesc Salvat. Barcelona. March, 2008.
C     Comentat i modificat per Daniel Izquierdo Juncas. 2014
C
      implicit none

      integer NPART,ISEED1,ISEED2
      integer NP, NSTEPS, NSHOTS
      integer I, ISHOT, ISTEP

      double precision RLAMB, TEMP, STEPT, PARET
      double precision XAV, SDEV, DIFFC, TIME, X, Y
      double precision R1, R2
      double precision XMEAN, X2MEAN, YMEAN, Y2MEAN, DIFFX, DIFFY

C     Maxim nombre de particules que accepta el programa
      PARAMETER (NPART=10000)

C     X(npart) i Y(npart) son les coordenades cartesianes de cadascuna de les particules
      DIMENSION X(NPART),Y(NPART)

C     Inicialitzacio del generador de nombres aleatoris
      COMMON/RSEED/ISEED1,ISEED2
      ISEED1=1
      ISEED2=1

    1 CONTINUE
      WRITE(6,'(/1X,
     1  '' Introdueix el nombre de particules: '')')
      READ(5,*) NP
C     Control de l'entrada de nombres invalids de particules
      IF(NP.LT.1.OR.NP.GT.NPART) GO TO 1

      WRITE(6,'(/1X,
     1  '' Introdueix el parametre de friccio lambda: '')')
      READ(5,*) RLAMB

    2 CONTINUE
      WRITE(6,'(/1X,
     1  '' Introdueix la temperatura: '')')
      READ(5,*) TEMP
      IF(TEMP.LT.0) GO TO 2

      WRITE(6,'(/1X,
     1  '' Introdueix el pas de temps: '')')
      READ(5,*) STEPT

      WRITE(6,'(/1X,
     1  '' Introdueix el nombre de passos de temps entre cada captura: '
     $')')
      READ(5,*) NSTEPS

      WRITE(6,'(/1X,
     1  '' Introdueix el nombre de captures: '')')
      READ(5,*) NSHOTS

C     Introduim la condicio de contorn de la paret
      write(6,*) "Introdueix la posicio de la paret perpendicular a l'ei
     $x x (no la situis massa lluny): "
      read(5,*) PARET


C     Format de presentacio per consola dels parametres introduits
      WRITE(6,1000) NP,RLAMB,TEMP,STEPT,NSTEPS,NSHOTS,PARET
 1000 FORMAT(
     1  /1X,'# Nombre de particules ....................',I6,
     1  /1X,'# Parametre de friccio ....................',1P,E14.6,
     1  /1X,'# Temperatura .............................',1P,E14.6,
     1  /1X,'# Pas de temps.............................',1P,E14.6,
     1  /1X,'# Nombre de passos de temps per captura ...',I6,
     1  /1X,'# Nombre de captures.......................',I6,
     1  /1X,'# Posicio de la paret......................',1P,E14.6)

C     Inicialitzacio de la mitjana (nula) i de la varian‡a (sense considerar
C     la constant de Boltzmann)
      XAV=0.0D0
      SDEV=SQRT(2.0D0*TEMP*STEPT/RLAMB)
      WRITE(6,'(/,'' # Desplacament mig ='',1P,E13.6)') SDEV

C     Coeficient de difusio d'Einstein (sense considerar la constant de Boltzmann)
      DIFFC=TEMP/RLAMB
      WRITE(6,'('' # Coeficient de difusio (teoric) ='',1P,E13.6,/)')
     $DIFFC

C     Inicialitzacio del temps i de les posicions de cada particula a 0
      TIME=0.0d0
      DO I=1,NPART
        X(I)= 0.0D0
        Y(I)= 0.0D0
      ENDDO

      OPEN(10,FILE='posicions.txt')
      OPEN(11,FILE='mitjanes.txt')

      DO ISHOT=1,NSHOTS

C       Calcul de la posicio de cada particula amb Box-Muller en una captura,
C       es te en compte el nombre de passos de temps que s'han realitzat
        DO ISTEP=1,NSTEPS
          DO I=1,NP
            CALL boxmuller(XAV,SDEV,R1,R2)
            X(I)=X(I)+R1
C           Introduccio de la condicio de la paret, quan la component x es
C           superior a la posicio de la paret s'inverteix la direccio del
C           moviment
            if(X(I).ge.PARET) then
              X(I)=PARET-(X(I)-PARET)
            end if
            
            Y(I)=Y(I)+R2
          ENDDO
        ENDDO

C       Temps en una captura
        TIME=TIME+STEPT*NSTEPS

C       Calcul de la posicio mitja i de la mitja quadratica de X i Y en
C       una captura
        XMEAN=0.0D0
        YMEAN=0.0D0
        X2MEAN=0.0D0
        Y2MEAN=0.0D0
        DO I=1,NP
          XMEAN=XMEAN+X(I)
          YMEAN=YMEAN+Y(I)
          X2MEAN=X2MEAN+X(I)*X(I)
          Y2MEAN=Y2MEAN+Y(I)*Y(I)
        ENDDO
        XMEAN=XMEAN/DBLE(NP)
        YMEAN=YMEAN/DBLE(NP)
        X2MEAN=X2MEAN/DBLE(NP)
        Y2MEAN=Y2MEAN/DBLE(NP)

C       Calcul dels coeficients de difusio de X i Y
        DIFFX=(X2MEAN-XMEAN**2)/(2.0D0*TIME)
        DIFFY=(Y2MEAN-YMEAN**2)/(2.0D0*TIME)

C       Es guarda en l'arxiu de posicions els parametres seguents que encap‡alen
C       cada captura
        WRITE(10,1001) ISHOT,TIME,XMEAN,YMEAN,DIFFX,DIFFY
 1001   FORMAT(
     1    /1X,'# Nombre de captures ...........',I6,
     1    /1X,'# Temps ........................',1P,E14.6,
     1    /1X,'# Posicio mitja ................',1P,
     2       '(',E14.6,',',E14.6,')',
     1    /1X,'# Coeficient de difusio ........',1P,
     2       '(',E14.6,',',E14.6,')')

        WRITE(6,1002) ISHOT,TIME,XMEAN,YMEAN,DIFFX,DIFFY

C       Les dades seguents es guarden a l'arxiu mitjana
        WRITE(11,1002) ISHOT,TIME,XMEAN,YMEAN,DIFFX,DIFFY
 1002   FORMAT(1X,I4,1P,7E13.5)

C       Es guarden les posicions a l'arxiu corresponent
        DO I=1,NP
          WRITE(10,'(1P,3E14.6)') X(I),Y(I)
        ENDDO
        write(10,*) ' '

      ENDDO

      CLOSE(10)
      CLOSE(11)

      WRITE(6,'(/,'' # Desplacament mig ='',1P,E13.6)') SDEV
      WRITE(6,'('' # Coeficient de difusio (teoric) ='',1P,E13.6,/)')
     $DIFFC
      WRITE(6,1001) ISHOT,TIME,XMEAN,YMEAN,DIFFX,DIFFY

      END
C--------------------------------------------------------------------------------
C     SUBRUTINA BOX-MULLER. Arguments:
C     mitja, sigma: variables d'entrada
C     x1,x2:variables de sortida

      subroutine boxmuller(mitja,sigma,x1,x2)

      implicit none
      double precision mitja,sigma
      double precision e1,e2,x1,x2,pi
      double precision RAND

      external RAND
      pi=acos(-1.0d0)

      e1=rand(1.0d0)
      e2=rand(2.0d0)

C     Valors independents de la variable normal estandard
      x1=sqrt(-2.0d0*dlog(e1))*cos(2.0d0*pi*e2)
      x2=sqrt(-2.0d0*dlog(e1))*sin(2.0d0*pi*e2)

C     Adaptant a la mitja i a la sigma de la gaussiana
      x1=mitja+sigma*x1
      x2=mitja+sigma*x2

      return
      end

C  *********************************************************************
C                         FUNCTION RAND (Random number generator)
C  *********************************************************************
C     NOTA: s'ha intentat utilitzar la funcio intrinseca rand() del Fortran
C     a la subrutina de Box-Muller, pero l'output d'algunes coordenades
C     anaven a infinit o a NaN

      FUNCTION RAND(DUMMY)
C  This is an adapted version of subroutine RANECU written by F. James
C  (Comput. Phys. Commun. 60 (1990) 329-344), which has been modified to
C  give a single random number at each call.
C
C  The 'seeds' ISEED1 and ISEED2 must be initialised in the main program
C  and transferred through the named common block /RSEED/.
C
C  Some compilers incorporate an intrinsic random number generator with
C  the same name (but with different argument lists). To avoid conflict,
C  it is advisable to declare RAND as an external function in all sub-
C  programs that call it.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      PARAMETER (USCALE=1.0D0/2.147483563D9)
      COMMON/RSEED/ISEED1,ISEED2
C
      I1=ISEED1/53668
      ISEED1=40014*(ISEED1-I1*53668)-I1*12211
      IF(ISEED1.LT.0) ISEED1=ISEED1+2147483563
C
      I2=ISEED2/52774
      ISEED2=40692*(ISEED2-I2*52774)-I2*3791
      IF(ISEED2.LT.0) ISEED2=ISEED2+2147483399
C
      IZ=ISEED1-ISEED2
      IF(IZ.LT.1) IZ=IZ+2147483562
      RAND=IZ*USCALE
C
      RETURN
      END
