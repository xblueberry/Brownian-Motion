c La linia que ve a continuacio facilita la identificacio de la posicio
c de cada caracter en la versio impresa del codi font:
c23456789012345678901234567890123456789012345678901234567890123456789012

C     BOX-MULLER. Daniel Izquierdo Juncas, 2014
C     Generacio de variables aleatories distribuides segons N(m,sigma) amb Box-Muller

      implicit none
      double precision m,sigma,x1
      integer i,N

      open(10, file="box-muller.txt")

      write(6,*) "METODE DE BOX-MULLER"
      
C     Tria de la mitjana i la variancia
      write(6,*) "Introdueix la mitjana:"
      read(5,*) m
      write(6,*) "Introdueix la variancia:"
      read(5,*) sigma

C     Nombres a generar (a mes nombres mes precisio) i reinicialitzacio
C     del generador de nombres aleatoris
      N=1000000
      call srand(0)

C     Nombres generats amb el metode de Box-Muller
      do i=1,N
        call boxmuller(m,sigma,x1)
        write(10,*) x1
      end do

      close (10)

      call histograma(N,m,sigma)

      write(6,*) "Les dades han estat guardades en un fitxer"

      end

C--------------------------------------------------------------------------------
C     SUBRUTINA BOX-MULLER. Arguments:
C     mitja, sigma: variables d'entrada
C     x1,x2:variables de sortida

      subroutine boxmuller(mitja,sigma,x1)
      
      double precision mitja,sigma
      double precision e1,e2,x1,pi
      pi=acos(-1.0d0)

      e1=rand()
      e2=rand()
      
C     Valors independents de la variable normal estandard
      x1=sqrt(-2.0d0*dlog(e1))*cos(2.0d0*pi*e2)
C     x2=sqrt(-2.0d0*dlog(e1))*sin(2.0d0*pi*e2)

C     Adaptant a la mitja i a la sigma de la gaussiana
      x1=mitja+sigma*x1

      return
      end


C---------------------------------------------------------------------------
C     SUBRUTINA HISTOGRAMA. Arguments:
C     N: nombre de numeros generats
C     m,sigma: parametres de la distribucio

      subroutine histograma(N,m,sigma)
      
      implicit none
      double precision m,sigma
      double precision inici,fi,h,nombre,hi,hf
      double precision  mitjana(1000),freq(1000)
      integer i,j,N,ncaixes

      open(10, file="box-muller.txt")
      open(11,file="histograma.txt")

      write(6,*) "Introdueix el nombre de caixes de l'histograma (fins a
     $ 1000):"
      read(5,*) ncaixes

C     Els extrems de l'histograma es delimiten amb el criteri de 3sigma
      inici=m-3.0d0*sigma
      fi=m+3.0d0*sigma
      h=(fi-inici)/dble(ncaixes)

C     Inicialitzacio dels vectors mitjana i frequencia
      do i=1,1000
        mitjana(i)=0.0d0
        freq(i)=0.0d0
      end do

C     Doble bucle que per a cada nombre del fitxer Box-Muller el classifica segons
C     l'interval on es troba i en suma 1 a la frequencia d'aquell
      do i=1,N
        read(10,*) nombre
        do j=1,ncaixes
C         Principi i final de cada interval
          hi=inici+(j-1)*h
          hf=inici+j*h

          if ((nombre.gt.hi).and.(nombre.lt.hf)) then
            mitjana(j)=mitjana(j)+nombre
            freq(j)=freq(j)+1.0d0
          end if
          
        end do
      end do

C     Es guarda la mitjana i la frequencia normalitzada al fitxer histograma
      do i=1,ncaixes
         mitjana(i)=mitjana(i)/freq(i)
         freq(i)=freq(i)/(N*h)
         write(11,*) mitjana(i), freq(i)
      end do

      close (10)
      close (11)
      
      end


