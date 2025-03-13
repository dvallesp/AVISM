

!********************************************************************************
     SUBROUTINE SMOOTH(NX,NY,NZ)
!Smooth quantities to a coarse grid: (NCOX*NCOY*NCOZ), with NCOX < NHYX
!********************************************************************************
     USE COMMONDATA

     IMPLICIT NONE
!input variables
     INTEGER NX,NY,NZ
!local variables
     INTEGER I, J, K
     INTEGER II, JJ, KK
     REAL*4 RADX1, RADXNX, RADY1, RADYNY, RADZ1, RADZNZ
     REAL*4 WEI, RX, RY, RZ

      !smallest/largest coordinates of the coarse grid
      RADX1=RADX(1)-0.5*DX
      RADXNX=RADX(NX)+0.5*DX
      RADY1=RADY(1)-0.5*DY
      RADYNY=RADY(NY)+0.5*DY
      RADZ1=RADZ(1)-0.5*DZ
      RADZNZ=RADZ(NZ)+0.5*DZ

      WEI=(DX0*DY0*DZ0)/(DX*DY*DZ) !dilution factor

      DO KK=1, NHYZ !level 0 grid
      DO JJ=1, NHYY
      DO II=1, NHYX

         RX=RADX0(II)
         RY=RADY0(JJ)
         RZ=RADZ0(KK)

     !move by LENG particles outside the grid, check if happens!
         IF(RX.LT.RADX1) WRITE(*,*) '      WARNING: RX<RADX1 in SMOOTH'
         IF(RX.GE.RADXNX) WRITE(*,*) '     WARNING: RX>=RADXNX in SMOOTH'
         IF(RY.LT.RADY1) WRITE(*,*) '      WARNING: RY<RADY1 in SMOOTH'
         IF(RY.GE.RADYNY) WRITE(*,*) '     WARNING: RY>=RADYNY in SMOOTH'
         IF(RZ.LT.RADZ1) WRITE(*,*) '      WARNING: RZ<RADZ1 in SMOOTH'
         IF(RZ.GE.RADZNZ) WRITE(*,*)'      WARNING: RZ>RADZNZ in SMOOTH'

         I=INT(((RX-RADX(1))/DX)+0.5) + 1 !coarse (NCOX^3) grid
         J=INT(((RY-RADY(1))/DY)+0.5) + 1
         K=INT(((RZ-RADZ(1))/DZ)+0.5) + 1

         U1DMCO(I,J,K)=U1DMCO(I,J,K)+U1DM(II,JJ,KK)*WEI
         U1GCO(I,J,K)=U1GCO(I,J,K)+U1G(II,JJ,KK)*WEI
         U1SCO(I,J,K)=U1SCO(I,J,K)+U1S(II,JJ,KK)*WEI
         U2GCO(I,J,K)=U2GCO(I,J,K)+U2G(II,JJ,KK)*WEI
         U3GCO(I,J,K)=U3GCO(I,J,K)+U3G(II,JJ,KK)*WEI
         U4GCO(I,J,K)=U4GCO(I,J,K)+U4G(II,JJ,KK)*WEI

      ENDDO
      ENDDO
      ENDDO

!*************************************************************************
END SUBROUTINE SMOOTH
!*************************************************************************


!********************************************************************************
     SUBROUTINE SMOOTH_2(NX,NY,NZ)
!Version just for U1G, U2G, U3G, U4G
!Smooth quantities to a coarse grid: (NCOX*NCOY*NCOZ), with NCOX < NHYX
!********************************************************************************
     USE COMMONDATA

     IMPLICIT NONE
!input variables
     INTEGER NX,NY,NZ
!local variables
     INTEGER I, J, K
     INTEGER II, JJ, KK
     REAL*4 RADX1, RADXNX, RADY1, RADYNY, RADZ1, RADZNZ
     REAL*4 WEI, RX, RY, RZ

      !smallest/largest coordinates of the coarse grid
      RADX1=RADX(1)-0.5*DX
      RADXNX=RADX(NX)+0.5*DX
      RADY1=RADY(1)-0.5*DY
      RADYNY=RADY(NY)+0.5*DY
      RADZ1=RADZ(1)-0.5*DZ
      RADZNZ=RADZ(NZ)+0.5*DZ

      WEI=(DX0*DY0*DZ0)/(DX*DY*DZ) !dilution factor

      DO KK=1, NHYZ !level 0 grid
      DO JJ=1, NHYY
      DO II=1, NHYX

         RX=RADX0(II)
         RY=RADY0(JJ)
         RZ=RADZ0(KK)

     !move by LENG particles outside the grid, check if happens!
         IF(RX.LT.RADX1) WRITE(*,*) '      WARNING: RX<RADX1 in SMOOTH'
         IF(RX.GE.RADXNX) WRITE(*,*) '     WARNING: RX>=RADXNX in SMOOTH'
         IF(RY.LT.RADY1) WRITE(*,*) '      WARNING: RY<RADY1 in SMOOTH'
         IF(RY.GE.RADYNY) WRITE(*,*) '     WARNING: RY>=RADYNY in SMOOTH'
         IF(RZ.LT.RADZ1) WRITE(*,*) '      WARNING: RZ<RADZ1 in SMOOTH'
         IF(RZ.GE.RADZNZ) WRITE(*,*)'      WARNING: RZ>RADZNZ in SMOOTH'

         I=INT(((RX-RADX(1))/DX)+0.5) + 1 !coarse (NCOX^3) grid
         J=INT(((RY-RADY(1))/DY)+0.5) + 1
         K=INT(((RZ-RADZ(1))/DZ)+0.5) + 1

         U1GCO(I,J,K)=U1GCO(I,J,K)+U1G(II,JJ,KK)*WEI
         U2GCO(I,J,K)=U2GCO(I,J,K)+U2G(II,JJ,KK)*WEI
         U3GCO(I,J,K)=U3GCO(I,J,K)+U3G(II,JJ,KK)*WEI
         U4GCO(I,J,K)=U4GCO(I,J,K)+U4G(II,JJ,KK)*WEI

      ENDDO
      ENDDO
      ENDDO

!*************************************************************************
      END SUBROUTINE SMOOTH_2
!*************************************************************************




!***************************************************************************
SUBROUTINE VMESH(IR, NX, NY, NZ, DXR, DYR, DZR, RX1, RY1, RZ1, NX0, NY0, NZ0)
!***************************************************************************
    USE COMMONDATA
    IMPLICIT NONE
    !input variables:
    INTEGER:: IR, NX0, NY0, NZ0
    INTEGER:: NX, NY, NZ 
    REAL*4:: DXR, DYR, DZR, RX1, RY1, RZ1

    !local variables
    INTEGER:: IX0, JY0, KZ0, IXR, JYR, KZR, IXR1, JYR1, KZR1
    INTEGER:: IPX, IPY, IPZ
    INTEGER:: N1, N2, N3, IPA, LOW1, LOW2
    INTEGER, ALLOCATABLE:: FLAGR(:,:,:), FLAGAMR(:,:,:)

    WRITE(*,*) '     Starting VMESH: IR=',IR
    WRITE(*,*) '        NX0, NY0, NZ0:', NX0, NY0, NZ0
    WRITE(*,*) '        NX, NY, NZ:',NX, NY, NZ

    ALLOCATE(UDMR(NX, NY, NZ))
    ALLOCATE(UGR(NX, NY, NZ))
    ALLOCATE(USR(NX, NY, NZ))
    ALLOCATE(DIVR(NX, NY, NZ))

    UDMR(:,:,:)=0.
    UGR(:,:,:)=0.
    USR(:,:,:)=0.
    DIVR(:,:,:)=0.
  
    !This variables are used for computing mean value of AMR variables when creating l>0 mesh
    ALLOCATE(FLAGAMR(NX, NY, NZ))
    FLAGAMR(:,:,:)=1 !updated only if VMESH is called

    ALLOCATE(FLAGR(NX0, NY0, NZ0))
    FLAGR(:,:,:)=0 ! =0-->cell in IR=1 is not refined; =1 --> cell is refined

    LOW1=SUM(NPATCH(0:IR-1))+1
    LOW2=SUM(NPATCH(0:IR))
    DO IPA=LOW1,LOW2

          N1 = PATCHNX(IPA)
          N2 = PATCHNY(IPA)
          N3 = PATCHNZ(IPA)

          IXR1=INT(((PATCHRX(IPA)-0.5*DXR)-RX1)/DXR+0.5)
          JYR1=INT(((PATCHRY(IPA)-0.5*DYR)-RY1)/DYR+0.5)
          KZR1=INT(((PATCHRZ(IPA)-0.5*DZR)-RZ1)/DZR+0.5)


          DO IPZ=1, N3
             DO IPY=1, N2
                DO IPX=1, N1

                   IXR=IXR1+IPX
                   JYR=JYR1+IPY
                   KZR=KZR1+IPZ

                   IF(IXR .LT.1 .OR. IXR .GT. NX) WRITE(*,*) '    WARNING ON IXR!', IXR, IXR1, IPX
                   IF(JYR .LT.1 .OR. JYR .GT. NY) WRITE(*,*) '    WARNING ON JYR!', JYR, JYR1, IPY
                   IF(KZR .LT.1 .OR. KZR .GT. NZ) WRITE(*,*) '    WARNING ON KZR!', KZR, KZR1, IPZ

                   !MEAN VALUE
                   FLAGAMR(IXR, JYR, KZR)=FLAGAMR(IXR, JYR, KZR)+1
                   UDMR(IXR, JYR, KZR)=UDMR(IXR, JYR, KZR)+U11DM(IPX, IPY, IPZ, IPA)
                   UGR(IXR, JYR, KZR)=UGR(IXR, JYR, KZR)+U11G(IPX, IPY, IPZ, IPA)
                   USR(IXR, JYR, KZR)=USR(IXR, JYR, KZR)+U11S(IPX, IPY, IPZ, IPA)
                   DIVR(IXR, JYR, KZR)=DIVR(IXR, JYR, KZR)+DIVER(IPX, IPY, IPZ, IPA)

                   IX0=INT((IXR-1)/2.)+1
                   JY0=INT((JYR-1)/2.)+1
                   KZ0=INT((KZR-1)/2.)+1
                   FLAGR(IX0, JY0, KZ0)=1 

                ENDDO
             ENDDO
          ENDDO

       ENDDO !loop on patches


       DO KZR=1, NZ
          DO JYR=1, NY
             DO IXR=1, NX
                IF(FLAGAMR(IXR,JYR,KZR) .GE. 1) THEN
                   UDMR(IXR, JYR, KZR)=UDMR(IXR, JYR, KZR)/REAL(FLAGAMR(IXR,JYR,KZR))
                   UGR(IXR, JYR, KZR)=UGR(IXR, JYR, KZR)/REAL(FLAGAMR(IXR,JYR,KZR))
                   USR(IXR, JYR, KZR)=USR(IXR, JYR, KZR)/REAL(FLAGAMR(IXR,JYR,KZR))
                   DIVR(IXR, JYR, KZR)=DIVR(IXR, JYR, KZR)/REAL(FLAGAMR(IXR,JYR,KZR))
                ENDIF
             ENDDO
          ENDDO
       ENDDO

       DO KZ0=1, NZ0
         DO JY0=1, NY0
            DO IX0=1, NX0
               IF(FLAGR(IX0, JY0, KZ0) .EQ. 0 ) THEN
                  IXR1=2*IX0-1
                  JYR1=2*JY0-1
                  KZR1=2*KZ0-1
                  DO KZR=KZR1, KZR1+1
                     DO JYR=JYR1, JYR1+1
                        DO IXR=IXR1, IXR1+1
                           UDMR(IXR, JYR, KZR)= U1DM(IX0, JY0, KZ0)
                           UGR(IXR, JYR, KZR)= U1G(IX0, JY0, KZ0)
                           USR(IXR, JYR, KZR)= U1S(IX0, JY0, KZ0)
                           DIVR(IXR, JYR, KZR)=DIVER0(IX0, JY0, KZ0)
                        ENDDO
                     ENDDO
                  ENDDO
                ENDIF
             ENDDO
          ENDDO
       ENDDO

    WRITE(*,*) '     End VMESH....'

    DEALLOCATE(FLAGR)
    DEALLOCATE(FLAGAMR)
    
!***************************************************************************
  END SUBROUTINE VMESH
!***************************************************************************