!***************************************************************************
   SUBROUTINE INDEXX(n,arr,indx)
!***************************************************************************
      INTEGER n,indx(n),M,NSTACK
      REAL arr(n)
      PARAMETER (M=7,NSTACK=50)
      INTEGER i,indxt,ir,itemp,j,jstack,k,l,istack(NSTACK)
      REAL a
      do 11 j=1,n
        indx(j)=j
11    continue
      jstack=0
      l=1
      ir=n
1     if(ir-l.lt.M)then
        do 13 j=l+1,ir
          indxt=indx(j)
          a=arr(indxt)
          do 12 i=j-1,1,-1
            if(arr(indx(i)).le.a)goto 2
            indx(i+1)=indx(i)
12        continue
          i=0
2         indx(i+1)=indxt
13      continue
        if(jstack.eq.0)return
        ir=istack(jstack)
        l=istack(jstack-1)
        jstack=jstack-2
      else
        k=(l+ir)/2
        itemp=indx(k)
        indx(k)=indx(l+1)
        indx(l+1)=itemp
        if(arr(indx(l+1)).gt.arr(indx(ir)))then
          itemp=indx(l+1)
          indx(l+1)=indx(ir)
          indx(ir)=itemp
        endif
        if(arr(indx(l)).gt.arr(indx(ir)))then
          itemp=indx(l)
          indx(l)=indx(ir)
          indx(ir)=itemp
        endif
        if(arr(indx(l+1)).gt.arr(indx(l)))then
          itemp=indx(l+1)
          indx(l+1)=indx(l)
          indx(l)=itemp
        endif
        i=l+1
        j=ir
        indxt=indx(l)
        a=arr(indxt)
3       continue
          i=i+1
        if(arr(indx(i)).lt.a)goto 3
4       continue
          j=j-1
        if(arr(indx(j)).gt.a)goto 4
        if(j.lt.i)goto 5
        itemp=indx(i)
        indx(i)=indx(j)
        indx(j)=itemp
        goto 3
5       indx(l)=indx(j)
        indx(j)=indxt
        jstack=jstack+2
        if(jstack.gt.NSTACK) then
         write(*,*) 'NSTACK too small in indexx'
         stop
        endif
        if(ir-i+1.ge.j-l)then
          istack(jstack)=ir
          istack(jstack-1)=i
          ir=j-1
        else
          istack(jstack)=j-1
          istack(jstack-1)=l
          l=i
        endif
      endif
      goto 1
!***************************************************************************
      END SUBROUTINE INDEXX
!***************************************************************************


!***************************************************************************
      SUBROUTINE REAL_MAKE_PERIODIC(NXX,NYY,NZZ,NPLUS,FIELD)
!***************************************************************************
      IMPLICIT NONE
      INTEGER NXX,NYY,NZZ,NPLUS
      REAL FIELD(-NPLUS+1:NXX+NPLUS,-NPLUS+1:NYY+NPLUS,-NPLUS+1:NZZ+NPLUS)
      INTEGER IX,JY,KZ

      DO KZ = -NPLUS+1, NZZ+NPLUS
        DO JY = -NPLUS+1, NYY+NPLUS
          DO IX = -NPLUS+1, NXX+NPLUS
            IF (IX .LT. 1) THEN
              FIELD(IX,JY,KZ) = FIELD(NXX+IX,JY,KZ)
            ELSE IF (IX .GT. NXX) THEN
              FIELD(IX,JY,KZ) = FIELD(IX-NXX,JY,KZ)
            END IF
          END DO
        END DO
      END DO  

      DO KZ = -NPLUS+1, NZZ+NPLUS
        DO JY = -NPLUS+1, NYY+NPLUS
          DO IX = -NPLUS+1, NXX+NPLUS
            IF (JY .LT. 1) THEN
              FIELD(IX,JY,KZ) = FIELD(IX,NYY+JY,KZ)
            ELSE IF (JY .GT. NYY) THEN
              FIELD(IX,JY,KZ) = FIELD(IX,JY-NYY,KZ)
            END IF
          END DO
        END DO
      END DO  

      DO KZ = -NPLUS+1, NZZ+NPLUS
        DO JY = -NPLUS+1, NYY+NPLUS
          DO IX = -NPLUS+1, NXX+NPLUS
            IF (KZ .LT. 1) THEN
              FIELD(IX,JY,KZ) = FIELD(IX,JY,NZZ+KZ)
            ELSE IF (KZ .GT. NZZ) THEN
              FIELD(IX,JY,KZ) = FIELD(IX,JY,KZ-NZZ)
            END IF
          END DO
        END DO
      END DO

!***************************************************************************
      END SUBROUTINE REAL_MAKE_PERIODIC
!***************************************************************************

!***************************************************************************
      SUBROUTINE INT_MAKE_PERIODIC(NXX,NYY,NZZ,NPLUS,FIELD)
!***************************************************************************
      IMPLICIT NONE
      INTEGER NXX,NYY,NZZ,NPLUS
      INTEGER FIELD(-NPLUS+1:NXX+NPLUS,-NPLUS+1:NYY+NPLUS,-NPLUS+1:NZZ+NPLUS)
      INTEGER IX,JY,KZ

      DO KZ = -NPLUS+1, NZZ+NPLUS
        DO JY = -NPLUS+1, NYY+NPLUS
          DO IX = -NPLUS+1, NXX+NPLUS
            IF (IX .LT. 1) THEN
              FIELD(IX,JY,KZ) = FIELD(NXX+IX,JY,KZ)
            ELSE IF (IX .GT. NXX) THEN
              FIELD(IX,JY,KZ) = FIELD(IX-NXX,JY,KZ)
            END IF
          END DO
        END DO
      END DO  

      DO KZ = -NPLUS+1, NZZ+NPLUS
        DO JY = -NPLUS+1, NYY+NPLUS
          DO IX = -NPLUS+1, NXX+NPLUS
            IF (JY .LT. 1) THEN
              FIELD(IX,JY,KZ) = FIELD(IX,NYY+JY,KZ)
            ELSE IF (JY .GT. NYY) THEN
              FIELD(IX,JY,KZ) = FIELD(IX,JY-NYY,KZ)
            END IF
          END DO
        END DO
      END DO  

      DO KZ = -NPLUS+1, NZZ+NPLUS
        DO JY = -NPLUS+1, NYY+NPLUS
          DO IX = -NPLUS+1, NXX+NPLUS
            IF (KZ .LT. 1) THEN
              FIELD(IX,JY,KZ) = FIELD(IX,JY,NZZ+KZ)
            ELSE IF (KZ .GT. NZZ) THEN
              FIELD(IX,JY,KZ) = FIELD(IX,JY,KZ-NZZ)
            END IF
          END DO
        END DO
      END DO

!***************************************************************************
      END SUBROUTINE INT_MAKE_PERIODIC
!***************************************************************************