! ======================================================================
! ======================================================================
!
SUBROUTINE mst(n,x,y,z,length,connections)

! generates a minimum spanning tree
! if 2d one of the x,y,z arrays must be a zero array
   IMPLICIT NONE
! inputs
   INTEGER, INTENT(in) :: n
   double precision, INTENT(in) :: x(1:n),y(1:n),z(1:n)
! outputs
   double precision, INTENT(out) :: length(1:n-1)      ! connections lengths
   double precision, intent(out) :: connections(1:6,1:n-1)
! internals
   integer, dimension(:), allocatable :: list
   integer, dimension(:,:), allocatable :: idents
   INTEGER ::i,j,memnum,told,ncount,imin,jmin,nlist
   INTEGER :: treemem(1:n)                             ! tree member
   double precision, dimension(:), allocatable :: sep

!
! how long is the list?
   nlist=0
   do i=1,n-1
     do j=i+1,n
       nlist=nlist + 1
     end do
   end do
!
   allocate(sep(1:nlist))
   allocate(list(1:nlist))
   allocate(idents(1:2,1:nlist))
!
! get distances between *all* points and put in a list
! I'm sure this could be clever and work-out what i and j are
! from the position in the array, but I can't be bothered
   nlist=0
   DO i=1,n-1
     DO j=i+1,n
       nlist=nlist + 1
       sep(nlist)=SQRT((x(i)-x(j))**2 + (y(i)-y(j))**2 + (z(i)-z(j))**2)
       if (sep(nlist)<=0.) then
          write (6,*) "One star needs shifting"
          write (6,*) i, x(i), y(i), j, x(j), y(j)
          stop
       end if
       list(nlist)=nlist
       idents(1,nlist)=i  !every sep gets an identity,
       idents(2,nlist)=j  !indicating included particle
     END DO
   END DO
!

! heapsort separations
   call heapsort(nlist,sep,list)
   if (sep(list(1))<=0.) then
!write(6,*) '******',idents(1,i),idents(2,j)
     write(6,*) 'zero length - same particles?',list(1); stop
   end if

!
! loop to allocate nodes
   length=0.
   connections=0. ! coordinates of connecting edges
   ncount=1
   treemem=0
   memnum=1
   nlist=0
   DO
! smallest distance
     nlist=nlist + 1
     length(ncount)=sep(list(nlist)) ! edge length
     imin=idents(1,list(nlist))      ! node i
     jmin=idents(2,list(nlist))      ! node j
     connections(1,ncount)=x(imin)
     connections(2,ncount)=y(imin)
     connections(3,ncount)=z(imin)
     connections(4,ncount)=x(jmin)
     connections(5,ncount)=y(jmin)
     connections(6,ncount)=z(jmin)
!
     IF (treemem(imin)/=0 .AND. treemem(imin)==treemem(jmin)) GOTO 23
!
! check if neither ptcl is already in a tree -------------------
     IF (treemem(imin)==0 .AND. treemem(jmin)==0) THEN
! not already part of a tree, so give a new treenum
       treemem(imin)=memnum
       treemem(jmin)=memnum
! advance memnum and ncount
       memnum=memnum + 1
       ncount=ncount + 1
       GOTO 23
     END IF
!
! check if one or the other is ---------------------------------
     IF (treemem(imin)==0 .AND. treemem(jmin)/=0) THEN
       treemem(imin)=treemem(jmin)
       ncount=ncount + 1
       GOTO 23
     END IF
!
     IF (treemem(imin)/=0 .AND. treemem(jmin)==0) THEN
       treemem(jmin)=treemem(imin)
       ncount=ncount + 1
       GOTO 23
     END IF
!
! are they in different trees? ---------------------------------
     IF (treemem(imin)/=0 .AND. treemem(jmin)/=treemem(imin)) THEN
! they'll all be part of this tree now
       told=treemem(jmin)
! loop over all current nodes
       DO i=1,n
         IF (treemem(i)==told) treemem(i)=treemem(imin)
       END DO
       ncount=ncount + 1
       GOTO 23
     END IF
!
23  IF (ncount==n) EXIT
!
   END DO
!

   RETURN
!
   deallocate(sep)
   deallocate(list)
   deallocate(idents)
   END SUBROUTINE mst
