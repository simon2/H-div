program test_bubun
! Description:
!
! Author: am
!
! Host: aofd30
! Directory: /work2/am/12.Work11/21.Climatology/32.Test2_Bi-linear/src
!
! Revision history:
!  This file is created by /usr/local/mybin/nff.sh at 11:07 on 07-16-2011.
!  use
!  implicit none
  integer A(3,2)
  write(*,'(a)')'Program test_bubun starts.'
  write(*,*)''
  A(1,1)=1
  A(2,1)=2
  A(3,1)=3
  A(1,2)=4
  A(2,2)=5
  A(3,2)=6
  print *,"2D Array:"
  do i=1,3
    print *,(A(i,j),j=1,2)
  enddo
  print *
  print *,"1D Array:"
  j=1
  print *,'j= ',j
  call bubun(A(:,j),3)
  print *
  print *,"1D Array:"
  i=1
  print *,'i= ',i
  call bubun(A(i,:),2)
  write(*,*)
  write(*,'(a)')'Done program test_bubun.'

  ! real x
  ! real y
  ! real z
  x = 2.5
  y = 1.5
  z = log(y)
  print *, 'z=', z
end program test_bubun

subroutine bubun(B,im)
  integer B(im)
  do i=1,im
    print *,B(i)
  enddo
end subroutine bubun