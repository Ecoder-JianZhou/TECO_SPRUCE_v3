! program main
!     implicit none
!     character(len=100) :: str
!     integer :: x,y,z,id
!     character(len=3) :: xyz
 
!      read(*,'(a)')  str
!      id=index(str,'xyz',BACK=.False.)
!      xyz=str(id:id+2)
!      str(id:id+2)=''
!      id=index(str,'x')
!      x=str(id)
!      str(id)=' '
!      id=index(str,'y')
!      y=str(id)
!      str(id)=' '
!      id=index(str,'z')
!      z=str(id)
!      str(id)=' '
!      read(str,*) x,y,z
 
! end program main

program test_str
    implicit none
    integer a(4), b(3), c(4)
    a = [1,2,3,4]
    b = a(:3)
    c = a(:)
    write(*,*)b,c
end program test_str