program main
    implicit none
    character(len=100) :: str
    integer :: x,y,z,id
    character(len=3) :: xyz
 
     read(*,'(a)')  str
     id=index(str,'xyz',BACK=.False.)
     xyz=str(id:id+2)
     str(id:id+2)=''
     id=index(str,'x')
     x=str(id)
     str(id)=' '
     id=index(str,'y')
     y=str(id)
     str(id)=' '
     id=index(str,'z')
     z=str(id)
     str(id)=' '
     read(str,*) x,y,z
 
end program main