program test_createFold
    implicit none
    character(len=100) :: cmdChar_folder, outdir, new_outdir, new_outdir1
    integer*4 :: i, j
    logical                      :: dirExists
    character(len=100) t

    outdir        = "../../outputs/new_version"
    write(new_outdir,*) adjustl(trim(outdir)),"/result_nc"
    write(cmdChar_folder,*) "mkdir ",adjustl(trim(new_outdir))
    ! write(cmdChar_folder,*) "shell [[ ! -e docs ]] && ",adjustl(trim(new_outdir1))

    inquire( file=trim(new_outdir)//'/.', exist=dirExists )  ! Works with gfortran, but not ifort
    ! inquire( directory=newDirPath, exist=dirExists )         ! Works with ifort, but not gfortran

    write(*,*) cmdChar_folder
    if (.not. dirExists) then
    i = system(adjustl(trim(cmdChar_folder)))
    endif
    t = "abcc"
    j = len(adjustl(trim(t)))
    new_outdir1 = adjustl(trim(t)) // adjustl(trim(t)) //"/hhhhh"
    ! call system('shell [[ ! -e docs ]] && mkdir docs')
    write(*,*)"test ...",i, dirExists, t, j, new_outdir1
end program test_createFold