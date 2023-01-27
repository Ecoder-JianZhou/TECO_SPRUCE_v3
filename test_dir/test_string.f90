program test_string
    implicit none
    character(len=1000) aaa

    aaa = "abcddfdfdfdf"
    write(*,*)trim(aaa), "ggggggggg"

end program test_string