
program test_type
    use mod_type
    implicit none
    real abc(4), bbb, book_id
    INTEGER ggg(4)
    call init_it()
    write(*,*)book1%book_id
    write(*,*)book1%year
    deallocate(book1%year) 
    allocate(book1%year(4))
    book1%year = [1,2,3,4] 
    write(*,*)book1%year
    abc = [1.2,2.2,3.3,4.4]
    ggg = int(abc)
    write(*,*)ggg
    bbb = 0.
    if (bbb .ne. 0. ) write(*,*)"yyyy"
    write(*,*)bbb
    book_id = 1.2
    write(*,*)book_id,book1%book_id
end 