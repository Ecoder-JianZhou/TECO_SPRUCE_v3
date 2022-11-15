module mod_type
    type Books
        character(len = 50) :: title
        character(len = 50) :: author
        character(len = 150) :: subject
        integer :: book_id
        INTEGER, dimension (:), allocatable :: year
    end type Books
    real abc1
    type(Books) :: book1 
    real bdc
    contains
    subroutine init_it()
        implicit none
        book1%book_id = 1
        book1%author = "Jian"
        allocate(book1%year(5))
        book1%year = [1,2,3,4,5] 
    end subroutine init_it
end module mod_type