module usleep
    use iso_c_binding, only: c_int

    ! Define interface for sleep function
    interface
        subroutine sleep(milliseconds) bind(C, name="sleep")
            import :: c_int
            integer(c_int), value :: milliseconds
        end subroutine sleep
    end interface

contains

    ! Subroutine to sleep for milliseconds
    subroutine sleep_milliseconds(milliseconds)
        integer(c_int), intent(in) :: milliseconds
        call sleep(milliseconds)
    end subroutine sleep_milliseconds

end module usleep