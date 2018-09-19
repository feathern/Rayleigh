Module Generic_Input
    Use Timers, Only : stopwatch
    Use ProblemSize
    Use Parallel_Framework
    Use Spherical_Buffer
    Use Linear_Solve, Only : get_all_rhs
    Use SendReceive
    Use ISendReceive
    Use Controls
    Use MPI_BASE
    Use Chebyshev_Polynomials_Alt

    Use BufferedOutput
    ! Simple Checkpointing Module
    ! Uses MPI-IO to split writing of files amongst rank zero processes from each row
    Implicit None

Contains

    Subroutine Read_Generic()
        Implicit None
        If (my_rank .eq. 0) Then
            Write(6,*)'Hello'
        Endif

    End Subroutine Read_Generic

End Module Generic_Input
