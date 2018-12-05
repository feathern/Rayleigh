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
        Character*120 :: tempfil='ginput'
        Integer :: gfunit   ! Generic file unit
        Integer :: etag     ! Endian check (should be 314)
        Integer :: fversion ! file version (for backwards compatibility)
        Integer :: imode    ! Input mode  (determines how modes are specified)
        Integer :: nlm      ! Number of l-m modes (imode = 1; assumed same ncheby_max)
        Integer :: inmax    ! Input max cheby degree (imode = 1; '')
        Integer :: nlmn     ! Number of l-m-n modes (imode = 2; completel arbitrary)
        Integer, Allocatable :: ilvalues(:)   ! Input l-values
        Integer, Allocatable :: imvalues(:)   ! Input m-values
        Integer, Allocatable :: invalues(:)   ! Input n-values

        Integer :: data_disp = 0    ! beginning of file's data section (in bytes)
            
        If (my_rank .eq. 0) Then
            Write(6,*)'Hello'
            Open(newunit=gfunit,file=tempfile,form='unformatted', status='old', access='stream')
            Read(gfunit)etag
            Read(gfunit)fversion
            Read(gfunit)imode

            If (imode .eq. 1) Then
                ! All l-m modes share the same nmax
                Read(gfunit)nlm
                Read(gfunit)inmax
                Allocate(ilvalues(1:nlm))
                Allocate(imvalues(1:nlm))
                Do i = 1, nlm
                    Read(gfunit)ilvalues(i)
                Enddo
                Do i = 1, nlm
                    Read(gfunit)imvalues(i)
                Enddo
                data_disp = 8 ! placeholder
                Close(gfunit)
            Endif
            If (imode .eq. 2) Then
                ! Each l-m-n mode is specified individually
                Read(gfunit)nlmn
                Allocate(ilvalues(1:nlmn))
                Allocate(imvalues(1:nlmn))
                Allocate(invalues(1:nlmn))
                Do i = 1, nlmn
                    Read(gfunit)ilvalues(i)
                Enddo
                Do i = 1, nlmn
                    Read(gfunit)imvalues(i)
                Enddo
                Do i = 1, nlmn
                    Read(gfunit)invalues(i)
                Enddo
                data_disp= 8
                Close(gfunit)
            Endif
        Endif

    End Subroutine Read_Generic

End Module Generic_Input
