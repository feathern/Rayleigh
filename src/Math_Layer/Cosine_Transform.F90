Module Cosine_Transform
    Use, intrinsic :: iso_c_binding
	Implicit None
	Include 'fftw3.f03'

Contains

    
    Subroutine r2r_4D_fftw_forward(x,y, p_r2r_f)
		Real*8, Intent(InOut) :: x(:,:,:,:)
		Real*8, Intent(InOut) :: y(:,:,:,:)
		Integer*8, Intent(In), Optional :: p_r2r_f ! plan for the forward transform
		Type(C_Ptr)  :: fresh_plan ! plan for the forward transform
        Type(C_Ptr) :: pdum, qdum
        Real(8), Pointer :: my_in (:,:,:,:)
        Real(8), Pointer :: my_out(:,:,:,:)
        Integer, Dimension(4) :: dims
        Integer :: n1, n2, n3, n4 ! the four dimensions of array x
        dims = shape(x)
        n1 = dims(1)
        n2 = dims(2)
        n3 = dims(3)
        n4 = dims(4)
        
        pdum = fftw_alloc_real(product(dims))
        Call C_F_Pointer(pdum, my_in, [n1, n2, n3, n4])
        qdum = fftw_alloc_real(product(dims))
        Call C_F_Pointer(qdum, my_out,[n1, n2, n3, n4])

        my_in  = x

        fresh_plan = fftw_plan_many_r2r(1, [n1], n2*n3*n4, &
                                my_in, [n1], 1, n1, &
                                my_out, [n1], 1, n1, &
                                [FFTW_REDFT10], FFTW_estimate)

        Call fftw_execute_r2r (fresh_plan, my_in, my_out)
        
        y = my_out/n1   ! can optimize this out later
        Call fftw_free(pdum)
        Call fftw_free(qdum)

    End Subroutine r2r_4D_fftw_forward

 
    Subroutine r2r_4D_fftw_inverse(x, y, p_r2r_i)
		Real*8, Intent(InOut) :: x(:,:,:,:)
		Real*8, Intent(InOut) :: y(:,:,:,:)
		Integer*8, Intent(In), Optional :: p_r2r_i ! plan for the inverse transform
        Type(C_Ptr) :: fresh_plan
        Type(C_Ptr) :: pdum, qdum
        Real(8), Pointer :: my_in (:,:,:,:)
        Real(8), Pointer :: my_out(:,:,:,:)
        Integer, Dimension(4) :: dims
        Integer :: n1, n2, n3, n4 ! the four dimensions of array x
        dims = shape(x)
        n1 = dims(1)
        n2 = dims(2)
        n3 = dims(3)
        n4 = dims(4)


        pdum = fftw_alloc_real(product(dims))
        Call C_F_Pointer(pdum, my_in, [n1, n2, n3, n4])
        qdum = fftw_alloc_real(product(dims))
        Call C_F_Pointer(qdum, my_out,[n1, n2, n3, n4])

        my_in  = x
        fresh_plan = fftw_plan_many_r2r(1, [n1], n2*n3*n4, &
                                my_in, [n1], 1, n1, &
                                my_out, [n1], 1, n1, &
                                [FFTW_REDFT01], FFTW_estimate)

        Call fftw_execute_r2r (fresh_plan, my_in, my_out)
        y = my_out/2  ! optimize this out later         
        Call fftw_free(pdum)
        Call fftw_free(qdum)

    End Subroutine r2r_4D_fftw_inverse

End Module Cosine_Transform
