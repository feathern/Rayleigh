Module Cosine_Transform
  Use, intrinsic :: iso_c_binding
  Implicit None
  Include 'fftw3.f03'

Contains



  Subroutine a_1d_dct_of(x,y) !< FIXME: change the ridiculous name
                              !! but get rid of the goto statement first
    Real*8, Intent(InOut) :: x(:)
    Real*8, Intent(InOut) :: y(:)
    Integer*8 :: p_r2r_f ! plan for the forward transform
    Type(C_Ptr)  :: fresh_plan ! plan for the forward transform
    Type(C_Ptr) :: pdum, qdum
    Real(8), Pointer :: my_in (:)
    Real(8), Pointer :: my_out(:)
    Integer :: n1! the four dimensions of array x
    n1 = size(x,1)

    pdum = fftw_alloc_real(int(n1, kind=C_SIZE_T))
    Call C_F_Pointer(pdum, my_in, [n1])
    qdum = fftw_alloc_real(int(n1, kind=C_SIZE_T))
    Call C_F_Pointer(qdum, my_out,[n1])

    my_in  = x

    fresh_plan = fftw_plan_many_r2r(1, [n1], 1, &
         my_in, [n1], 1, n1, &
         my_out, [n1], 1, n1, &
         [FFTW_REDFT10], FFTW_estimate) !< FIXME: call the regular planner, not the general one

    Call fftw_execute_r2r (fresh_plan, my_in, my_out)

    y = my_out/n1   ! can optimize this out later
    Call fftw_free(pdum)
    Call fftw_free(qdum)
   
    
  End Subroutine 


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

    pdum = fftw_alloc_real(int(product(dims), kind=C_SIZE_T))
    Call C_F_Pointer(pdum, my_in, [n1, n2, n3, n4])
    qdum = fftw_alloc_real(int(product(dims), kind=C_SIZE_T))
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


    pdum = fftw_alloc_real(int(product(dims), kind=C_SIZE_T))
    Call C_F_Pointer(pdum, my_in, [n1, n2, n3, n4])
    qdum = fftw_alloc_real(int(product(dims), kind=C_SIZE_T))
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

  Subroutine r2r_3D_nVar_dct_forward(x,y, howmany_vars, p_r2r_f)
    ! takes the DCT transform of equation_set(:,:)%rhs(:,:,:), 
    ! i.e. the DCT of an real(8) rank-3 array where the fastest index
    ! contains howmany_vars variables. 
    ! Therefore the size of the DCT is size(rhs,1)/howmany_vars
    ! and product(size(rhs))*howmany_vars transforms are taken.
    ! The result is dealiased by a 2/3 factor and stored in y.
    !
    ! OPTIMIZATIONS: 
    ! pdum/my_in and qdum/my_out are allocated everytime this routine is called
    ! extraction of y from my_out at the end
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Implicit None
    Real*8, Intent(InOut) :: x(:,:,:)
    Real*8, Intent(InOut) :: y(:,:,:)
    Integer, Intent(In) :: howmany_vars
    Integer*8, Intent(In), Optional :: p_r2r_f ! plan for the forward transform
    Type(C_Ptr)  :: fresh_plan ! plan for the forward transform
    Type(C_Ptr) :: pdum, qdum
    Real(8), Pointer :: my_in (:,:,:)
    Real(8), Pointer :: my_out(:,:,:,:)
    Integer, Dimension(3) :: dims
    Integer :: n1, n2, n3, n1_dealias, n1_all, ista
    Integer :: ivar, j,k
    dims = shape(x)
    n1_all = dims(1)
    n1 = n1_all/howmany_vars
    n1_dealias = (n1*2)/3
    n2 = dims(2)
    n3 = dims(3)

    pdum = fftw_alloc_real(int(product(dims), kind=C_SIZE_T)) ! C_SIZE_T probably kind=8
    Call C_F_Pointer(pdum, my_in,  [n1_all, n2, n3])
    qdum = fftw_alloc_real(int(product(dims), kind=C_SIZE_T))
    Call C_F_Pointer(qdum, my_out, [n1, howmany_vars, n2, n3])

    ! copy x to pdum
    my_in = x

    fresh_plan = fftw_plan_many_r2r(1, [n1], howmany_vars*n2*n3, &
         my_in,  [n1], 1, n1, &
         my_out, [n1], 1, n1, &
         [FFTW_REDFT10], FFTW_estimate)

    Call fftw_execute_r2r (fresh_plan, my_in, my_out)


    Do k = 1, n3
       Do j = 1, n2
          Do ivar = 1,howmany_vars
             ista = (ivar-1)*n1_dealias
             y(ista+1:ista+n1_dealias,j,k) = my_out(1:n1_dealias,ivar,j,k)/n1
          End Do
       End Do
    End Do

    Call fftw_free(pdum)
    Call fftw_free(qdum)

  End Subroutine r2r_3D_nVar_dct_forward
End Module Cosine_Transform
