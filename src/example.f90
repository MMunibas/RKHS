program example
use RKHS            ! This module needs to be used by your code
implicit none

!Some variables are defined here
integer, parameter :: D = 2       ! number of dimensions
real(kind(0d0))    :: x(D)        ! coordinate array (even when D = 1, this needs to be an array)
real(kind(0d0))    :: f           ! will store the value of the model function 
real(kind(0d0))    :: dfdx(D)     ! will store partial first derivatives of the model function with respect to any component of x
real(kind(0d0))    :: d2fdx2(D,D) ! will store the Hessian matrix of the model function
type(kernel)       :: K           ! The kernel type is needed to set up and evaluate a RKHS model

!for calculating only some of the partial derivatives/some entries of the Hessian
logical            :: dfdx_mask(D)
logical            :: d2fdx2_mask(D,D)

! we choose some arbitrary point at which the model function should be evaluated
x = (/1.5d0,0.5d0/)

!Below follow the 3 (or 4) steps necessary to set up and evaluate a RKHS model. 
!======================================================================================================================
! 1) read training data
!======================================================================================================================
call K%read_grid("multidimensional-grid.csv") 
!   Each line of the grid file must contain 
!   x1 , x2 , x3 , ... , xn , f 
!   where the xi are the individual coordinates and f the reference function values. 
!   Further, the file needs to be sorted ascendingly in this order: x1>x2>x3>...>xn
!   For example (3 dimensions):
!
!   x1_1 , x2_1 , x3_1 , f_1_1_1
!   x1_1 , x2_1 , x3_2 , f_1_1_2
!   x1_1 , x2_1 , x3_3 , f_1_1_3
!     .  ,   .  ,   .  ,    . 
!     .  ,   .  ,   .  ,    . 
!     .  ,   .  ,   .  ,    . 
!   x1_1 , x2_1 , x3_N3, f_1_1_N3
!   x1_1 , x2_2 , x3_1 , f_1_2_1
!   x1_1 , x2_2 , x3_2 , f_1_2_2
!     .  ,   .  ,   .  ,    . 
!     .  ,   .  ,   .  ,    . 
!     .  ,   .  ,   .  ,    . 
!   x1_1 , x2_2 , x3_N3 , f_1_2_N3
!     .  ,   .  ,   .  ,    . 
!     .  ,   .  ,   .  ,    . 
!     .  ,   .  ,   .  ,    . 
!   x1_N1, x2_N2, x3_N3, f_N1_N2_N3
!
!   where xi_1 < xi_2 < xi_3 < ... < xi_Ni and Ni is the amount of grid points in dimension i
!   Note that the grid file MUST CONTAIN THE FULL GRID SPECIFICATION, even if data is not available
!   at some grid points (else there will be reading issues). If data is not available, place the string
!   "NaN" in the place of f for the associated grid point.
!   See also the provided grid file "multidimensional-grid.csv" as a reference!


!======================================================================================================================
! 2) choose one-dimensional kernel functions 
!======================================================================================================================
call K%k1d(1)%init(EXPONENTIAL_DECAY_N3_KERNEL)   ! choose one-dimensional kernel for dimension 1
K%k1d(1)%par(1) = 2d0 !OPTIONAL: If the kernel has parameters, you can set them like this (parameters default to 1)
call K%k1d(2)%init(TAYLOR_SPLINE_N3_KERNEL)       ! choose one-dimensional kernel for dimension 2
!   For higher-dimensional kernels, this pattern is simply repeated, e.g.
!   call K%k1d(d)%init(KERNEL_TYPE) ! choose one-dimensional kernel for dimension d
!   The available KERNEL_TYPEs are:
!
!   RECIPROCAL_POWER_N2_M6_KERNEL   (no parameters)
!   RECIPROCAL_POWER_N3_M6_KERNEL   (no parameters)
!   RECIPROCAL_POWER_N2_M5_KERNEL   (no parameters)
!   RECIPROCAL_POWER_N3_M5_KERNEL   (no parameters)
!   RECIPROCAL_POWER_N2_M4_KERNEL   (no parameters)
!   RECIPROCAL_POWER_N3_M4_KERNEL   (no parameters)
!   RECIPROCAL_POWER_N2_M3_KERNEL   (no parameters)
!   RECIPROCAL_POWER_N3_M3_KERNEL   (no parameters)
!   RECIPROCAL_POWER_N2_M2_KERNEL   (no parameters)
!   RECIPROCAL_POWER_N3_M2_KERNEL   (no parameters)
!   RECIPROCAL_POWER_N2_M1_KERNEL   (no parameters)
!   RECIPROCAL_POWER_N3_M1_KERNEL   (no parameters)
!   RECIPROCAL_POWER_N2_M0_KERNEL   (no parameters)
!   RECIPROCAL_POWER_N3_M0_KERNEL   (no parameters) 
!   EXPONENTIAL_DECAY_N2_KERNEL     ( 1 parameter )  
!   EXPONENTIAL_DECAY_N3_KERNEL     ( 1 parameter )   
!   TAYLOR_SPLINE_N2_KERNEL         (no parameters)   
!   TAYLOR_SPLINE_N3_KERNEL         (no parameters)
!   LAPLACIAN_KERNEL                ( 1 parameter )        
!   BERNOULLI_N2_KERNEL             (no parameters)
! 


!======================================================================================================================
! 3) calculate kernel coefficients
!======================================================================================================================
call K%calculate_coefficients_fast() ! once all kernels are initialized, the coefficients can be calculated (uses Tensor product)
!   ATTENTION: If the kernel matrix is ill-conditioned, the code will crash at this point and print the error message:
!
!   "Matrix M is not positive definite! If you are sure M is correct, maybe the matrix is ill conditioned. 
!   Try to use the Tikhonov regularization procedure."
!
!   In that case, the fast method for calculating coefficients unfortunately is not applicable. However, it is still possible
!   to calculate the coefficients using Tikhonov regularization with:
!
!   call K%calculate_coefficients_slow(lambda) 
!
!   where lambda is typically a very small value, e.g. 1d-16. If you still receive the same error message, try to increase lambda!
!
!   At this stage, the model function can already be evaluated using the slow evaluation method:
call K%evaluate_slow(x,f)
write(*,*) "Slow evaluation at point ", x, " gives ", f
!   Or, if derivatives are also desired, we just need to supply an appropriate array for storing them
call K%evaluate_slow(x,f,dfdx)
write(*,*) "Slow evaluation at point ", x, " gives derivative of f with respect to x(1): ", dfdx(1)
write(*,*) "Slow evaluation at point ", x, " gives derivative of f with respect to x(2): ", dfdx(2)
!   Or, if the Hessian is also desired, we just need to supply an appropriate array for storing them
call K%evaluate_slow(x,f,dfdx,d2fdx2)
write(*,*) "Slow evaluation at point ", x, " gives the Hessian: "
write(*,*) d2fdx2(1,:)
write(*,*) d2fdx2(2,:)

!   The fast evaluation method needs additional setup:
!======================================================================================================================
! 4) calculate lookup table for fast evaluation
!======================================================================================================================
call K%calculate_sums() ! this call will calculate the complete lookup table
!   Now, it is also possible to use the fast evaluation method:
call K%evaluate_fast(x,f)
write(*,*) "Fast evaluation at point ", x, " gives ", f
!   Or, if derivatives are also desired, we just need to supply an appropriate array for storing them
call K%evaluate_fast(x,f,dfdx)
write(*,*) "Fast evaluation at point ", x, " gives derivative of f with respect to x(1): ", dfdx(1)
write(*,*) "Fast evaluation at point ", x, " gives derivative of f with respect to x(2): ", dfdx(2)
!   Or, if the Hessian is also desired, we just need to supply an appropriate array for storing them
call K%evaluate_fast(x,f,dfdx,d2fdx2)
write(*,*) "Fast evaluation at point ", x, " gives the Hessian: "
write(*,*) d2fdx2(1,:)
write(*,*) d2fdx2(2,:)

!   Note that the results from the slow evaluation and fast evaluation method might disagree slightly due to rounding errors


!======================================================================================================================
! Although the set up is complete at this stage, there is still some useful additional functionality implemented:
!======================================================================================================================
!   At any stage of the setup, the kernel (and all data that has been calculated up to this point, e.g. coefficients,
!   lookup tables, etc.) can be saved to a binary file
call K%save_to_file("test.kernel")
!   If the kernel is not needed anymore, all associated memory can be freed with
call K%free()

!   If a kernel is available as a binary file, it can be loaded and directly evaluated without repetition of the 
!   previous setup steps
call K%load_from_file("test.kernel")
call K%evaluate_fast(x,f,dfdx,d2fdx2)
write(*,*) "Fast evaluation at point ", x, " gives ", f
write(*,*) "Fast evaluation at point ", x, " gives derivative of f with respect to x(1): ", dfdx(1)
write(*,*) "Fast evaluation at point ", x, " gives derivative of f with respect to x(2): ", dfdx(2)
write(*,*) "Fast evaluation at point ", x, " gives the Hessian: "
write(*,*) d2fdx2(1,:)
write(*,*) d2fdx2(2,:)
!   This allows easy sharing of PES or models for other functions with other people: Simply send them the
!   .kernel file! It also means that e.g. in a Molecular Dynamics code, the PES can be loaded from a .kernel
!   file without the need to perform the complete setup.
!   The training data can always be recovered in human-readable format from a .kernel file:
call K%write_grid("multidimensional-grid-RECOVERED.csv") 

!   Note that for efficiency reasons, it is possible to provide a "mask" when calculating 
!   partial derivatives or the Hessian of the model function. In case only certain derivatives
!   or entries of the Hessian are required, only those entries are calculated, which saves
!   computational time
dfdx   = 0d0 !reset to 0
d2fdx2 = 0d0 !reset to 0
!only the derivative with respect to x(1) will be calculated
dfdx_mask        = (/.true.,.false./) 
!only the diagonal entries of the Hessian will be calculated
d2fdx2_mask(1,:) = (/.true.  , .false./)
d2fdx2_mask(2,:) = (/.false. , .true. /)
call K%evaluate_fast(x,f,dfdx,d2fdx2,dfdx_mask,d2fdx2_mask)
write(*,*) 
write(*,*) "Evaluation with supplied mask:"
write(*,*) "Fast evaluation at point ", x, " gives ", f
write(*,*) "Fast evaluation at point ", x, " gives derivative of f with respect to x(1): ", dfdx(1)
write(*,*) "Fast evaluation at point ", x, " gives derivative of f with respect to x(2): ", dfdx(2)
write(*,*) "Fast evaluation at point ", x, " gives the Hessian: "
write(*,*) d2fdx2(1,:)
write(*,*) d2fdx2(2,:)

! Note: masks can also be supplied in the same fashion for the slow evaluation method

call K%free() ! free resources

!======================================================================================================================
! IMPORTANT: 
!======================================================================================================================
!   The RKHS module uses file unit 30 for all file i/o by default. If this causes problems with your existing code,
!   you can modify the i/o default unit. Go to the source file RKHS.f90 and search for 
!   "integer, parameter :: io_unit = 30" 
!   in "module file_io" (it's located at the very top of the RKHS.f90 file) and modify the file unit to whatever unit you want.
end program example
