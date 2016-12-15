!###############################################################################
!                                VERSION 1.0
! Copyright (C) 2016 Oliver T. Unke <oliver.unke@googlemail.com>
!
! This program is free software. It comes without any warranty, to
! the extent permitted by applicable law. You can redistribute it
! and/or modify it under the terms of the MIT license (see below).
!###############################################################################
!-------------------------------------------------------------------------------
!MIT License
!
!Copyright (c) 2016 Oliver T. Unke
!
!Permission is hereby granted, free of charge, to any person obtaining a copy
!of this software and associated documentation files (the "Software"), to deal
!in the Software without restriction, including without limitation the rights
!to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
!copies of the Software, and to permit persons to whom the Software is
!furnished to do so, subject to the following conditions:
!
!The above copyright notice and this permission notice shall be included in all
!copies or substantial portions of the Software.
!
!THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
!IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
!FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
!AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
!LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
!OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
!SOFTWARE.
!-------------------------------------------------------------------------------


!-----------------------------------------------------------------------
!> @author
!> Oliver T. Unke, University of Basel
!
!> @brief
!> File input/output unit
!
!> @details
!> This module contains the standard file i/o unit that is used by the RKHS
!! module. In case the file unit conflicts with any existing code, it can
!! be easily changed here.
!-----------------------------------------------------------------------
module file_io
    implicit none
    integer, parameter :: io_unit = 30
end module file_io


!-----------------------------------------------------------------------
!> @author
!> Oliver T. Unke, University of Basel
!
!> @brief
!> Laplacian kernel
!
!> @details
!> This module contains all functions necessary to evaluate the Laplacian kernel
!! given by exp(||x-x'||/sigma). Note that sigma can be set to
!! any positive value, but is initialized automatically with the value 1.
!! It is not recommended for interpolating potential energy surfaces, but it can
!! be useful for machine learning purposes. Note however, that in order for the
!! decomposition into f2k and f3k functions to work, it is necessary that all
!! input variables are positive. Therefore, the kernel is only defined in the
!! interval [0,infinity), whereas normally a Laplacian kernel is valid in the
!! interval (-infinity,infinity).
!-----------------------------------------------------------------------
module kernel_laplacian
    implicit none
    integer, parameter :: M2   = 1
    integer, parameter :: Npar = 1
    real(kind(0d0)), parameter :: p21 = 1d0
    
    contains        
    !direct kernel impelementation
    pure real(kind(0d0)) function k(x1,x2,par)
        implicit none
        real(kind(0d0)), intent(in) :: x1,x2
        real(kind(0d0))             :: xs,xl
        real(kind(0d0)), dimension(:), intent(in) :: par 
        !find larger/smaller of x1 and x2
        if(x1 <= x2) then
            xs = x1
            xl = x2
        else
            xs = x2
            xl = x1
        end if
        k = exp(-(xl-xs)/par(1))
    end function k
    
    pure real(kind(0d0)) function dk(x1,x2,par)
        implicit none
        real(kind(0d0)), intent(in) :: x1,x2
        real(kind(0d0)), dimension(:), intent(in) :: par 
        !find larger/smaller of x1 and x2
        if(x1 <= x2) then
            dk = exp(-(x2-x1)/par(1))/par(1)
        else
            dk = -exp(-(x1-x2)/par(1))/par(1)
        end if
    end function dk
    
    pure real(kind(0d0)) function d2k(x1,x2,par)
        implicit none
        real(kind(0d0)), intent(in) :: x1,x2
        real(kind(0d0)), dimension(:), intent(in) :: par 
        !find larger/smaller of x1 and x2
        if(x1 <= x2) then
            d2k = exp(-(x2-x1)/par(1))/par(1)**2
        else
            d2k = exp(-(x1-x2)/par(1))/par(1)**2
        end if
    end function d2k
        
    !f2 and f3 functions and their derivatives
    pure real(kind(0d0)) function f21(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        f21 = exp(x/par(1))
    end function f21
    
    pure real(kind(0d0)) function df21(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        df21 = exp(x/par(1))/par(1)
    end function df21
    
    pure real(kind(0d0)) function d2f21(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        d2f21 = exp(x/par(1))/par(1)**2
    end function d2f21
    
    pure real(kind(0d0)) function f31(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        f31 = exp(-x/par(1))
    end function f31
    
    pure real(kind(0d0)) function df31(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        df31 = -exp(-x/par(1))/par(1)
    end function df31 
    
    pure real(kind(0d0)) function d2f31(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        d2f31 = exp(-x/par(1))/par(1)**2
    end function d2f31 
       
end module kernel_laplacian


!-----------------------------------------------------------------------
!> @author
!> Oliver T. Unke, University of Basel
!
!> @brief
!> Bernoulli polynomial kernel with n = 2
!
!> @details
!> This module contains all functions necessary to evaluate the Bernoulli
!! kernel with n = 2. Note that this kernel is only defined for values in the interval [0,1].
!! It is only applicable for periodic functions, such as sin(2*pi*x).
!-----------------------------------------------------------------------
module kernel_bernoulli_2
    implicit none
    integer, parameter :: M2   = 12
    integer, parameter :: Npar = 0
    real(kind(0d0)), parameter :: p21  =  1d0/12d0,&
                                  p22  = -1d0/24d0,&
                                  p23  = -1d0/24d0,&
                                  p24  = -1d0/12d0,&
                                  p25  =  1d0/12d0,&
                                  p26  = -1d0/4d0 ,&
                                  p27  =  1d0/4d0 ,&
                                  p28  = -1d0/24d0,&
                                  p29  = -1d0/24d0,&
                                  p210 =  1d0/6d0,&
                                  p211 =  1d0/6d0,&
                                  p212 = -1d0/4d0
                                  
    contains        
    !direct kernel impelementation
    pure real(kind(0d0)) function k(x1,x2,par)
        implicit none
        real(kind(0d0)), intent(in) :: x1,x2
        real(kind(0d0))             :: xs,xl
        real(kind(0d0)), dimension(:), intent(in) :: par 
        !find larger/smaller of x1 and x2
        if(x1 <= x2) then
            xs = x1
            xl = x2
        else
            xs = x2
            xl = x1
        end if
        k = 0.1D1 / 0.720D3 + xs * xl / 0.12D2 - xs ** 2 / 0.24D2 - xl ** 2 &
            / 0.24D2 - xs ** 3 / 0.12D2 + xl ** 3 / 0.12D2 - xs * xl ** 2 / &
            0.4D1 + xs ** 2 * xl / 0.4D1 - xs ** 4 / 0.24D2 - xl ** 4 / 0.24D2 &
            + xs ** 3 * xl / 0.6D1 + xs * xl ** 3 / 0.6D1 - xs ** 2 * xl ** 2 / 0.4D1
    end function k
    
    pure real(kind(0d0)) function dk(x1,x2,par)
        implicit none
        real(kind(0d0)), intent(in) :: x1,x2
        real(kind(0d0)), dimension(:), intent(in) :: par 
        !find larger/smaller of x1 and x2
        if(x1 <= x2) then
            dk = -x1 / 0.12D2 + x2 / 0.12D2 - x1 ** 3 / 0.6D1 + x1 ** 2 * x2 &
                 / 0.2D1 - x1 * x2 ** 2 / 0.2D1 + x2 ** 3 / 0.6D1 - x1 ** 2 / 0.4D1 &
                  + x2 * x1 / 0.2D1 - x2 ** 2 / 0.4D1  
        else
            dk = x2 / 0.12D2 - x1 / 0.12D2 + x1 ** 2 / 0.4D1 - x2 * x1 / 0.2D1 + &
                 x2 ** 2 / 0.4D1 - x1 ** 3 / 0.6D1 + x2 ** 3 / 0.6D1 + x1 ** 2 * x2 &
                 / 0.2D1 - x1 * x2 ** 2 / 0.2D1
        end if
    end function dk
    
    pure real(kind(0d0)) function d2k(x1,x2,par)
        implicit none
        real(kind(0d0)), intent(in) :: x1,x2
        real(kind(0d0)), dimension(:), intent(in) :: par 
        !find larger/smaller of x1 and x2
        if(x1 <= x2) then
            d2k = -0.1D1 / 0.12D2 - x1 ** 2 / 0.2D1 + x2 * x1 - x2 ** 2 / 0.2D1 &
                  - x1 / 0.2D1 + x2 / 0.2D1
  
        else
            d2k = -0.1D1 / 0.12D2 + x1 / 0.2D1 - x2 / 0.2D1 - x1 ** 2 / 0.2D1 &
                  + x2 * x1 - x2 ** 2 / 0.2D1  
        end if
    end function d2k
        
    !f2 and f3 functions and their derivatives
    pure real(kind(0d0)) function f21(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        f21 = x
    end function f21
    
    pure real(kind(0d0)) function df21(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        df21 = 1d0
    end function df21
    
    pure real(kind(0d0)) function d2f21(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        d2f21 = 0d0
    end function d2f21
    
    pure real(kind(0d0)) function f31(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        f31 = x
    end function f31
    
    pure real(kind(0d0)) function df31(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        df31 = 1d0
    end function df31 
    
    pure real(kind(0d0)) function d2f31(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        d2f31 = 0d0
    end function d2f31 
    
    pure real(kind(0d0)) function f22(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        f22 = x**2 - 1d0/30d0
    end function f22
    
    pure real(kind(0d0)) function df22(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        df22 = 2d0*x
    end function df22
    
    pure real(kind(0d0)) function d2f22(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        d2f22 = 2d0
    end function d2f22
    
    pure real(kind(0d0)) function f32(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        f32 = 1d0
    end function f32
    
    pure real(kind(0d0)) function df32(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        df32 = 0d0
    end function df32
    
    pure real(kind(0d0)) function d2f32(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        d2f32 = 0d0
    end function d2f32
    
    pure real(kind(0d0)) function f23(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        f23 = 1d0
    end function f23
    
    pure real(kind(0d0)) function df23(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        df23 = 0d0
    end function df23
    
    pure real(kind(0d0)) function d2f23(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        d2f23 = 0d0
    end function d2f23
    
    pure real(kind(0d0)) function f33(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        f33 = x**2
    end function f33
    
    pure real(kind(0d0)) function df33(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par
        df33 = 2d0*x
    end function df33
    
    pure real(kind(0d0)) function d2f33(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par
        d2f33 = 2d0
    end function d2f33
    
    pure real(kind(0d0)) function f24(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        f24 = x**3
    end function f24
    
    pure real(kind(0d0)) function df24(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        df24 = 3d0*x**2
    end function df24
    
    pure real(kind(0d0)) function d2f24(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        d2f24 = 6d0*x
    end function d2f24
    
    pure real(kind(0d0)) function f34(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        f34 = 1d0
    end function f34
    
    pure real(kind(0d0)) function df34(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par
        df34 = 0d0
    end function df34
    
    pure real(kind(0d0)) function d2f34(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par
        d2f34 = 0d0
    end function d2f34
    
    pure real(kind(0d0)) function f25(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        f25 = 1d0
    end function f25
    
    pure real(kind(0d0)) function df25(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        df25 = 0d0
    end function df25
    
    pure real(kind(0d0)) function d2f25(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        d2f25 = 0d0
    end function d2f25
    
    pure real(kind(0d0)) function f35(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        f35 = x**3
    end function f35
    
    pure real(kind(0d0)) function df35(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        df35 = 3d0*x**2
    end function df35
    
    pure real(kind(0d0)) function d2f35(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        d2f35 = 6d0*x
    end function d2f35
    
    pure real(kind(0d0)) function f26(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        f26 = x
    end function f26
    
    pure real(kind(0d0)) function df26(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par        
        df26 = 1d0
    end function df26
    
    pure real(kind(0d0)) function d2f26(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par        
        d2f26 = 0d0
    end function d2f26
    
    pure real(kind(0d0)) function f36(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        f36 = x**2
    end function f36
    
    pure real(kind(0d0)) function df36(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        df36 = 2d0*x
    end function df36
    
    pure real(kind(0d0)) function d2f36(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        d2f36 = 2d0
    end function d2f36
    
    pure real(kind(0d0)) function f27(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        f27 = x**2
    end function f27
    
    pure real(kind(0d0)) function df27(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        df27 = 2d0*x
    end function df27
    
    pure real(kind(0d0)) function d2f27(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        d2f27 = 2d0
    end function d2f27
    
    pure real(kind(0d0)) function f37(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        f37 = x
    end function f37
    
    pure real(kind(0d0)) function df37(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par         
        df37 = 1d0
    end function df37
    
    pure real(kind(0d0)) function d2f37(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par         
        d2f37 = 0d0
    end function d2f37
    
    pure real(kind(0d0)) function f28(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        f28 = x**4
    end function f28
    
    pure real(kind(0d0)) function df28(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        df28 = 4d0*x**3
    end function df28
    
    pure real(kind(0d0)) function d2f28(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        d2f28 = 12d0*x**2
    end function d2f28
    
    pure real(kind(0d0)) function f38(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        f38 = 1d0
    end function f38
    
    pure real(kind(0d0)) function df38(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par         
        df38 = 0d0
    end function df38
    
    pure real(kind(0d0)) function d2f38(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par         
        d2f38 = 0d0
    end function d2f38
    
    pure real(kind(0d0)) function f29(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        f29 = 1d0
    end function f29
    
    pure real(kind(0d0)) function df29(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        df29 = 0d0
    end function df29
    
    pure real(kind(0d0)) function d2f29(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        d2f29 = 0d0
    end function d2f29
    
    pure real(kind(0d0)) function f39(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        f39 = x**4
    end function f39
    
    pure real(kind(0d0)) function df39(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        df39 = 4d0*x**3
    end function df39
    
    pure real(kind(0d0)) function d2f39(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        d2f39 = 12d0*x**2
    end function d2f39
    
    pure real(kind(0d0)) function f210(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        f210 = x**3
    end function f210
    
    pure real(kind(0d0)) function df210(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        df210 = 3d0*x**2
    end function df210
    
    pure real(kind(0d0)) function d2f210(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        d2f210 = 6d0*x
    end function d2f210
    
    pure real(kind(0d0)) function f310(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        f310 = x
    end function f310
    
    pure real(kind(0d0)) function df310(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        df310 = 1d0
    end function df310
    
    pure real(kind(0d0)) function d2f310(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        d2f310 = 0d0
    end function d2f310
    
    pure real(kind(0d0)) function f211(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        f211 = x
    end function f211
    
    pure real(kind(0d0)) function df211(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        df211 = 1d0
    end function df211
    
    pure real(kind(0d0)) function d2f211(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        d2f211 = 0d0
    end function d2f211
    
    pure real(kind(0d0)) function f311(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        f311 = x**3
    end function f311
    
    pure real(kind(0d0)) function df311(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        df311 = 3d0*x**2
    end function df311
    
    pure real(kind(0d0)) function d2f311(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        d2f311 = 6d0*x
    end function d2f311
    
    pure real(kind(0d0)) function f212(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        f212 = x**2
    end function f212
    
    pure real(kind(0d0)) function df212(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        df212 = 2d0*x
    end function df212
    
    pure real(kind(0d0)) function d2f212(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        d2f212 = 2d0
    end function d2f212
    
    pure real(kind(0d0)) function f312(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        f312 = x**2
    end function f312
    
    pure real(kind(0d0)) function df312(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        df312 = 2d0*x
    end function df312
    
    pure real(kind(0d0)) function d2f312(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        d2f312 = 2d0
    end function d2f312
     
end module kernel_bernoulli_2


!-----------------------------------------------------------------------
!> @author
!> Oliver T. Unke, University of Basel
!
!> @brief
!> Exponential decay kernel with n = 2
!
!> @details
!> This module contains all functions necessary to evaluate the exponential
!! decay kernel with n = 2. For values larger than the greatest point in the grid,
!! the kernel decays exponentially with exp(-beta*x). Note that beta can be set to
!! any positive value, but is initialized automatically with the value 1.
!! This type of kernel is recommended for short-range intermolecular interactions, which
!! often decay exponentially. The kernel is only defined for values in the interval [0,infinity).
!-----------------------------------------------------------------------
module kernel_exp_2
    implicit none
    integer, parameter :: M2   = 2
    integer, parameter :: Npar = 1
    real(kind(0d0)), parameter :: p21 = 4d0,&
                                  p22 = 4d0
    contains        
    !direct kernel impelementation
    pure real(kind(0d0)) function k(x1,x2,par)
        implicit none
        real(kind(0d0)), intent(in) :: x1,x2
        real(kind(0d0))             :: xs,xl
        real(kind(0d0)), dimension(:), intent(in) :: par 
        !find larger/smaller of x1 and x2
        if(x1 <= x2) then
            xs = x1
            xl = x2
        else
            xs = x2
            xl = x1
        end if
        k = 4d0*exp(-par(1)*xl)/par(1)**3 * &
            (par(1)*(xl - xs) + 2d0)
    end function k
    
    pure real(kind(0d0)) function dk(x1,x2,par)
        implicit none
        real(kind(0d0)), intent(in) :: x1,x2
        real(kind(0d0)), dimension(:), intent(in) :: par 
        !find larger/smaller of x1 and x2
        if(x1 <= x2) then
            dk = -4d0*exp(-par(1)*x2)/par(1)**2
        else
            dk = -4d0*exp(-par(1)*x1)/par(1)**2 * &
                  (par(1)*(x1-x2) + 1d0)
        end if
    end function dk
    
    pure real(kind(0d0)) function d2k(x1,x2,par)
        implicit none
        real(kind(0d0)), intent(in) :: x1,x2
        real(kind(0d0)), dimension(:), intent(in) :: par 
        !find larger/smaller of x1 and x2
        if(x1 <= x2) then
            d2k = 0d0
        else
            d2k = 4d0*exp(-par(1)*x1) * (x1 - x2)
        end if
    end function d2k
        
    !f2 and f3 functions and their derivatives
    pure real(kind(0d0)) function f21(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        f21 = 2d0-par(1)*x
    end function f21
    
    pure real(kind(0d0)) function df21(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        df21 = -par(1)
    end function df21
    
    pure real(kind(0d0)) function d2f21(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        d2f21 = 0d0
    end function d2f21
    
    pure real(kind(0d0)) function f22(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        f22 = par(1)
    end function f22
    
    pure real(kind(0d0)) function df22(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        df22 = 0d0
    end function df22
    
    pure real(kind(0d0)) function d2f22(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        d2f22 = 0d0
    end function d2f22
    
    pure real(kind(0d0)) function f31(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        f31 = exp(-par(1)*x)/par(1)**3
    end function f31
    
    pure real(kind(0d0)) function df31(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        df31 = -exp(-par(1)*x)/par(1)**2
    end function df31
    
    pure real(kind(0d0)) function d2f31(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        d2f31 = exp(-par(1)*x)/par(1)
    end function d2f31
    
    pure real(kind(0d0)) function f32(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        f32 = x*exp(-par(1)*x)/par(1)**3
    end function f32

    pure real(kind(0d0)) function df32(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        df32 = exp(-par(1)*x)/par(1)**2 * (1d0/par(1) - x)
    end function df32
    
    pure real(kind(0d0)) function d2f32(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        d2f32 = exp(-par(1)*x)/par(1) * (x - 2d0/par(1))
    end function d2f32
    
end module kernel_exp_2


!-----------------------------------------------------------------------
!> @author
!> Oliver T. Unke, University of Basel
!
!> @brief
!> Exponential decay kernel with n = 3
!
!> @details
!> This module contains all functions necessary to evaluate the exponential
!! decay kernel with n = 3. For values larger than the greatest point in the grid,
!! the kernel decays exponentially with exp(-beta*x). Note that beta can be set to
!! any positive value, but is initialized automatically with the value 1.
!! This type of kernel is recommended for short-range intermolecular interactions, which
!! often decay exponentially. The kernel is only defined for values in the interval [0,infinity).
!-----------------------------------------------------------------------
module kernel_exp_3
    implicit none
    integer, parameter :: M2   = 5
    integer, parameter :: Npar = 1
    real(kind(0d0)), parameter :: p21 = 18d0,&
                                  p22 = 18d0,&
                                  p23 = 18d0,&
                                  p24 = 18d0,&
                                  p25 = 18d0
    contains
            
    pure real(kind(0d0)) function k(x1,x2,par)
        implicit none
        real(kind(0d0)), intent(in) :: x1,x2
        real(kind(0d0))             :: xs,xl
        real(kind(0d0)), dimension(:), intent(in) :: par 
        !find larger/smaller of x1 and x2
        if(x1 <= x2) then
            xs = x1
            xl = x2
        else
            xs = x2
            xl = x1
        end if
        k = 18d0*exp(-par(1)*xl)/par(1)**5 * (par(1)**2*(xl**2 - 2d0*xl*xs + xs**2) &
                                              + 6d0*par(1)*(xl - xs) + 12d0)
    end function k
    
    !direct kernel impelementation
    pure real(kind(0d0)) function dk(x1,x2,par)
        implicit none
        real(kind(0d0)), intent(in) :: x1,x2
        real(kind(0d0)), dimension(:), intent(in) :: par 
        !find larger/smaller of x1 and x2
        if(x1 <= x2) then
            dk = -36d0*exp(-par(1)*x2)/par(1)**4 * (par(1)*(x2-x1)+3d0)
        else
            dk = -18d0*exp(-par(1)*x1)/par(1)**4 * (par(1)**2 * (x1**2 - 2d0*x1*x2 + x2**2) &
                                                    + 4d0*par(1)*(x1 - x2) + 6d0)
        end if
    end function dk
            
    pure real(kind(0d0)) function d2k(x1,x2,par)
        implicit none
        real(kind(0d0)), intent(in) :: x1,x2
        real(kind(0d0)), dimension(:), intent(in) :: par 
        !find larger/smaller of x1 and x2
        if(x1 <= x2) then
            d2k = 36d0*exp(-par(1)*x2)/par(1)**3 
        else
            d2k = 18d0*exp(-par(1)*x1)/par(1)**3 * ((par(1)*(x1-x2))**2 + 2d0*par(1)*(x1-x2) + 2d0)
        end if
    end function d2k
        
    !f2 and f3 functions and their derivatives
    pure real(kind(0d0)) function f21(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        f21 = 12d0 - 6d0*par(1)*x
    end function f21
    
    pure real(kind(0d0)) function df21(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        df21 = -6d0*par(1)
    end function df21
    
    pure real(kind(0d0)) function d2f21(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        d2f21 = 0d0
    end function d2f21
    
    pure real(kind(0d0)) function f22(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        f22 = 6d0*par(1)
    end function f22
    
    pure real(kind(0d0)) function df22(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        df22 = 0d0
    end function df22
    
    pure real(kind(0d0)) function d2f22(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        d2f22 = 0d0
    end function d2f22
    
    pure real(kind(0d0)) function f23(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        f23 = (par(1)*x)**2
    end function f23
    
    pure real(kind(0d0)) function df23(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        df23 = 2d0*par(1)**2*x
    end function df23
    
    pure real(kind(0d0)) function d2f23(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        d2f23 = 2d0*par(1)**2
    end function d2f23
    
    pure real(kind(0d0)) function f24(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        f24 = -2d0*par(1)**2*x
    end function f24
    
    pure real(kind(0d0)) function df24(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        df24 = -2d0*par(1)**2
    end function df24
    
    pure real(kind(0d0)) function d2f24(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        d2f24 = 0d0
    end function d2f24
    
    pure real(kind(0d0)) function f25(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        f25 = par(1)**2
    end function f25
    
    pure real(kind(0d0)) function df25(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        df25 = 0d0
    end function df25
    
    pure real(kind(0d0)) function d2f25(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        d2f25 = 0d0
    end function d2f25
    
    pure real(kind(0d0)) function f31(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        f31 = exp(-par(1)*x)/par(1)**5
    end function f31
    
    pure real(kind(0d0)) function df31(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        df31 = -exp(-par(1)*x)/par(1)**4
    end function df31
    
    pure real(kind(0d0)) function d2f31(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        d2f31 = exp(-par(1)*x)/par(1)**3
    end function d2f31
    
    pure real(kind(0d0)) function f32(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        f32 = x*exp(-par(1)*x)/par(1)**5
    end function f32

    pure real(kind(0d0)) function df32(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        df32 = exp(-par(1)*x)/par(1)**4 * (1d0/par(1) - x)
    end function df32
    
    pure real(kind(0d0)) function d2f32(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        d2f32 = exp(-par(1)*x)/par(1)**3 * (x - 2d0/par(1))
    end function d2f32
    
    pure real(kind(0d0)) function f33(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        f33 = exp(-par(1)*x)/par(1)**5
    end function f33

    pure real(kind(0d0)) function df33(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        df33 = -exp(-par(1)*x)/par(1)**4
    end function df33
    
    pure real(kind(0d0)) function d2f33(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        d2f33 = exp(-par(1)*x)/par(1)**3
    end function d2f33
    
    pure real(kind(0d0)) function f34(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        f34 = x*exp(-par(1)*x)/par(1)**5
    end function f34

    pure real(kind(0d0)) function df34(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        df34 = exp(-par(1)*x)/par(1)**4 * (1d0/par(1) - x)
    end function df34
    
    pure real(kind(0d0)) function d2f34(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        d2f34 = exp(-par(1)*x)/par(1)**3 * (x - 2d0/par(1))
    end function d2f34
    
    pure real(kind(0d0)) function f35(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        f35 = x**2*exp(-par(1)*x)/par(1)**5
    end function f35

    pure real(kind(0d0)) function df35(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        df35 = exp(-par(1)*x)/par(1)**4 * (2d0*x/par(1) - x**2)
    end function df35
    
    pure real(kind(0d0)) function d2f35(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        d2f35 = exp(-par(1)*x)/par(1)**3 * (x**2 - 4d0*x/par(1) + 2d0/par(1)**2)
    end function d2f35
    
end module kernel_exp_3


!-----------------------------------------------------------------------
!> @author
!> Oliver T. Unke, University of Basel
!
!> @brief
!> Reciprocal power decay kernel with n = 2 and m = 0
!
!> @details
!> This module contains all functions necessary to evaluate the reciprocal power
!! decay kernel with n = 2 and m = 0. This corresponds to a 1/r decay for
!! values larger than the greatest point in the grid. This kernel is recommended
!! for modelling charge-charge interactions. The kernel is only
!! defined for values in the interval [0,infinity).
!-----------------------------------------------------------------------
module kernel_rp_2_0
    implicit none
    integer, parameter :: M2   = 2
    integer, parameter :: Npar = 0
    real(kind(0d0)), parameter :: p21 = 2d0  ,&
                                  p22 = -2d0/3d0
    contains   
         
    pure real(kind(0d0)) function k(x1,x2,par)
        implicit none
        real(kind(0d0)), intent(in) :: x1,x2
        real(kind(0d0))             :: xs,xl
        real(kind(0d0)), dimension(:), intent(in) :: par 
        !find larger/smaller of x1 and x2
        if(x1 <= x2) then
            xs = x1
            xl = x2
        else
            xs = x2
            xl = x1
        end if
        k = 2d0/xl - 2d0/3d0 * xs/xl**2
    end function k
    
    !direct kernel impelementation
    pure real(kind(0d0)) function dk(x1,x2,par)
        implicit none
        real(kind(0d0)), intent(in) :: x1,x2
        real(kind(0d0)), dimension(:), intent(in) :: par 
        !find larger/smaller of x1 and x2
        if(x1 <= x2) then
            dk = -2d0/(3d0*x2**2)
        else
            dk = 4d0/3d0 * x2/x1**3 - 2d0/x1**2
        end if
    end function dk
    
    !direct kernel impelementation
    pure real(kind(0d0)) function d2k(x1,x2,par)
        implicit none
        real(kind(0d0)), intent(in) :: x1,x2
        real(kind(0d0)), dimension(:), intent(in) :: par 
        !find larger/smaller of x1 and x2
        if(x1 <= x2) then
            d2k = 0d0
        else
            d2k = 4d0/x1**3 - 4d0*x2/x1**4
        end if
    end function d2k
        
    !f2 and f3 functions and their derivatives
    pure real(kind(0d0)) function f21(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        f21 = 1d0
    end function f21
    
    pure real(kind(0d0)) function df21(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        df21 = 0d0
    end function df21
    
    pure real(kind(0d0)) function d2f21(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        d2f21 = 0d0
    end function d2f21
    
    pure real(kind(0d0)) function f22(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        f22 = x
    end function f22
    
    pure real(kind(0d0)) function df22(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        df22 = 1d0
    end function df22
    
    pure real(kind(0d0)) function d2f22(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        d2f22 = 0d0
    end function d2f22
    
    pure real(kind(0d0)) function f31(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        f31 = 1d0/x
    end function f31
    
    pure real(kind(0d0)) function df31(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        df31 = -1d0/x**2
    end function df31
    
    pure real(kind(0d0)) function d2f31(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        d2f31 = 2d0/x**3
    end function d2f31
    
    pure real(kind(0d0)) function f32(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        f32 = 1d0/x**2
    end function f32

    pure real(kind(0d0)) function df32(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        df32 = -2d0/x**3
    end function df32
    
    pure real(kind(0d0)) function d2f32(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        d2f32 = 6d0/x**4
    end function d2f32
    
end module kernel_rp_2_0


!-----------------------------------------------------------------------
!> @author
!> Oliver T. Unke, University of Basel
!
!> @brief
!> Reciprocal power decay kernel with n = 3 and m = 0
!
!> @details
!> This module contains all functions necessary to evaluate the reciprocal power
!! decay kernel with n = 3 and m = 0. This corresponds to a 1/r decay for
!! values larger than the greatest point in the grid. This kernel is recommended
!! for modelling charge-charge interactions. The kernel is only
!! defined for values in the interval [0,infinity).
!-----------------------------------------------------------------------
module kernel_rp_3_0
    implicit none
    integer, parameter :: M2   = 3
    integer, parameter :: Npar = 0
    real(kind(0d0)), parameter :: p21 = 3d0,&
                                  p22 = -3d0/2d0,&
                                  p23 = 3d0/10d0
    contains        
    !direct kernel impelementation
    pure real(kind(0d0)) function k(x1,x2,par)
        implicit none
        real(kind(0d0)), intent(in) :: x1,x2
        real(kind(0d0))             :: xs,xl
        real(kind(0d0)), dimension(:), intent(in) :: par 
        !find larger/smaller of x1 and x2
        if(x1 <= x2) then
            xs = x1
            xl = x2
        else
            xs = x2
            xl = x1
        end if
        k = 3d0/xl - 3d0/2d0 * xs/xl**2 + 3d0/10d0 * xs**2/xl**3
    end function k
    
    !direct kernel impelementation
    pure real(kind(0d0)) function dk(x1,x2,par)
        implicit none
        real(kind(0d0)), intent(in) :: x1,x2
        real(kind(0d0)), dimension(:), intent(in) :: par 
        !find larger/smaller of x1 and x2
        if(x1 <= x2) then
            dk = 3d0/5d0 * x1/x2**3 - 3d0/(2d0*x2**2)
        else
            dk = -3d0/x1**2 + 3d0*x2/x1**3 - 9d0/10d0 * x2**2/x1**4
        end if
    end function dk
    
    !direct kernel impelementation
    pure real(kind(0d0)) function d2k(x1,x2,par)
        implicit none
        real(kind(0d0)), intent(in) :: x1,x2
        real(kind(0d0)), dimension(:), intent(in) :: par 
        !find larger/smaller of x1 and x2
        if(x1 <= x2) then
            d2k = 3d0/(5d0*x2**3)
        else
            d2k = 6d0/x1**3 - 9d0*x2/x1**4 + 18d0/5d0*x2**2/x1**5
        end if
    end function d2k
        
    !f2 and f3 functions and their derivatives
    pure real(kind(0d0)) function f21(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        f21 = 1d0
    end function f21
    
    pure real(kind(0d0)) function df21(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        df21 = 0d0
    end function df21
    
    pure real(kind(0d0)) function d2f21(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        d2f21 = 0d0
    end function d2f21
    
    pure real(kind(0d0)) function f22(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        f22 = x
    end function f22
    
    pure real(kind(0d0)) function df22(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        df22 = 1d0
    end function df22
    
    pure real(kind(0d0)) function d2f22(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        d2f22 = 0d0
    end function d2f22
    
    pure real(kind(0d0)) function f23(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        f23 = x**2
    end function f23
    
    pure real(kind(0d0)) function df23(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        df23 = 2d0*x
    end function df23
    
    pure real(kind(0d0)) function d2f23(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        d2f23 = 2d0
    end function d2f23
    
    pure real(kind(0d0)) function f31(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        f31 = 1d0/x
    end function f31
    
    pure real(kind(0d0)) function df31(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        df31 = -1d0/x**2
    end function df31
    
    pure real(kind(0d0)) function d2f31(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        d2f31 = 2d0/x**3
    end function d2f31
    
    pure real(kind(0d0)) function f32(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        f32 = 1d0/x**2
    end function f32

    pure real(kind(0d0)) function df32(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        df32 = -2d0/x**3
    end function df32
    
    pure real(kind(0d0)) function d2f32(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        d2f32 = 6d0/x**4
    end function d2f32
    
    pure real(kind(0d0)) function f33(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        f33 = 1d0/x**3
    end function f33

    pure real(kind(0d0)) function df33(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        df33 = -3d0/x**4
    end function df33
    
    pure real(kind(0d0)) function d2f33(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        d2f33 = 12d0/x**5
    end function d2f33
    
end module kernel_rp_3_0


!-----------------------------------------------------------------------
!> @author
!> Oliver T. Unke, University of Basel
!
!> @brief
!> Reciprocal power decay kernel with n = 2 and m = 1
!
!> @details
!> This module contains all functions necessary to evaluate the reciprocal power
!! decay kernel with n = 2 and m = 1. This corresponds to a 1/r^2 decay for
!! values larger than the greatest point in the grid. This kernel is recommended
!! for modelling charge-dipole interactions. The kernel is only
!! defined for values in the interval [0,infinity).
!-----------------------------------------------------------------------
module kernel_rp_2_1
    implicit none
    integer, parameter :: M2   = 2
    integer, parameter :: Npar = 0
    real(kind(0d0)), parameter :: p21 = 2d0/3d0,&
                                  p22 = -1d0/3d0
    contains        
    !direct kernel impelementation
    pure real(kind(0d0)) function k(x1,x2,par)
        implicit none
        real(kind(0d0)), intent(in) :: x1,x2
        real(kind(0d0))             :: xs,xl
        real(kind(0d0)), dimension(:), intent(in) :: par 
        !find larger/smaller of x1 and x2
        if(x1 <= x2) then
            xs = x1
            xl = x2
        else
            xs = x2
            xl = x1
        end if
        k = 2d0/(3d0*xl**2) - 1d0/3d0 * xs/xl**3
    end function k
    
    !direct kernel impelementation
    pure real(kind(0d0)) function dk(x1,x2,par)
        implicit none
        real(kind(0d0)), intent(in) :: x1,x2
        real(kind(0d0)), dimension(:), intent(in) :: par 
        !find larger/smaller of x1 and x2
        if(x1 <= x2) then
            dk = -1d0/(3d0*x2**3)
        else
            dk = x2/x1**4 - 4d0/(3d0*x1**3)
        end if
    end function dk
    
    !direct kernel impelementation
    pure real(kind(0d0)) function d2k(x1,x2,par)
        implicit none
        real(kind(0d0)), intent(in) :: x1,x2
        real(kind(0d0)), dimension(:), intent(in) :: par 
        !find larger/smaller of x1 and x2
        if(x1 <= x2) then
            d2k = 0d0
        else
            d2k = 4d0/x1**4 - 4d0*x2/x1**5
        end if
    end function d2k
        
    !f2 and f3 functions and their derivatives
    pure real(kind(0d0)) function f21(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        f21 = 1d0
    end function f21
    
    pure real(kind(0d0)) function df21(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        df21 = 0d0
    end function df21
    
    pure real(kind(0d0)) function d2f21(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        d2f21 = 0d0
    end function d2f21
    
    pure real(kind(0d0)) function f22(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        f22 = x
    end function f22
    
    pure real(kind(0d0)) function df22(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        df22 = 1d0
    end function df22
    
    pure real(kind(0d0)) function d2f22(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        d2f22 = 0d0
    end function d2f22
    
    pure real(kind(0d0)) function f31(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        f31 = 1d0/x**2
    end function f31
    
    pure real(kind(0d0)) function df31(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        df31 = -2d0/x**3
    end function df31
    
    pure real(kind(0d0)) function d2f31(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        d2f31 = 6d0/x**4
    end function d2f31
    
    pure real(kind(0d0)) function f32(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        f32 = 1d0/x**3
    end function f32

    pure real(kind(0d0)) function df32(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        df32 = -3d0/x**4
    end function df32
    
    pure real(kind(0d0)) function d2f32(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        d2f32 = 12d0/x**5
    end function d2f32
    
end module kernel_rp_2_1


!-----------------------------------------------------------------------
!> @author
!> Oliver T. Unke, University of Basel
!
!> @brief
!> Reciprocal power decay kernel with n = 3 and m = 1
!
!> @details
!> This module contains all functions necessary to evaluate the reciprocal power
!! decay kernel with n = 3 and m = 1. This corresponds to a 1/r^2 decay for
!! values larger than the greatest point in the grid. This kernel is recommended
!! for modelling charge-dipole interactions. The kernel is only
!! defined for values in the interval [0,infinity).
!-----------------------------------------------------------------------
module kernel_rp_3_1
    implicit none
    integer, parameter :: M2   = 3
    integer, parameter :: Npar = 0
    real(kind(0d0)), parameter :: p21 = 3d0/4d0,&
                                  p22 = -3d0/5d0,&
                                  p23 = 3d0/20d0
    contains        
    !direct kernel impelementation
    pure real(kind(0d0)) function k(x1,x2,par)
        implicit none
        real(kind(0d0)), intent(in) :: x1,x2
        real(kind(0d0))             :: xs,xl
        real(kind(0d0)), dimension(:), intent(in) :: par 
        !find larger/smaller of x1 and x2
        if(x1 <= x2) then
            xs = x1
            xl = x2
        else
            xs = x2
            xl = x1
        end if
        k = 3d0/(4d0*xl**2) - 3d0/5d0 * xs/xl**3 + 3d0/20d0 * xs**2/xl**4
    end function k
    
    !direct kernel impelementation
    pure real(kind(0d0)) function dk(x1,x2,par)
        implicit none
        real(kind(0d0)), intent(in) :: x1,x2
        real(kind(0d0)), dimension(:), intent(in) :: par 
        !find larger/smaller of x1 and x2
        if(x1 <= x2) then
            dk = 3d0/10d0 * x1/x2**4 - 3d0/(5d0*x2**3)
        else
            dk = -3d0/(2d0*x1**3) + 9d0/5d0 * x2/x1**4 - 3d0/5d0 * x2**2/x1**5
        end if
    end function dk
    
    !direct kernel impelementation
    pure real(kind(0d0)) function d2k(x1,x2,par)
        implicit none
        real(kind(0d0)), intent(in) :: x1,x2
        real(kind(0d0)), dimension(:), intent(in) :: par 
        !find larger/smaller of x1 and x2
        if(x1 <= x2) then
            d2k = 3d0/(10d0*x2**4)
        else
            d2k = 9d0/(2d0*x1**4) - 36d0/5d0*x2/x1**5 + 3d0*x2**2/x1**6
        end if
    end function d2k
        
    !f2 and f3 functions and their derivatives
    pure real(kind(0d0)) function f21(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        f21 = 1d0
    end function f21
    
    pure real(kind(0d0)) function df21(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        df21 = 0d0
    end function df21
    
    pure real(kind(0d0)) function d2f21(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        d2f21 = 0d0
    end function d2f21
    
    pure real(kind(0d0)) function f22(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        f22 = x
    end function f22
    
    pure real(kind(0d0)) function df22(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        df22 = 1d0
    end function df22
    
    pure real(kind(0d0)) function d2f22(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        d2f22 = 0d0
    end function d2f22
    
    pure real(kind(0d0)) function f23(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        f23 = x**2
    end function f23
    
    pure real(kind(0d0)) function df23(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        df23 = 2d0*x
    end function df23
    
    pure real(kind(0d0)) function d2f23(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        d2f23 = 2d0
    end function d2f23
    
    pure real(kind(0d0)) function f31(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        f31 = 1d0/x**2
    end function f31
    
    pure real(kind(0d0)) function df31(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        df31 = -2d0/x**3
    end function df31
    
    pure real(kind(0d0)) function d2f31(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        d2f31 = 6d0/x**4
    end function d2f31
    
    pure real(kind(0d0)) function f32(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        f32 = 1d0/x**3
    end function f32

    pure real(kind(0d0)) function df32(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        df32 = -3d0/x**4
    end function df32
    
    pure real(kind(0d0)) function d2f32(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        d2f32 = 12d0/x**5
    end function d2f32
    
    pure real(kind(0d0)) function f33(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        f33 = 1d0/x**4
    end function f33

    pure real(kind(0d0)) function df33(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        df33 = -4d0/x**5
    end function df33
    
    pure real(kind(0d0)) function d2f33(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        d2f33 = 20d0/x**6
    end function d2f33
    
end module kernel_rp_3_1


!-----------------------------------------------------------------------
!> @author
!> Oliver T. Unke, University of Basel
!
!> @brief
!> Reciprocal power decay kernel with n = 2 and m = 2
!
!> @details
!> This module contains all functions necessary to evaluate the reciprocal power
!! decay kernel with n = 2 and m = 2. This corresponds to a 1/r^3 decay for
!! values larger than the greatest point in the grid. This kernel is recommended
!! for modelling dipole-dipole interactions. The kernel is only
!! defined for values in the interval [0,infinity).
!-----------------------------------------------------------------------
module kernel_rp_2_2
    implicit none
    integer, parameter :: M2   = 2
    integer, parameter :: Npar = 0
    real(kind(0d0)), parameter :: p21 = 1d0/3d0,&
                                  p22 = -1d0/5d0
    contains        
    !direct kernel impelementation
    pure real(kind(0d0)) function k(x1,x2,par)
        implicit none
        real(kind(0d0)), intent(in) :: x1,x2
        real(kind(0d0))             :: xs,xl
        real(kind(0d0)), dimension(:), intent(in) :: par 
        !find larger/smaller of x1 and x2
        if(x1 <= x2) then
            xs = x1
            xl = x2
        else
            xs = x2
            xl = x1
        end if
        k = 1d0/(3d0*xl**3) - 1d0/5d0 * xs/xl**4
    end function k
    
    !direct kernel impelementation
    pure real(kind(0d0)) function dk(x1,x2,par)
        implicit none
        real(kind(0d0)), intent(in) :: x1,x2
        real(kind(0d0)), dimension(:), intent(in) :: par 
        !find larger/smaller of x1 and x2
        if(x1 <= x2) then
            dk = -1d0/(5d0*x2**4)
        else
            dk = 4d0/5d0*(x2/x1**5) - 1d0/x1**4
        end if
    end function dk
    
    !direct kernel impelementation
    pure real(kind(0d0)) function d2k(x1,x2,par)
        implicit none
        real(kind(0d0)), intent(in) :: x1,x2
        real(kind(0d0)), dimension(:), intent(in) :: par 
        !find larger/smaller of x1 and x2
        if(x1 <= x2) then
            d2k = 0d0
        else
            d2k = 4d0/x1**5 - 4d0*x2/x1**6
        end if
    end function d2k
        
    !f2 and f3 functions and their derivatives
    pure real(kind(0d0)) function f21(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        f21 = 1d0
    end function f21
    
    pure real(kind(0d0)) function df21(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        df21 = 0d0
    end function df21
    
    pure real(kind(0d0)) function d2f21(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        d2f21 = 0d0
    end function d2f21
    
    pure real(kind(0d0)) function f22(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        f22 = x
    end function f22
    
    pure real(kind(0d0)) function df22(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        df22 = 1d0
    end function df22
    
    pure real(kind(0d0)) function d2f22(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        d2f22 = 0d0
    end function d2f22
    
    pure real(kind(0d0)) function f31(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        f31 = 1d0/x**3
    end function f31
    
    pure real(kind(0d0)) function df31(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        df31 = -3d0/x**4
    end function df31
    
    pure real(kind(0d0)) function d2f31(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        d2f31 = 12d0/x**5
    end function d2f31
    
    pure real(kind(0d0)) function f32(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        f32 = 1d0/x**4
    end function f32

    pure real(kind(0d0)) function df32(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        df32 = -4d0/x**5
    end function df32
    
    pure real(kind(0d0)) function d2f32(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        d2f32 = 20d0/x**6
    end function d2f32
    
end module kernel_rp_2_2


!-----------------------------------------------------------------------
!> @author
!> Oliver T. Unke, University of Basel
!
!> @brief
!> Reciprocal power decay kernel with n = 3 and m = 2
!
!> @details
!> This module contains all functions necessary to evaluate the reciprocal power
!! decay kernel with n = 3 and m = 2. This corresponds to a 1/r^3 decay for
!! values larger than the greatest point in the grid. This kernel is recommended
!! for modelling dipole-dipole interactions. The kernel is only
!! defined for values in the interval [0,infinity).
!-----------------------------------------------------------------------
module kernel_rp_3_2
    implicit none
    integer, parameter :: M2   = 3
    integer, parameter :: Npar = 0
    real(kind(0d0)), parameter :: p21 = 3d0/10d0,&
                                  p22 = -3d0/10d0,&
                                  p23 = 3d0/35d0
    contains        
    !direct kernel impelementation
    pure real(kind(0d0)) function k(x1,x2,par)
        implicit none
        real(kind(0d0)), intent(in) :: x1,x2
        real(kind(0d0))             :: xs,xl
        real(kind(0d0)), dimension(:), intent(in) :: par 
        !find larger/smaller of x1 and x2
        if(x1 <= x2) then
            xs = x1
            xl = x2
        else
            xs = x2
            xl = x1
        end if
        k = 3d0/(10d0*xl**3) - 3d0/10d0 * xs/xl**4 + 3d0/35d0 * xs**2/xl**5
    end function k
    
    !direct kernel impelementation
    pure real(kind(0d0)) function dk(x1,x2,par)
        implicit none
        real(kind(0d0)), intent(in) :: x1,x2
        real(kind(0d0)), dimension(:), intent(in) :: par 
        !find larger/smaller of x1 and x2
        if(x1 <= x2) then
            dk = 6d0/35d0 * x1/x2**5 - 3d0/(10d0*x2**4)
        else
            dk = -9d0/(10d0*x1**4) + 6d0/5d0 * x2/x1**5 - 3d0/7d0 * x2**2/x1**6
        end if
    end function dk
    
    !direct kernel impelementation
    pure real(kind(0d0)) function d2k(x1,x2,par)
        implicit none
        real(kind(0d0)), intent(in) :: x1,x2
        real(kind(0d0)), dimension(:), intent(in) :: par 
        !find larger/smaller of x1 and x2
        if(x1 <= x2) then
            d2k = 6d0/(35d0*x2**5)
        else
            d2k = 18d0/(5d0*x1**5) - 6d0*x2/x1**6 + 18d0/7d0*x2**2/x1**7
        end if
    end function d2k
        
    !f2 and f3 functions and their derivatives
    pure real(kind(0d0)) function f21(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        f21 = 1d0
    end function f21
    
    pure real(kind(0d0)) function df21(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        df21 = 0d0
    end function df21
    
    pure real(kind(0d0)) function d2f21(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        d2f21 = 0d0
    end function d2f21
    
    pure real(kind(0d0)) function f22(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        f22 = x
    end function f22
    
    pure real(kind(0d0)) function df22(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        df22 = 1d0
    end function df22
    
    pure real(kind(0d0)) function d2f22(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        d2f22 = 0d0
    end function d2f22
    
    pure real(kind(0d0)) function f23(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        f23 = x**2
    end function f23
    
    pure real(kind(0d0)) function df23(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        df23 = 2d0*x
    end function df23
    
    pure real(kind(0d0)) function d2f23(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        d2f23 = 2d0
    end function d2f23
    
    pure real(kind(0d0)) function f31(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        f31 = 1d0/x**3
    end function f31
    
    pure real(kind(0d0)) function df31(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        df31 = -3d0/x**4
    end function df31
    
    pure real(kind(0d0)) function d2f31(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        d2f31 = 12d0/x**5
    end function d2f31
    
    pure real(kind(0d0)) function f32(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        f32 = 1d0/x**4
    end function f32

    pure real(kind(0d0)) function df32(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        df32 = -4d0/x**5
    end function df32
    
    pure real(kind(0d0)) function d2f32(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        d2f32 = 20d0/x**6
    end function d2f32
    
    pure real(kind(0d0)) function f33(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        f33 = 1d0/x**5
    end function f33

    pure real(kind(0d0)) function df33(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        df33 = -5d0/x**6
    end function df33
    
    pure real(kind(0d0)) function d2f33(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        d2f33 = 30d0/x**7
    end function d2f33
    
end module kernel_rp_3_2


!-----------------------------------------------------------------------
!> @author
!> Oliver T. Unke, University of Basel
!
!> @brief
!> Reciprocal power decay kernel with n = 2 and m = 3
!
!> @details
!> This module contains all functions necessary to evaluate the reciprocal power
!! decay kernel with n = 2 and m = 3. This corresponds to a 1/r^4 decay for
!! values larger than the greatest point in the grid. This kernel is recommended
!! for modelling charge-induced dipole interactions. The kernel is only
!! defined for values in the interval [0,infinity).
!-----------------------------------------------------------------------
module kernel_rp_2_3
    implicit none
    integer, parameter :: M2   = 2
    integer, parameter :: Npar = 0
    real(kind(0d0)), parameter :: p21 = 1d0/5d0  ,&
                                  p22 = -2d0/15d0
    contains        
    !direct kernel impelementation
    pure real(kind(0d0)) function k(x1,x2,par)
        implicit none
        real(kind(0d0)), intent(in) :: x1,x2
        real(kind(0d0))             :: xs,xl
        real(kind(0d0)), dimension(:), intent(in) :: par 
        !find larger/smaller of x1 and x2
        if(x1 <= x2) then
            xs = x1
            xl = x2
        else
            xs = x2
            xl = x1
        end if
        k = 1d0/(5d0*xl**4) - 2d0/15d0 * xs/xl**5
    end function k
    
    !direct kernel impelementation
    pure real(kind(0d0)) function dk(x1,x2,par)
        implicit none
        real(kind(0d0)), intent(in) :: x1,x2
        real(kind(0d0)), dimension(:), intent(in) :: par 
        !find larger/smaller of x1 and x2
        if(x1 <= x2) then
            dk = -2d0/(15d0*x2**5)
        else
            dk = 2d0/3d0*x2/x1**6 - 4d0/(5d0*x1**5)
        end if
    end function dk
    
    !direct kernel impelementation
    pure real(kind(0d0)) function d2k(x1,x2,par)
        implicit none
        real(kind(0d0)), intent(in) :: x1,x2
        real(kind(0d0)), dimension(:), intent(in) :: par 
        !find larger/smaller of x1 and x2
        if(x1 <= x2) then
            d2k = 0d0
        else
            d2k = 4d0/x1**6 - 4d0*x2/x1**7
        end if
    end function d2k
        
    !f2 and f3 functions and their derivatives
    pure real(kind(0d0)) function f21(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        f21 = 1d0
    end function f21
    
    pure real(kind(0d0)) function df21(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        df21 = 0d0
    end function df21
    
    pure real(kind(0d0)) function d2f21(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        d2f21 = 0d0
    end function d2f21
    
    pure real(kind(0d0)) function f22(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        f22 = x
    end function f22
    
    pure real(kind(0d0)) function df22(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        df22 = 1d0
    end function df22
    
    pure real(kind(0d0)) function d2f22(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        d2f22 = 0d0
    end function d2f22
    
    pure real(kind(0d0)) function f31(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        f31 = 1d0/x**4
    end function f31
    
    pure real(kind(0d0)) function df31(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        df31 = -4d0/x**5
    end function df31
    
    pure real(kind(0d0)) function d2f31(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        d2f31 = 20d0/x**6
    end function d2f31
    
    pure real(kind(0d0)) function f32(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        f32 = 1d0/x**5
    end function f32

    pure real(kind(0d0)) function df32(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        df32 = -5d0/x**6
    end function df32
    
    pure real(kind(0d0)) function d2f32(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        d2f32 = 30d0/x**7
    end function d2f32
    
end module kernel_rp_2_3


!-----------------------------------------------------------------------
!> @author
!> Oliver T. Unke, University of Basel
!
!> @brief
!> Reciprocal power decay kernel with n = 3 and m = 3
!
!> @details
!> This module contains all functions necessary to evaluate the reciprocal power
!! decay kernel with n = 3 and m = 3. This corresponds to a 1/r^4 decay for
!! values larger than the greatest point in the grid. This kernel is recommended
!! for modelling charge-induced dipole interactions. The kernel is only
!! defined for values in the interval [0,infinity).
!-----------------------------------------------------------------------
module kernel_rp_3_3
    implicit none
    integer, parameter :: M2   = 3
    integer, parameter :: Npar = 0
    real(kind(0d0)), parameter :: p21 = 3d0/20,&
                                  p22 = -6d0/35d0,&
                                  p23 = 3d0/56d0
    contains        
    !direct kernel impelementation
    pure real(kind(0d0)) function k(x1,x2,par)
        implicit none
        real(kind(0d0)), intent(in) :: x1,x2
        real(kind(0d0))             :: xs,xl
        real(kind(0d0)), dimension(:), intent(in) :: par 
        !find larger/smaller of x1 and x2
        if(x1 <= x2) then
            xs = x1
            xl = x2
        else
            xs = x2
            xl = x1
        end if
        k = 3d0/(20d0*xl**4) - 6d0/35d0 * xs/xl**5 + 3d0/56d0 * xs**2/xl**6
    end function k
    
    !direct kernel impelementation
    pure real(kind(0d0)) function dk(x1,x2,par)
        implicit none
        real(kind(0d0)), intent(in) :: x1,x2
        real(kind(0d0)), dimension(:), intent(in) :: par 
        !find larger/smaller of x1 and x2
        if(x1 <= x2) then
            dk = 3d0/28d0 * x1/x2**6 - 6d0/(35d0*x2**5)
        else
            dk = -3d0/(5d0*x1**5) + 6d0/7d0 * x2/x1**6 - 9d0/28d0 * x2**2/x1**7
        end if 
    end function dk
    
    !direct kernel impelementation
    pure real(kind(0d0)) function d2k(x1,x2,par)
        implicit none
        real(kind(0d0)), intent(in) :: x1,x2
        real(kind(0d0)), dimension(:), intent(in) :: par 
        !find larger/smaller of x1 and x2
        if(x1 <= x2) then
            d2k = 3d0/(28d0*x2**6)
        else
            d2k = 3d0/x1**6 - 36d0/7d0*x2/x1**7 + 9d0/4d0*x2**2/x1**8
        end if 
    end function d2k
        
    !f2 and f3 functions and their derivatives
    pure real(kind(0d0)) function f21(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        f21 = 1d0
    end function f21
    
    pure real(kind(0d0)) function df21(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        df21 = 0d0
    end function df21
    
    pure real(kind(0d0)) function d2f21(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        d2f21 = 0d0
    end function d2f21
    
    pure real(kind(0d0)) function f22(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        f22 = x
    end function f22
    
    pure real(kind(0d0)) function df22(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        df22 = 1d0
    end function df22
    
    pure real(kind(0d0)) function d2f22(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        d2f22 = 0d0
    end function d2f22
    
    pure real(kind(0d0)) function f23(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        f23 = x**2
    end function f23
    
    pure real(kind(0d0)) function df23(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        df23 = 2d0*x
    end function df23
    
    pure real(kind(0d0)) function d2f23(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        d2f23 = 2d0
    end function d2f23
    
    pure real(kind(0d0)) function f31(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        f31 = 1d0/x**4
    end function f31
    
    pure real(kind(0d0)) function df31(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        df31 = -4d0/x**5
    end function df31
    
    pure real(kind(0d0)) function d2f31(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        d2f31 = 20d0/x**6
    end function d2f31
    
    pure real(kind(0d0)) function f32(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        f32 = 1d0/x**5
    end function f32

    pure real(kind(0d0)) function df32(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        df32 = -5d0/x**6
    end function df32
    
    pure real(kind(0d0)) function d2f32(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        d2f32 = 30d0/x**7
    end function d2f32
    
    pure real(kind(0d0)) function f33(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        f33 = 1d0/x**6
    end function f33

    pure real(kind(0d0)) function df33(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        df33 = -6d0/x**7
    end function df33
    
    pure real(kind(0d0)) function d2f33(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        d2f33 = 42d0/x**8
    end function d2f33
    
end module kernel_rp_3_3


!-----------------------------------------------------------------------
!> @author
!> Oliver T. Unke, University of Basel
!
!> @brief
!> Reciprocal power decay kernel with n = 2 and m = 4
!
!> @details
!> This module contains all functions necessary to evaluate the reciprocal power
!! decay kernel with n = 2 and m = 4. This corresponds to a 1/r^5 decay for
!! values larger than the greatest point in the grid. The kernel is only
!! defined for values in the interval [0,infinity).
!-----------------------------------------------------------------------
module kernel_rp_2_4
    implicit none
    integer, parameter :: M2   = 2
    integer, parameter :: Npar = 0
    real(kind(0d0)), parameter :: p21 = 2d0/15d0,&
                                  p22 = -2d0/21d0
    contains        
    !direct kernel impelementation
    pure real(kind(0d0)) function k(x1,x2,par)
        implicit none
        real(kind(0d0)), intent(in) :: x1,x2
        real(kind(0d0))             :: xs,xl
        real(kind(0d0)), dimension(:), intent(in) :: par 
        !find larger/smaller of x1 and x2
        if(x1 <= x2) then
            xs = x1
            xl = x2
        else
            xs = x2
            xl = x1
        end if
        k = 2d0/(15d0*xl**5) - 2d0/21d0*xs/xl**6
    end function k
    
    !direct kernel impelementation
    pure real(kind(0d0)) function dk(x1,x2,par)
        implicit none
        real(kind(0d0)), intent(in) :: x1,x2
        real(kind(0d0)), dimension(:), intent(in) :: par 
        !find larger/smaller of x1 and x2
        if(x1 <= x2) then
            dk = -2d0/(21d0*x2**6)
        else
            dk = -2d0/(3d0*x1**6) + 4d0/7d0 * x2/x1**7
        end if
    end function dk
    
    !direct kernel impelementation
    pure real(kind(0d0)) function d2k(x1,x2,par)
        implicit none
        real(kind(0d0)), intent(in) :: x1,x2
        real(kind(0d0)), dimension(:), intent(in) :: par 
        !find larger/smaller of x1 and x2
        if(x1 <= x2) then
            d2k = 0d0
        else
            d2k = 4d0/x1**7 - 4d0*x2/x1**8
        end if
    end function d2k
        
    !f2 and f3 functions and their derivatives
    pure real(kind(0d0)) function f21(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        f21 = 1d0
    end function f21
    
    pure real(kind(0d0)) function df21(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        df21 = 0d0
    end function df21
    
    pure real(kind(0d0)) function d2f21(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        d2f21 = 0d0
    end function d2f21
    
    pure real(kind(0d0)) function f22(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        f22 = x
    end function f22
    
    pure real(kind(0d0)) function df22(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        df22 = 1d0
    end function df22
    
    pure real(kind(0d0)) function d2f22(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        d2f22 = 0d0
    end function d2f22
    
    pure real(kind(0d0)) function f31(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        f31 = 1d0/x**5
    end function f31
    
    pure real(kind(0d0)) function df31(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        df31 = -5d0/x**6
    end function df31
    
    pure real(kind(0d0)) function d2f31(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        d2f31 = 30d0/x**7
    end function d2f31
    
    pure real(kind(0d0)) function f32(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        f32 = 1d0/x**6
    end function f32

    pure real(kind(0d0)) function df32(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        df32 = -6d0/x**7
    end function df32
    
    pure real(kind(0d0)) function d2f32(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        d2f32 = 42d0/x**8
    end function d2f32
    
end module kernel_rp_2_4


!-----------------------------------------------------------------------
!> @author
!> Oliver T. Unke, University of Basel
!
!> @brief
!> Reciprocal power decay kernel with n = 3 and m = 4
!
!> @details
!> This module contains all functions necessary to evaluate the reciprocal power
!! decay kernel with n = 3 and m = 4. This corresponds to a 1/r^5 decay for
!! values larger than the greatest point in the grid. The kernel is only
!! defined for values in the interval [0,infinity).
!-----------------------------------------------------------------------
module kernel_rp_3_4
    implicit none
    integer, parameter :: M2   = 3
    integer, parameter :: Npar = 0
    real(kind(0d0)), parameter :: p21 = 3d0/35d0,&
                                  p22 = -3d0/28d0,&
                                  p23 = 1d0/28d0
    contains        
    !direct kernel impelementation
    pure real(kind(0d0)) function k(x1,x2,par)
        implicit none
        real(kind(0d0)), intent(in) :: x1,x2
        real(kind(0d0))             :: xs,xl
        real(kind(0d0)), dimension(:), intent(in) :: par 
        !find larger/smaller of x1 and x2
        if(x1 <= x2) then
            xs = x1
            xl = x2
        else
            xs = x2
            xl = x1
        end if
        k = 3d0/(35d0*xl**5) - 3d0/28d0 * xs/xl**6 + 1d0/28d0 * xs**2/xl**7
    end function k
    
    !direct kernel impelementation
    pure real(kind(0d0)) function dk(x1,x2,par)
        implicit none
        real(kind(0d0)), intent(in) :: x1,x2
        real(kind(0d0)), dimension(:), intent(in) :: par 
        !find larger/smaller of x1 and x2
        if(x1 <= x2) then
            dk = 1d0/14d0 * x1/x2**7 - 3d0/(28d0*x2**6)
        else
            dk = -3d0/(7d0*x1**6) + 9d0/14d0 * x2/x1**7 - 1d0/4d0 * x2**2/x1**8
        end if
    end function dk
    
    !direct kernel impelementation
    pure real(kind(0d0)) function d2k(x1,x2,par)
        implicit none
        real(kind(0d0)), intent(in) :: x1,x2
        real(kind(0d0)), dimension(:), intent(in) :: par 
        !find larger/smaller of x1 and x2
        if(x1 <= x2) then
            d2k = 1d0/(14d0*x2**7)
        else
            d2k = 18d0/(7d0*x1**7) - 9d0/2d0*x2/x1**8 + 2d0*x2**2/x1**9 
        end if
    end function d2k
        
    !f2 and f3 functions and their derivatives
    pure real(kind(0d0)) function f21(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        f21 = 1d0
    end function f21
    
    pure real(kind(0d0)) function df21(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        df21 = 0d0
    end function df21
    
    pure real(kind(0d0)) function d2f21(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        d2f21 = 0d0
    end function d2f21
    
    pure real(kind(0d0)) function f22(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        f22 = x
    end function f22
    
    pure real(kind(0d0)) function df22(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        df22 = 1d0
    end function df22
    
    pure real(kind(0d0)) function d2f22(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        d2f22 = 0d0
    end function d2f22
    
    pure real(kind(0d0)) function f23(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        f23 = x**2
    end function f23
    
    pure real(kind(0d0)) function df23(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        df23 = 2d0*x
    end function df23
    
    pure real(kind(0d0)) function d2f23(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        d2f23 = 2d0
    end function d2f23
    
    pure real(kind(0d0)) function f31(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        f31 = 1d0/x**5
    end function f31
    
    pure real(kind(0d0)) function df31(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        df31 = -5d0/x**6
    end function df31
    
    pure real(kind(0d0)) function d2f31(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        d2f31 = 30d0/x**7
    end function d2f31
    
    pure real(kind(0d0)) function f32(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        f32 = 1d0/x**6
    end function f32

    pure real(kind(0d0)) function df32(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        df32 = -6d0/x**7
    end function df32
    
    pure real(kind(0d0)) function d2f32(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        d2f32 = 42d0/x**8
    end function d2f32
    
    pure real(kind(0d0)) function f33(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        f33 = 1d0/x**7
    end function f33

    pure real(kind(0d0)) function df33(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        df33 = -7d0/x**8
    end function df33
    
    pure real(kind(0d0)) function d2f33(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        d2f33 = 56d0/x**9
    end function d2f33
    
end module kernel_rp_3_4


!-----------------------------------------------------------------------
!> @author
!> Oliver T. Unke, University of Basel
!
!> @brief
!> Reciprocal power decay kernel with n = 2 and m = 5
!
!> @details
!> This module contains all functions necessary to evaluate the reciprocal power
!! decay kernel with n = 2 and m = 5. This corresponds to a 1/r^6 decay for
!! values larger than the greatest point in the grid. This kernel is recommended
!! for modelling dipole-induced dipole and dispersion interactions. The kernel is only
!! defined for values in the interval [0,infinity).
!-----------------------------------------------------------------------
module kernel_rp_2_5
    implicit none 
    integer, parameter :: M2   = 2
    integer, parameter :: Npar = 0
    real(kind(0d0)), parameter :: p21 = 2d0/21d0  ,&
                                  p22 = -1d0/14d0
    contains        
    !direct kernel impelementation
    pure real(kind(0d0)) function k(x1,x2,par)
        implicit none
        real(kind(0d0)), intent(in) :: x1,x2
        real(kind(0d0))             :: xs,xl
        real(kind(0d0)), dimension(:), intent(in) :: par 
        !find larger/smaller of x1 and x2
        if(x1 <= x2) then
            xs = x1
            xl = x2
        else
            xs = x2
            xl = x1
        end if
        k = 2d0/(21d0 * xl**6) - 1d0/14d0 * xs/xl**7
    end function k
    
    !direct kernel impelementation
    pure real(kind(0d0)) function dk(x1,x2,par)
        implicit none
        real(kind(0d0)), intent(in) :: x1,x2
        real(kind(0d0)), dimension(:), intent(in) :: par 
        !find larger/smaller of x1 and x2
        if(x1 <= x2) then
            dk = -1d0/(14d0*x2**7)   
        else
            dk = 1d0*x2/(2d0*x1**8) - 4d0/(7d0*x1**7) 
        end if
    end function dk
    
    !direct kernel impelementation
    pure real(kind(0d0)) function d2k(x1,x2,par)
        implicit none
        real(kind(0d0)), intent(in) :: x1,x2
        real(kind(0d0)), dimension(:), intent(in) :: par 
        !find larger/smaller of x1 and x2
        if(x1 <= x2) then
            d2k = 0d0   
        else
            d2k = 4d0/x1**8 - 4d0*x2/x1**9
        end if
    end function d2k
        
    !f2 and f3 functions and their derivatives
    pure real(kind(0d0)) function f21(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        f21 = 1d0
    end function f21
    
    pure real(kind(0d0)) function df21(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        df21 = 0d0
    end function df21
    
    pure real(kind(0d0)) function d2f21(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        d2f21 = 0d0
    end function d2f21
    
    pure real(kind(0d0)) function f22(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        f22 = x
    end function f22
    
    pure real(kind(0d0)) function df22(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        df22 = 1d0
    end function df22
    
    pure real(kind(0d0)) function d2f22(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        d2f22 = 0d0
    end function d2f22
    
    pure real(kind(0d0)) function f31(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        f31 = 1d0/x**6
    end function f31
    
    pure real(kind(0d0)) function df31(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        df31 = -6d0/x**7
    end function df31
    
    pure real(kind(0d0)) function d2f31(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        d2f31 = 42d0/x**8
    end function d2f31
    
    pure real(kind(0d0)) function f32(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        f32 = 1d0/x**7
    end function f32

    pure real(kind(0d0)) function df32(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        df32 = -7d0/x**8
    end function df32
    
    pure real(kind(0d0)) function d2f32(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        d2f32 = 56d0/x**9
    end function d2f32
    
end module kernel_rp_2_5


!-----------------------------------------------------------------------
!> @author
!> Oliver T. Unke, University of Basel
!
!> @brief
!> Reciprocal power decay kernel with n = 3 and m = 5
!
!> @details
!> This module contains all functions necessary to evaluate the reciprocal power
!! decay kernel with n = 3 and m = 5. This corresponds to a 1/r^6 decay for
!! values larger than the greatest point in the grid. This kernel is recommended
!! for modelling dipole-induced dipole and dispersion interactions. The kernel is only
!! defined for values in the interval [0,infinity).
!-----------------------------------------------------------------------
module kernel_rp_3_5
    implicit none
    integer, parameter :: M2   = 3
    integer, parameter :: Npar = 0
    real(kind(0d0)), parameter :: p21 = 3d0/56d0,&
                                  p22 = -1d0/14d0,&
                                  p23 = 1d0/40d0
    contains        
    !direct kernel impelementation
    pure real(kind(0d0)) function k(x1,x2,par)
        implicit none
        real(kind(0d0)), intent(in) :: x1,x2
        real(kind(0d0))             :: xs,xl
        real(kind(0d0)), dimension(:), intent(in) :: par 
        !find larger/smaller of x1 and x2
        if(x1 <= x2) then
            xs = x1
            xl = x2
        else
            xs = x2
            xl = x1
        end if
        k = 3d0/(56d0*xl**6) - 1d0/14d0 * xs/xl**7 + 1d0/40d0 * xs**2/xl**8
    end function k
    
    !direct kernel impelementation
    pure real(kind(0d0)) function dk(x1,x2,par)
        implicit none
        real(kind(0d0)), intent(in) :: x1,x2
        real(kind(0d0)), dimension(:), intent(in) :: par 
        !find larger/smaller of x1 and x2
        if(x1 <= x2) then
            dk = 1d0/20d0 * x1/x2**8 - 1d0/(14d0*x2**7)
        else
            dk = -1d0/5d0 * x2**2/x1**9 + 0.5d0 * x2/x1**8 - 9d0/(28d0*x1**7)
        end if
    end function dk
    
    !direct kernel impelementation
    pure real(kind(0d0)) function d2k(x1,x2,par)
        implicit none
        real(kind(0d0)), intent(in) :: x1,x2
        real(kind(0d0)), dimension(:), intent(in) :: par 
        !find larger/smaller of x1 and x2
        if(x1 <= x2) then
            d2k = 1d0/(20d0*x2**8)
        else
            d2k = 9d0/(4d0*x1**8) - 4d0*x2/x1**9 + 9d0/5d0*x2**2/x1**10
        end if
    end function d2k
        
    !f2 and f3 functions and their derivatives
    pure real(kind(0d0)) function f21(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        f21 = 1d0
    end function f21
    
    pure real(kind(0d0)) function df21(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        df21 = 0d0
    end function df21
    
    pure real(kind(0d0)) function d2f21(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        d2f21 = 0d0
    end function d2f21
    
    pure real(kind(0d0)) function f22(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        f22 = x
    end function f22
    
    pure real(kind(0d0)) function df22(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        df22 = 1d0
    end function df22
    
    pure real(kind(0d0)) function d2f22(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        d2f22 = 0d0
    end function d2f22
    
    pure real(kind(0d0)) function f23(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        f23 = x**2
    end function f23
    
    pure real(kind(0d0)) function df23(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        df23 = 2d0*x
    end function df23
    
    pure real(kind(0d0)) function d2f23(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        d2f23 = 2d0
    end function d2f23
    
    pure real(kind(0d0)) function f31(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        f31 = 1d0/x**6
    end function f31
    
    pure real(kind(0d0)) function df31(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        df31 = -6d0/x**7
    end function df31
    
    pure real(kind(0d0)) function d2f31(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        d2f31 = 42d0/x**8
    end function d2f31
    
    pure real(kind(0d0)) function f32(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        f32 = 1d0/x**7
    end function f32

    pure real(kind(0d0)) function df32(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        df32 = -7d0/x**8
    end function df32
    
    pure real(kind(0d0)) function d2f32(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        d2f32 = 56d0/x**9
    end function d2f32
    
    pure real(kind(0d0)) function f33(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        f33 = 1d0/x**8
    end function f33

    pure real(kind(0d0)) function df33(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        df33 = -8d0/x**9
    end function df33
    
    pure real(kind(0d0)) function d2f33(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        d2f33 = 72d0/x**10
    end function d2f33
    
end module kernel_rp_3_5


!-----------------------------------------------------------------------
!> @author
!> Oliver T. Unke, University of Basel
!
!> @brief
!> Reciprocal power decay kernel with n = 2 and m = 6
!
!> @details
!> This module contains all functions necessary to evaluate the reciprocal power
!! decay kernel with n = 2 and m = 6. This corresponds to a 1/r^7 decay for
!! values larger than the greatest point in the grid. The kernel is only
!! defined for values in the interval [0,infinity).
!-----------------------------------------------------------------------
module kernel_rp_2_6
    implicit none
    integer, parameter :: M2   = 2
    integer, parameter :: Npar = 0
    real(kind(0d0)), parameter :: p21 = 1d0/14d0,&
                                  p22 = -1d0/18d0
    contains        
    !direct kernel impelementation
    pure real(kind(0d0)) function k(x1,x2,par)
        implicit none
        real(kind(0d0)), intent(in) :: x1,x2
        real(kind(0d0))             :: xs,xl
        real(kind(0d0)), dimension(:), intent(in) :: par 
        !find larger/smaller of x1 and x2
        if(x1 <= x2) then
            xs = x1
            xl = x2
        else
            xs = x2
            xl = x1
        end if
        k = 1d0/(14d0*xl**7) - 1d0/18d0 * xs/xl**8
    end function k
    
    !direct kernel impelementation
    pure real(kind(0d0)) function dk(x1,x2,par)
        implicit none
        real(kind(0d0)), intent(in) :: x1,x2
        real(kind(0d0)), dimension(:), intent(in) :: par 
        !find larger/smaller of x1 and x2
        if(x1 <= x2) then
            dk = -1d0/(18d0*x2**8)
        else
            dk = 4d0/9d0 * x2/x1**9 - 1d0/(2d0*x1**8)
        end if
    end function dk
    
    !direct kernel impelementation
    pure real(kind(0d0)) function d2k(x1,x2,par)
        implicit none
        real(kind(0d0)), intent(in) :: x1,x2
        real(kind(0d0)), dimension(:), intent(in) :: par 
        !find larger/smaller of x1 and x2
        if(x1 <= x2) then
            d2k = 0d0
        else
            d2k = 4d0/x1**9 - 4d0*x2/x1**10
        end if
    end function d2k
        
    !f2 and f3 functions and their derivatives
    pure real(kind(0d0)) function f21(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        f21 = 1d0
    end function f21
    
    pure real(kind(0d0)) function df21(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        df21 = 0d0
    end function df21
    
    pure real(kind(0d0)) function d2f21(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        d2f21 = 0d0
    end function d2f21
    
    pure real(kind(0d0)) function f22(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        f22 = x
    end function f22
    
    pure real(kind(0d0)) function df22(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        df22 = 1d0
    end function df22
    
    pure real(kind(0d0)) function d2f22(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        d2f22 = 0d0
    end function d2f22
    
    pure real(kind(0d0)) function f31(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        f31 = 1d0/x**7
    end function f31
    
    pure real(kind(0d0)) function df31(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        df31 = -7d0/x**8
    end function df31
    
    pure real(kind(0d0)) function d2f31(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        d2f31 = 56d0/x**9
    end function d2f31
    
    pure real(kind(0d0)) function f32(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        f32 = 1d0/x**8
    end function f32

    pure real(kind(0d0)) function df32(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        df32 = -8d0/x**9
    end function df32
    
    pure real(kind(0d0)) function d2f32(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        d2f32 = 72d0/x**10
    end function d2f32
    
end module kernel_rp_2_6


!-----------------------------------------------------------------------
!> @author
!> Oliver T. Unke, University of Basel
!
!> @brief
!> Reciprocal power decay kernel with n = 3 and m = 6
!
!> @details
!> This module contains all functions necessary to evaluate the reciprocal power
!! decay kernel with n = 3 and m = 6. This corresponds to a 1/r^7 decay for
!! values larger than the greatest point in the grid. The kernel is only
!! defined for values in the interval [0,infinity).
!-----------------------------------------------------------------------
module kernel_rp_3_6
    implicit none
    integer, parameter :: M2   = 3
    integer, parameter :: Npar = 0
    real(kind(0d0)), parameter :: p21 = 1d0/28d0,&
                                  p22 = -1d0/20d0,&
                                  p23 = 1d0/55d0
    contains        
    !direct kernel impelementation
    pure real(kind(0d0)) function k(x1,x2,par)
        implicit none
        real(kind(0d0)), intent(in) :: x1,x2
        real(kind(0d0))             :: xs,xl
        real(kind(0d0)), dimension(:), intent(in) :: par 
        !find larger/smaller of x1 and x2
        if(x1 <= x2) then
            xs = x1
            xl = x2
        else
            xs = x2
            xl = x1
        end if
        k = 1d0/(28d0*xl**7) - 1d0/20d0 * xs/xl**8 + 1d0/55d0 * xs**2/xl**9
    end function k
    
    !direct kernel impelementation
    pure real(kind(0d0)) function dk(x1,x2,par)
        implicit none
        real(kind(0d0)), intent(in) :: x1,x2
        real(kind(0d0)), dimension(:), intent(in) :: par 
        !find larger/smaller of x1 and x2
        if(x1 <= x2) then
            dk = 2d0/55d0 * x1/x2**9 - 1d0/(20d0*x2**8)
        else
            dk = -1d0/(4d0*x1**8) + 2d0/5d0 * x2/x1**9 - 9d0/55d0 * x2**2/x1**10
        end if
    end function dk
    
    !direct kernel impelementation
    pure real(kind(0d0)) function d2k(x1,x2,par)
        implicit none
        real(kind(0d0)), intent(in) :: x1,x2
        real(kind(0d0)), dimension(:), intent(in) :: par 
        !find larger/smaller of x1 and x2
        if(x1 <= x2) then
            d2k = 2d0/(55d0*x2**9)
        else
            d2k = 2d0/x1**9 - 18d0/5d0*x2/x1**10 + 18d0/11d0*x2**2/x1**11
        end if
    end function d2k
        
    !f2 and f3 functions and their derivatives
    pure real(kind(0d0)) function f21(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        f21 = 1d0
    end function f21
    
    pure real(kind(0d0)) function df21(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        df21 = 0d0
    end function df21
    
    pure real(kind(0d0)) function d2f21(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        d2f21 = 0d0
    end function d2f21
    
    pure real(kind(0d0)) function f22(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        f22 = x
    end function f22
    
    pure real(kind(0d0)) function df22(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        df22 = 1d0
    end function df22
    
    pure real(kind(0d0)) function d2f22(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        d2f22 = 0d0
    end function d2f22
    
    pure real(kind(0d0)) function f23(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        f23 = x**2
    end function f23
    
    pure real(kind(0d0)) function df23(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        df23 = 2d0*x
    end function df23
    
    pure real(kind(0d0)) function d2f23(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        d2f23 = 2d0
    end function d2f23
    
    pure real(kind(0d0)) function f31(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        f31 = 1d0/x**7
    end function f31
    
    pure real(kind(0d0)) function df31(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        df31 = -7d0/x**8
    end function df31
    
    pure real(kind(0d0)) function d2f31(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        d2f31 = 56d0/x**9
    end function d2f31
    
    pure real(kind(0d0)) function f32(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        f32 = 1d0/x**8
    end function f32

    pure real(kind(0d0)) function df32(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        df32 = -8d0/x**9
    end function df32
    
    pure real(kind(0d0)) function d2f32(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        d2f32 = 72d0/x**10
    end function d2f32
    
    pure real(kind(0d0)) function f33(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        f33 = 1d0/x**9
    end function f33

    pure real(kind(0d0)) function df33(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        df33 = -9d0/x**10
    end function df33
    
    pure real(kind(0d0)) function d2f33(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        d2f33 = 90d0/x**11
    end function d2f33
    
end module kernel_rp_3_6


!-----------------------------------------------------------------------
!> @author
!> Oliver T. Unke, University of Basel
!
!> @brief
!> Taylor spline kernel with n = 2
!
!> @details
!> This module contains all functions necessary to evaluate the Taylor spline kernel
!! with n = 2. Note that this kernel is only defined for values in the interval [0,1].
!! It is useful for example in describing angular coordinates, as long as a new coordinate
!! is introduced which scales the angular coordinate to the interval [0,1]. 
!! An example for such a coordinate would be y = (1d0-cos(alpha))/2d0.
!! Note that by clever choice of y, it is possible to capture the symmetry inherent in the system.
!! The kernel is also applicable to general machine learning problems to interpolate any arbitrary 
!! function defined in a finite interval.
!-----------------------------------------------------------------------
module kernel_ts_2
    implicit none
    integer, parameter :: M2   = 3
    integer, parameter :: Npar = 0
    real(kind(0d0)), parameter :: p21 = -2d0/3d0,&
                                  p22 = 1d0,&
                                  p23 = 2d0
    contains        
    !direct kernel impelementation
    pure real(kind(0d0)) function k(x1,x2,par)
        implicit none
        real(kind(0d0)), intent(in) :: x1,x2
        real(kind(0d0))             :: xs,xl
        real(kind(0d0)), dimension(:), intent(in) :: par 
        !find larger/smaller of x1 and x2
        if(x1 <= x2) then
            xs = x1
            xl = x2
        else
            xs = x2
            xl = x1
        end if
        k = 1d0 + xs*xl + 2d0*xs**2*xl - 2d0/3d0*xs**3
    end function k
    
    !direct kernel impelementation
    pure real(kind(0d0)) function dk(x1,x2,par)
        implicit none
        real(kind(0d0)), intent(in) :: x1,x2
        real(kind(0d0)), dimension(:), intent(in) :: par 
        !find larger/smaller of x1 and x2
        if(x1 <= x2) then
            dk = x2 + 4d0*x2*x1 - 2d0*x1**2
        else
            dk = x2*(2d0*x2 + 1d0)
        end if
    end function dk
    
    !direct kernel impelementation
    pure real(kind(0d0)) function d2k(x1,x2,par)
        implicit none
        real(kind(0d0)), intent(in) :: x1,x2
        real(kind(0d0)), dimension(:), intent(in) :: par 
        !find larger/smaller of x1 and x2
        if(x1 <= x2) then
            d2k = 4d0*x2 - 4d0*x1
        else
            d2k = 0d0
        end if
    end function d2k
        
    !f2 and f3 functions and their derivatives
    pure real(kind(0d0)) function f21(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        f21 = x**3 - 3d0/2d0
    end function f21
    
    pure real(kind(0d0)) function df21(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        df21 = 3d0*x**2
    end function df21
    
    pure real(kind(0d0)) function d2f21(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        d2f21 = 6d0*x
    end function d2f21
    
    pure real(kind(0d0)) function f22(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        f22 = x
    end function f22
    
    pure real(kind(0d0)) function df22(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        df22 = 1d0
    end function df22
    
    pure real(kind(0d0)) function d2f22(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        d2f22 = 0d0
    end function d2f22
    
    pure real(kind(0d0)) function f23(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        f23 = x**2
    end function f23
    
    pure real(kind(0d0)) function df23(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        df23 = 2d0*x
    end function df23
    
    pure real(kind(0d0)) function d2f23(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        d2f23 = 2d0
    end function d2f23
    
    pure real(kind(0d0)) function f31(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        f31 = 1d0
    end function f31
    
    pure real(kind(0d0)) function df31(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        df31 = 0d0
    end function df31
    
    pure real(kind(0d0)) function d2f31(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        d2f31 = 0d0
    end function d2f31
    
    pure real(kind(0d0)) function f32(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        f32 = x
    end function f32

    pure real(kind(0d0)) function df32(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        df32 = 1d0
    end function df32
    
    pure real(kind(0d0)) function d2f32(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        d2f32 = 0d0
    end function d2f32
    
    pure real(kind(0d0)) function f33(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        f33 = x
    end function f33

    pure real(kind(0d0)) function df33(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        df33 = 1d0
    end function df33
    
    pure real(kind(0d0)) function d2f33(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        d2f33 = 0d0
    end function d2f33
    
end module kernel_ts_2


!-----------------------------------------------------------------------
!> @author
!> Oliver T. Unke, University of Basel
!
!> @brief
!> Taylor spline kernel with n = 3
!
!> @details
!> This module contains all functions necessary to evaluate the Taylor spline kernel
!! with n = 3. Note that this kernel is only defined for values in the interval [0,1].
!! It is useful for example in describing angular coordinates, as long as a new coordinate
!! is introduced which scales the angular coordinate to the interval [0,1].
!! An example for such a coordinate would be y = (1d0-cos(alpha))/2d0.
!! Note that by clever choice of y, it is possible to capture the symmetry inherent in the system.
!! The kernel is also applicable to general machine learning problems to interpolate any arbitrary 
!! function defined in a finite interval.
!-----------------------------------------------------------------------
module kernel_ts_3
    implicit none
    integer, parameter :: M2   = 5
    integer, parameter :: Npar = 0
    real(kind(0d0)), parameter :: p21 = 3d0/10d0,&
                                  p22 = -3d0/2d0,&
                                  p23 = 3d0,&
                                  p24 = 1d0,&
                                  p25 = 1d0
    contains        
    !direct kernel impelementation
    pure real(kind(0d0)) function k(x1,x2,par)
        implicit none
        real(kind(0d0)), intent(in) :: x1,x2
        real(kind(0d0))             :: xs,xl
        real(kind(0d0)), dimension(:), intent(in) :: par 
        !find larger/smaller of x1 and x2
        if(x1 <= x2) then
            xs = x1
            xl = x2
        else
            xs = x2
            xl = x1
        end if
        k = 3d0/10d0*xs**5 - 3d0/2d0*xl*xs**4 + 3d0*xl**2*xs**3 &
            + xs*xl + (xs*xl)**2 + 1d0
    end function k
    
    !direct kernel impelementation
    pure real(kind(0d0)) function dk(x1,x2,par)
        implicit none
        real(kind(0d0)), intent(in) :: x1,x2
        real(kind(0d0)), dimension(:), intent(in) :: par 
        !find larger/smaller of x1 and x2
        if(x1 <= x2) then
            dk = 2d0*x1*x2**2 + x2 + 9d0 * (x1*x2)**2 &
                 - 6d0*x1**3*x2 + 3d0/2d0 * x1**4
        else
            dk = 0.5d0*x2 * (12d0*x1*x2**2 - 3d0*x2**3 &
                            + 4d0*x1*x2 + 2d0)
        end if
    end function dk
    
    !direct kernel impelementation
    pure real(kind(0d0)) function d2k(x1,x2,par)
        implicit none
        real(kind(0d0)), intent(in) :: x1,x2
        real(kind(0d0)), dimension(:), intent(in) :: par 
        !find larger/smaller of x1 and x2
        if(x1 <= x2) then
            d2k = 6d0*x1**3 - 18d0*x1**2*x2 + 18d0*x1*x2**2 + 2d0*x2**2
        else
            d2k = 2d0*x2**2*(3d0*x2+1d0)
        end if
    end function d2k
        
    !f2 and f3 functions and their derivatives
    pure real(kind(0d0)) function f21(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        f21 = x**5 + 10d0/3d0 
    end function f21
    
    pure real(kind(0d0)) function df21(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        df21 = 5d0*x**4 
    end function df21
    
    pure real(kind(0d0)) function d2f21(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        d2f21 = 20d0*x**3
    end function d2f21
    
    pure real(kind(0d0)) function f22(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        f22 = x**4
    end function f22
    
    pure real(kind(0d0)) function df22(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        df22 = 4d0*x**3
    end function df22
    
    pure real(kind(0d0)) function d2f22(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        d2f22 = 12d0*x**2
    end function d2f22
    
    pure real(kind(0d0)) function f23(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        f23 = x**3
    end function f23
    
    pure real(kind(0d0)) function df23(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        df23 = 3d0*x**2
    end function df23
    
    pure real(kind(0d0)) function d2f23(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        d2f23 = 6d0*x
    end function d2f23
    
    pure real(kind(0d0)) function f24(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        f24 = x
    end function f24
    
    pure real(kind(0d0)) function df24(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        df24 = 1d0
    end function df24
    
    pure real(kind(0d0)) function d2f24(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        d2f24 = 0d0
    end function d2f24
    
    pure real(kind(0d0)) function f25(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        f25 = x**2
    end function f25
    
    pure real(kind(0d0)) function df25(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        df25 = 2d0*x
    end function df25
    
    pure real(kind(0d0)) function d2f25(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        d2f25 = 2d0
    end function d2f25
    
    pure real(kind(0d0)) function f31(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        f31 = 1d0
    end function f31
    
    pure real(kind(0d0)) function df31(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        df31 = 0d0
    end function df31
    
    pure real(kind(0d0)) function d2f31(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        d2f31 = 0d0
    end function d2f31
    
    pure real(kind(0d0)) function f32(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        f32 = x
    end function f32

    pure real(kind(0d0)) function df32(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        df32 = 1d0
    end function df32
    
    pure real(kind(0d0)) function d2f32(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        d2f32 = 0d0
    end function d2f32
    
    pure real(kind(0d0)) function f33(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        f33 = x**2
    end function f33

    pure real(kind(0d0)) function df33(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        df33 = 2d0*x
    end function df33
    
    pure real(kind(0d0)) function d2f33(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        d2f33 = 2d0
    end function d2f33
    
    pure real(kind(0d0)) function f34(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        f34 = x
    end function f34

    pure real(kind(0d0)) function df34(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        df34 = 1d0
    end function df34
    
    pure real(kind(0d0)) function d2f34(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        d2f34 = 0d0
    end function d2f34
    
    pure real(kind(0d0)) function f35(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        f35 = x**2
    end function f35

    pure real(kind(0d0)) function df35(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        df35 = 2d0*x
    end function df35
    
    pure real(kind(0d0)) function d2f35(x,par)
        implicit none
        real(kind(0d0)), intent(in) :: x
        real(kind(0d0)), dimension(:), intent(in) :: par 
        d2f35 = 2d0
    end function d2f35
    
end module kernel_ts_3


!-----------------------------------------------------------------------
!> @author
!> Oliver T. Unke, University of Basel
!
!> @brief
!> Some routines needed for solving matrix equations.
!
!> @details
!> This module contains all methods which are necessary to solve a matrix
!! equation by Cholesky decomposition. Namely, it contains a method for 
!! performing the decomposition and a method to solve an equation 
!! using forward and backward substitution.
!-----------------------------------------------------------------------
module LinearAlgebra
implicit none

private
public :: cholesky_decomposition, cholesky_solve, cholesky_inverse, &
          forward_substitution, backward_substitution

contains
!-----------------------------------------------------------------------
!> @author
!> Oliver T. Unke, University of Basel
!
!> @brief
!> Throws error messages if problems are encountered
!-----------------------------------------------------------------------
subroutine throw_error(proc,errormsg)
    implicit none
    character(len=*), intent(in) :: proc, errormsg
    write(*,*) "ERROR in module LinearAlgebra: "//proc//": "//errormsg
    stop
end subroutine throw_error

!-----------------------------------------------------------------------
!> @author
!> Oliver T. Unke, University of Basel
!
!> @brief
!> Solves L*y = b for y by forward substitution
!
!> @details
!> Here, L is a lower triangular matrix.
!-----------------------------------------------------------------------
subroutine forward_substitution(L,y,b)
    implicit none
    real(kind(0d0)), dimension(:,:), intent(in)  :: L
    real(kind(0d0)), dimension(:),   intent(out) :: y
    real(kind(0d0)), dimension(:),   intent(in)  :: b    
    real(kind(0d0)) :: s
    integer :: i,k,n
    
    n = size(L,dim=1)
    if(n /= size(y,dim=1) .or. n /= size(b,dim=1)) then
        call throw_error("forward_substitution","Matrix L is not the appropriate size for y and/or b!")
    end if
    
    !forward substitution
    y(1) = b(1)/L(1,1)
    do i = 2,n,1
        s = 0d0
        do k = 1,i-1,1
            s = s + L(i,k)*y(k)
        end do
        y(i) = (b(i)-s)/L(i,i)  
    end do
    return
end subroutine forward_substitution

!-----------------------------------------------------------------------
!> @author
!> Oliver T. Unke, University of Basel
!
!> @brief
!> Solves U*x = y for x by backward substitution
!
!> @details
!> Here, U is an upper triangular matrix.
!-----------------------------------------------------------------------
subroutine backward_substitution(U,x,y)
    implicit none
    real(kind(0d0)), dimension(:,:), intent(in)  :: U
    real(kind(0d0)), dimension(:),   intent(out) :: x
    real(kind(0d0)), dimension(:),   intent(in)  :: y    
    real(kind(0d0)) :: s
    integer :: i,k,n
    
    n = size(U,dim=1)
    if(n /= size(x,dim=1) .or. n /= size(y,dim=1)) then
        call throw_error("forward_substitution","Matrix U is not the appropriate size for x and/or y!")
    end if
    
    !back substitution
    x(n) = y(n)/U(n,n)
    do i = n-1,1,-1
        s = 0d0
        do k = i+1,n,1
            s = s + U(i,k)*x(k) 
        end do
        x(i) = (y(i)-s)/U(i,i)
    end do  
    return
end subroutine backward_substitution

!-----------------------------------------------------------------------
!> @author
!> Oliver T. Unke, University of Basel
!
!> @brief
!> does the Cholesky decomposition of a positive-definite symmetric square matrix M
!
!> @details
!> Note: Only the lower triagonal part of M is needed as input and only that 
!! part of the matrix is modified.
!-----------------------------------------------------------------------
subroutine cholesky_decomposition(M)
    implicit none
    real(kind(0d0)), dimension(:,:), intent(out) :: M
    integer :: i
    
    if(size(M,dim=1) /= size(M,dim=2)) then
        call throw_error("cholesky_decomposition","Matrix M is not a square matrix!")
    end if
    
    !do the Cholesky decomposition 
    do i = 1,size(M,dim=1)
        !diagonal component
        M(i,i) = M(i,i) - dot_product(M(i,1:i-1),M(i,1:i-1))
        if(M(i,i) > 0d0) then
            M(i,i) = sqrt(M(i,i))
        else
            call throw_error("cholesky_decomposition","Matrix M is not positive definite! "//&
                             "If you are sure M is correct, maybe the matrix is ill conditioned. "//&
                             "Try to use the Tikhonov regularization procedure.")
        end if
        
        !off-diagonal component
        if(i < size(M,dim=1)) M(i+1:size(M,dim=1),i) = (M(i+1:size(M,dim=1),i) - &
            matmul(M(i+1:size(M,dim=1),1:i-1),M(i,1:i-1))) / M(i,i)
    end do
    return 
end subroutine cholesky_decomposition

!-----------------------------------------------------------------------
!> @author
!> Oliver T. Unke, University of Basel
!
!> @brief
!> solves the matrix equation M*x = b for x
!
!> @details
!> The input matrix needs to be the lower tridiagonal matrix, which is obtained by
!! Cholesky decomposition. The solution is the obtained by forward and backward
!! substitution.
!-----------------------------------------------------------------------
subroutine cholesky_solve(L,x,b)
    implicit none
    real(kind(0d0)), dimension(:,:), intent(in)  :: L
    real(kind(0d0)), dimension(:),   intent(out) :: x
    real(kind(0d0)), dimension(:),   intent(in)  :: b
    real(kind(0d0)), dimension(size(x)) :: y
    integer :: n
    
    n = size(L,dim=1)
    if(n /= size(x,dim=1) .or. n /= size(b,dim=1)) then
        call throw_error("cholesky_solve","Matrix L is not the appropriate size for x and/or b!")
    end if
    
    call forward_substitution(L,y,b)
    call backward_substitution(transpose(L),x,y)
    return
end subroutine cholesky_solve

!-----------------------------------------------------------------------
!> @author
!> Oliver T. Unke, University of Basel
!
!> @brief
!> computes the inverse of a matrix
!
!> @details
!> The input matrix needs to be the lower tridiagonal matrix, which is obtained by
!! Cholesky decomposition. 
!-----------------------------------------------------------------------
subroutine cholesky_inverse(L,invL)
    implicit none
    real(kind(0d0)), dimension(:,:), intent(in)  :: L
    real(kind(0d0)), dimension(:,:), intent(out) :: invL
    real(kind(0d0)), dimension(size(L,dim=1))    :: b 
    integer :: i,n    
    n = size(L,dim=1)
    if(n /= size(invL,dim=1) .or. n /= size(invL,dim=2) .or. n /= size(L,dim=2)) then
        call throw_error("cholesky_inverse","Matrix L/invL is not the appropriate size!")
    end if    
    do i = 1,n
        b    = 0d0
        b(i) = 1d0
        call cholesky_solve(L,invL(i,:),b)     
    end do   
    return
end subroutine cholesky_inverse

end module LinearAlgebra


!-----------------------------------------------------------------------
!> @author
!> Oliver T. Unke, University of Basel
!
!> @brief
!> Contains all the definitions for the 1-dimensional kernel functions
!
!> @details
!> Note that for the naive implementation of the derivative of the kernel function, 
!! when calling dk(x1,x2), it is always assumed that x1 is the variable with respect
!! to which the kernel function should be derived. It is extremely easy to add new
!! 1-dimensional kernels simply by adding them in this module.
!-----------------------------------------------------------------------
module reproducing_kernels
    implicit none
    public
    
    !> enumerator for kernel types
    integer, parameter :: RECIPROCAL_POWER_N2_M6_KERNEL =  0, &  
                          RECIPROCAL_POWER_N3_M6_KERNEL =  1, &
                          RECIPROCAL_POWER_N2_M5_KERNEL =  2, &
                          RECIPROCAL_POWER_N3_M5_KERNEL =  3, &
                          RECIPROCAL_POWER_N2_M4_KERNEL =  4, &
                          RECIPROCAL_POWER_N3_M4_KERNEL =  5, &
                          RECIPROCAL_POWER_N2_M3_KERNEL =  6, &
                          RECIPROCAL_POWER_N3_M3_KERNEL =  7, &
                          RECIPROCAL_POWER_N2_M2_KERNEL =  8, &
                          RECIPROCAL_POWER_N3_M2_KERNEL =  9, &
                          RECIPROCAL_POWER_N2_M1_KERNEL = 10, &
                          RECIPROCAL_POWER_N3_M1_KERNEL = 11, &
                          RECIPROCAL_POWER_N2_M0_KERNEL = 12, &
                          RECIPROCAL_POWER_N3_M0_KERNEL = 13, &
                          EXPONENTIAL_DECAY_N2_KERNEL   = 14, &
                          EXPONENTIAL_DECAY_N3_KERNEL   = 15, &
                          TAYLOR_SPLINE_N2_KERNEL       = 16, &
                          TAYLOR_SPLINE_N3_KERNEL       = 17, &                          
                          LAPLACIAN_KERNEL              = 18, &
                          BERNOULLI_N2_KERNEL           = 19

    !> interface for the f2k/f3k functions that are stored in function pointer arrays
    interface func        
        pure real(kind(0d0)) function func(x,par)
            real(kind(0d0)), intent(in)               :: x   
            real(kind(0d0)), dimension(:), intent(in) :: par !optional parameters, e.g. beta for exp kernel
        end function func
    end interface func
    
    !> interface for naive implementation of kernels (instead of using the fast evaluation)
    interface kdirect
        pure real(kind(0d0)) function kdirect(x1,x2,par)
            real(kind(0d0)), intent(in)               :: x1,x2 
            real(kind(0d0)), dimension(:), intent(in) :: par
        end function kdirect
    end interface kdirect
        
    !> for wrapping function pointers, such that it is possible to have arrays of function pointers
    type f_ptr
        procedure(func), pointer, nopass :: f
    end type f_ptr
    
    !> defines 1-dimensional kernels
    type kernel_1d
        !> stores information about what type of kernel is used (-1 signals uninitialized kernels)
        integer :: kernel_type = -1
        !> how many f2k/f3k functions and p2k coefficients need to be used in the kernel decomposition
        integer :: M2 = 0
        !> pointer to naive (slow) implementation of kernel function
        procedure(kdirect), pointer, nopass :: k 
        !> pointer to naive (slow) implementation of first derivative of kernel function  
        procedure(kdirect), pointer, nopass :: dk  
        !> pointer to naive (slow) implementation of second derivative of kernel function  
        procedure(kdirect), pointer, nopass :: d2k  
        !> array of the p2k coefficients
        real(kind(0d0)), dimension(:), allocatable :: p2
        !> array of the f2k function pointers  
        type(f_ptr),     dimension(:), allocatable :: f2 
        !> array of the f3k function pointers 
        type(f_ptr),     dimension(:), allocatable :: f3 
        !> array of the first derivative of f2k function pointers 
        type(f_ptr),     dimension(:), allocatable :: df2 
        !> array of the first derivative of f3k function pointers
        type(f_ptr),     dimension(:), allocatable :: df3 
        !> array of the first derivative of f2k function pointers 
        type(f_ptr),     dimension(:), allocatable :: d2f2 
        !> array of the first derivative of f3k function pointers
        type(f_ptr),     dimension(:), allocatable :: d2f3 
        !> needed for supporting kernels with parameters
        integer :: Npar = 0
        !> array of parameters (only needed when kernel function has parameters)
        real(kind(0d0)), dimension(:), allocatable :: par 
        contains
            !> needs to be called to initialize the kernel and set it to a certain type
            procedure :: init => init_kernel
    end type kernel_1d
                                 
contains 

    !> test kernel implementation (fast vs normal) -> dval_slow and dval_fast,  d2val_slow and d2val_fast
    !! as well as fval_fast and fval_slow should be equal after calling this subroutine.
    !! If not, there must be some bug in the kernel implementation!
    subroutine debug_kernel(kernel,x1,x2,fval_slow,dval_slow,d2val_slow,fval_fast,dval_fast,d2val_fast)
        implicit none
        type(kernel_1d), intent(in)  :: kernel
        real(kind(0d0)), intent(in)  :: x1,x2     !input variables
        real(kind(0d0)), intent(out) :: fval_slow,dval_slow,d2val_slow,&
                                        fval_fast,dval_fast,d2val_fast    !function value and derivative values of 
                                                                          !slow and fast implementation
        real(kind(0d0)) :: xs,xl
        integer :: i
        
        !determine larger and smaller value (needed for fast implementation)
        if(x1 <= x2) then
            xs = x1
            xl = x2
        else
            xs = x2
            xl = x1
        end if        
        
        fval_slow  = kernel%  k(x1,x2,kernel%par) 
        dval_slow  = kernel% dk(x1,x2,kernel%par) 
        d2val_slow = kernel%d2k(x1,x2,kernel%par)
        fval_fast = 0d0
        do i = 1,kernel%M2
            fval_fast = fval_fast + kernel%p2(i) *  kernel% f2(i)%f(xs,kernel%par)*kernel% f3(i)%f(xl,kernel%par)
        end do  
        !depending on which value x1 is, we need to take the derivative of a different function
        dval_fast  = 0d0
        d2val_fast = 0d0
        if(x1 <= x2) then
            do i = 1,kernel%M2
                dval_fast  =  dval_fast + kernel%p2(i) *  kernel% df2(i)%f(xs,kernel%par)*kernel% f3(i)%f(xl,kernel%par)
                d2val_fast = d2val_fast + kernel%p2(i) *  kernel%d2f2(i)%f(xs,kernel%par)*kernel% f3(i)%f(xl,kernel%par)
            end do 
        else
            do i = 1,kernel%M2
                dval_fast  =  dval_fast + kernel%p2(i) *  kernel% f2(i)%f(xs,kernel%par)*kernel% df3(i)%f(xl,kernel%par)
                d2val_fast = d2val_fast + kernel%p2(i) *  kernel% f2(i)%f(xs,kernel%par)*kernel%d2f3(i)%f(xl,kernel%par)
            end do 
        end if

        return
    end subroutine debug_kernel 
    
    !> needs to be called to initialize the kernel and set it to a certain type
    subroutine init_kernel(kernel,kernel_type)
        use kernel_rp_2_5        
        implicit none
        class(kernel_1d) :: kernel
        integer, intent(in) :: kernel_type
                                  
        select case (kernel_type)
            case (RECIPROCAL_POWER_N2_M6_KERNEL)
                call init_rp_2_6(kernel)
            case (RECIPROCAL_POWER_N3_M6_KERNEL)
                call init_rp_3_6(kernel)
            case (RECIPROCAL_POWER_N2_M5_KERNEL)
                call init_rp_2_5(kernel)
            case (RECIPROCAL_POWER_N3_M5_KERNEL)
                call init_rp_3_5(kernel)
            case (RECIPROCAL_POWER_N2_M4_KERNEL)
                call init_rp_2_4(kernel)
            case (RECIPROCAL_POWER_N3_M4_KERNEL)
                call init_rp_3_4(kernel)
            case (RECIPROCAL_POWER_N2_M3_KERNEL)
                call init_rp_2_3(kernel)
            case (RECIPROCAL_POWER_N3_M3_KERNEL)
                call init_rp_3_3(kernel)
            case (RECIPROCAL_POWER_N2_M2_KERNEL)
                call init_rp_2_2(kernel)
            case (RECIPROCAL_POWER_N3_M2_KERNEL)
                call init_rp_3_2(kernel)
            case (RECIPROCAL_POWER_N2_M1_KERNEL)
                call init_rp_2_1(kernel)
            case (RECIPROCAL_POWER_N3_M1_KERNEL)
                call init_rp_3_1(kernel)
            case (RECIPROCAL_POWER_N2_M0_KERNEL)
                call init_rp_2_0(kernel)
            case (RECIPROCAL_POWER_N3_M0_KERNEL)
                call init_rp_3_0(kernel)
            case (EXPONENTIAL_DECAY_N2_KERNEL)
                call init_exp_2(kernel)
            case (EXPONENTIAL_DECAY_N3_KERNEL)
                call init_exp_3(kernel)
            case (TAYLOR_SPLINE_N2_KERNEL)
                call init_ts_2(kernel)
            case (TAYLOR_SPLINE_N3_KERNEL)
                call init_ts_3(kernel)
            case (LAPLACIAN_KERNEL)
                call init_laplacian(kernel)
            case (BERNOULLI_N2_KERNEL)
                call init_bernoulli_2(kernel)
            case default
                write(*,*) "ERROR in module reproducing_kernels: "//&
                           "Unknown kernel_type"
                stop
        end select
        return
    end subroutine init_kernel

    !> initializes the kernel to reciprocal decay kernel with n = 2 and m = 6
    subroutine init_rp_2_6(kernel)
        use kernel_rp_2_6        
        implicit none
        type(kernel_1d), intent(out)  :: kernel
        kernel%kernel_type = RECIPROCAL_POWER_N2_M6_KERNEL
        kernel%M2   =  M2  !M2
        kernel%Npar =  Npar!Npar                             
        kernel%k    => k   !k
        kernel%dk   => dk  !dk
        kernel%d2k  => d2k !d2k
     
        !allocate memory
        allocate(kernel%p2  (M2),&
                 kernel%f2  (M2),&
                 kernel%f3  (M2),&
                 kernel%df2 (M2),&
                 kernel%df3 (M2),&
                 kernel%d2f2(M2),&
                 kernel%d2f3(M2))
                 
        !set p2 coefficients
        kernel%p2(1) = p21
        kernel%p2(2) = p22
        
        !set f2 functions
        kernel%f2(1)%f => f21
        kernel%f2(2)%f => f22
        
        !set f3 functions
        kernel%f3(1)%f => f31
        kernel%f3(2)%f => f32        
        
        !set df2 functions
        kernel%df2(1)%f => df21
        kernel%df2(2)%f => df22
        
        !set df3 functions
        kernel%df3(1)%f => df31
        kernel%df3(2)%f => df32 
        
        !set d2f2 functions
        kernel%d2f2(1)%f => d2f21
        kernel%d2f2(2)%f => d2f22
        
        !set d2f3 functions
        kernel%d2f3(1)%f => d2f31
        kernel%d2f3(2)%f => d2f32 
        return    
    end subroutine init_rp_2_6
    
    !> initializes the kernel to reciprocal decay kernel with n = 3 and m = 6
    subroutine init_rp_3_6(kernel)
        use kernel_rp_3_6        
        implicit none
        type(kernel_1d), intent(out)  :: kernel
        kernel%kernel_type = RECIPROCAL_POWER_N3_M6_KERNEL
        kernel%M2   =  M2  !M2
        kernel%Npar =  Npar!Npar                             
        kernel%k    => k   !k
        kernel%dk   => dk  !dk
        kernel%d2k  => d2k !d2k
        
        !allocate memory
        allocate(kernel%p2  (M2),&
                 kernel%f2  (M2),&
                 kernel%f3  (M2),&
                 kernel%df2 (M2),&
                 kernel%df3 (M2),&
                 kernel%d2f2(M2),&
                 kernel%d2f3(M2))
                 
        !set p2 coefficients
        kernel%p2(1) = p21
        kernel%p2(2) = p22
        kernel%p2(3) = p23
        
        !set f2 functions
        kernel%f2(1)%f => f21
        kernel%f2(2)%f => f22
        kernel%f2(3)%f => f23 
        
        !set f3 functions
        kernel%f3(1)%f => f31
        kernel%f3(2)%f => f32        
        kernel%f3(3)%f => f33  
        
        !set df2 functions
        kernel%df2(1)%f => df21
        kernel%df2(2)%f => df22
        kernel%df2(3)%f => df23
        
        !set df3 functions
        kernel%df3(1)%f => df31
        kernel%df3(2)%f => df32  
        kernel%df3(3)%f => df33 
        
        !set d2f2 functions
        kernel%d2f2(1)%f => d2f21
        kernel%d2f2(2)%f => d2f22
        kernel%d2f2(3)%f => d2f23
        
        !set d2f3 functions
        kernel%d2f3(1)%f => d2f31
        kernel%d2f3(2)%f => d2f32  
        kernel%d2f3(3)%f => d2f33 
        return    
    end subroutine init_rp_3_6 
    
    !> initializes the kernel to reciprocal decay kernel with n = 2 and m = 5
    subroutine init_rp_2_5(kernel)
        use kernel_rp_2_5        
        implicit none
        type(kernel_1d), intent(out)  :: kernel
        kernel%kernel_type = RECIPROCAL_POWER_N2_M5_KERNEL
        kernel%M2   =  M2  !M2
        kernel%Npar =  Npar!Npar                             
        kernel%k    => k   !k
        kernel%dk   => dk  !dk
        kernel%d2k  => d2k !d2k
     
        !allocate memory
        allocate(kernel%p2  (M2),&
                 kernel%f2  (M2),&
                 kernel%f3  (M2),&
                 kernel%df2 (M2),&
                 kernel%df3 (M2),&
                 kernel%d2f2(M2),&
                 kernel%d2f3(M2))
                 
        !set p2 coefficients
        kernel%p2(1) = p21
        kernel%p2(2) = p22
        
        !set f2 functions
        kernel%f2(1)%f => f21
        kernel%f2(2)%f => f22
        
        !set f3 functions
        kernel%f3(1)%f => f31
        kernel%f3(2)%f => f32        
        
        !set df2 functions
        kernel%df2(1)%f => df21
        kernel%df2(2)%f => df22
        
        !set df3 functions
        kernel%df3(1)%f => df31
        kernel%df3(2)%f => df32 
        
        !set d2f2 functions
        kernel%d2f2(1)%f => d2f21
        kernel%d2f2(2)%f => d2f22
        
        !set d2f3 functions
        kernel%d2f3(1)%f => d2f31
        kernel%d2f3(2)%f => d2f32   
        return   
    end subroutine init_rp_2_5
    
    !> initializes the kernel to reciprocal decay kernel with n = 3 and m = 5
    subroutine init_rp_3_5(kernel)
        use kernel_rp_3_5        
        implicit none
        type(kernel_1d), intent(out)  :: kernel
        kernel%kernel_type = RECIPROCAL_POWER_N3_M5_KERNEL
        kernel%M2   =  M2  !M2
        kernel%Npar =  Npar!Npar                             
        kernel%k    => k   !k
        kernel%dk   => dk  !dk
        kernel%d2k  => d2k !d2k
        
        !allocate memory
        allocate(kernel%p2  (M2),&
                 kernel%f2  (M2),&
                 kernel%f3  (M2),&
                 kernel%df2 (M2),&
                 kernel%df3 (M2),&
                 kernel%d2f2(M2),&
                 kernel%d2f3(M2))
                 
        !set p2 coefficients
        kernel%p2(1) = p21
        kernel%p2(2) = p22
        kernel%p2(3) = p23
        
        !set f2 functions
        kernel%f2(1)%f => f21
        kernel%f2(2)%f => f22
        kernel%f2(3)%f => f23 
        
        !set f3 functions
        kernel%f3(1)%f => f31
        kernel%f3(2)%f => f32        
        kernel%f3(3)%f => f33  
        
        !set df2 functions
        kernel%df2(1)%f => df21
        kernel%df2(2)%f => df22
        kernel%df2(3)%f => df23
        
        !set df3 functions
        kernel%df3(1)%f => df31
        kernel%df3(2)%f => df32  
        kernel%df3(3)%f => df33
               
        !set d2f2 functions
        kernel%d2f2(1)%f => d2f21
        kernel%d2f2(2)%f => d2f22
        kernel%d2f2(3)%f => d2f23
        
        !set d2f3 functions
        kernel%d2f3(1)%f => d2f31
        kernel%d2f3(2)%f => d2f32  
        kernel%d2f3(3)%f => d2f33 
        return     
    end subroutine init_rp_3_5
    
    !> initializes the kernel to reciprocal decay kernel with n = 2 and m = 4
    subroutine init_rp_2_4(kernel)
        use kernel_rp_2_4        
        implicit none
        type(kernel_1d), intent(out)  :: kernel
        kernel%kernel_type = RECIPROCAL_POWER_N2_M4_KERNEL
        kernel%M2   =  M2  !M2
        kernel%Npar =  Npar!Npar                             
        kernel%k    => k   !k
        kernel%dk   => dk  !dk
        kernel%d2k  => d2k !d2k
        
        !allocate memory
        allocate(kernel%p2  (M2),&
                 kernel%f2  (M2),&
                 kernel%f3  (M2),&
                 kernel%df2 (M2),&
                 kernel%df3 (M2),&
                 kernel%d2f2(M2),&
                 kernel%d2f3(M2))
                 
        !set p2 coefficients
        kernel%p2(1) = p21
        kernel%p2(2) = p22
        
        !set f2 functions
        kernel%f2(1)%f => f21
        kernel%f2(2)%f => f22
        
        !set f3 functions
        kernel%f3(1)%f => f31
        kernel%f3(2)%f => f32        
        
        !set df2 functions
        kernel%df2(1)%f => df21
        kernel%df2(2)%f => df22
        
        !set df3 functions
        kernel%df3(1)%f => df31
        kernel%df3(2)%f => df32 
        
        !set d2f2 functions
        kernel%d2f2(1)%f => d2f21
        kernel%d2f2(2)%f => d2f22
        
        !set d2f3 functions
        kernel%d2f3(1)%f => d2f31
        kernel%d2f3(2)%f => d2f32 
        return     
    end subroutine init_rp_2_4
    
    !> initializes the kernel to reciprocal decay kernel with n = 3 and m = 4
    subroutine init_rp_3_4(kernel)
        use kernel_rp_3_4       
        implicit none
        type(kernel_1d), intent(out)  :: kernel
        kernel%kernel_type = RECIPROCAL_POWER_N3_M4_KERNEL
        kernel%M2   =  M2  !M2
        kernel%Npar =  Npar!Npar                             
        kernel%k    => k   !k
        kernel%dk   => dk  !dk
        kernel%d2k  => d2k !d2k
        
        !allocate memory
        allocate(kernel%p2  (M2),&
                 kernel%f2  (M2),&
                 kernel%f3  (M2),&
                 kernel%df2 (M2),&
                 kernel%df3 (M2),&
                 kernel%d2f2(M2),&
                 kernel%d2f3(M2))
                 
        !set p2 coefficients
        kernel%p2(1) = p21
        kernel%p2(2) = p22
        kernel%p2(3) = p23
        
        !set f2 functions
        kernel%f2(1)%f => f21
        kernel%f2(2)%f => f22
        kernel%f2(3)%f => f23 
        
        !set f3 functions
        kernel%f3(1)%f => f31
        kernel%f3(2)%f => f32        
        kernel%f3(3)%f => f33  
        
        !set df2 functions
        kernel%df2(1)%f => df21
        kernel%df2(2)%f => df22
        kernel%df2(3)%f => df23
        
        !set df3 functions
        kernel%df3(1)%f => df31
        kernel%df3(2)%f => df32  
        kernel%df3(3)%f => df33
        
        !set d2f2 functions
        kernel%d2f2(1)%f => d2f21
        kernel%d2f2(2)%f => d2f22
        kernel%d2f2(3)%f => d2f23
        
        !set d2f3 functions
        kernel%d2f3(1)%f => d2f31
        kernel%d2f3(2)%f => d2f32  
        kernel%d2f3(3)%f => d2f33 
        return     
    end subroutine init_rp_3_4
    
    !> initializes the kernel to reciprocal decay kernel with n = 2 and m = 3
    subroutine init_rp_2_3(kernel)
        use kernel_rp_2_3        
        implicit none
        type(kernel_1d), intent(out)  :: kernel
        kernel%kernel_type = RECIPROCAL_POWER_N2_M3_KERNEL
        kernel%M2   =  M2  !M2
        kernel%Npar =  Npar!Npar                             
        kernel%k    => k   !k
        kernel%dk   => dk  !dk
        kernel%d2k  => d2k !d2k
        
        !allocate memory
        allocate(kernel%p2  (M2),&
                 kernel%f2  (M2),&
                 kernel%f3  (M2),&
                 kernel%df2 (M2),&
                 kernel%df3 (M2),&
                 kernel%d2f2(M2),&
                 kernel%d2f3(M2))
                 
        !set p2 coefficients
        kernel%p2(1) = p21
        kernel%p2(2) = p22
        
        !set f2 functions
        kernel%f2(1)%f => f21
        kernel%f2(2)%f => f22
        
        !set f3 functions
        kernel%f3(1)%f => f31
        kernel%f3(2)%f => f32        
        
        !set df2 functions
        kernel%df2(1)%f => df21
        kernel%df2(2)%f => df22
        
        !set df3 functions
        kernel%df3(1)%f => df31
        kernel%df3(2)%f => df32
        
        !set d2f2 functions
        kernel%d2f2(1)%f => d2f21
        kernel%d2f2(2)%f => d2f22
        
        !set d2f3 functions
        kernel%d2f3(1)%f => d2f31
        kernel%d2f3(2)%f => d2f32 
        return      
    end subroutine init_rp_2_3
    
    !> initializes the kernel to reciprocal decay kernel with n = 3 and m = 3
    subroutine init_rp_3_3(kernel)
        use kernel_rp_3_3       
        implicit none
        type(kernel_1d), intent(out)  :: kernel
        kernel%kernel_type = RECIPROCAL_POWER_N3_M3_KERNEL
        kernel%M2   =  M2  !M2
        kernel%Npar =  Npar!Npar                             
        kernel%k    => k   !k
        kernel%dk   => dk  !dk
        kernel%d2k  => d2k !d2k
        
        !allocate memory
        allocate(kernel%p2  (M2),&
                 kernel%f2  (M2),&
                 kernel%f3  (M2),&
                 kernel%df2 (M2),&
                 kernel%df3 (M2),&
                 kernel%d2f2(M2),&
                 kernel%d2f3(M2))
                 
        !set p2 coefficients
        kernel%p2(1) = p21
        kernel%p2(2) = p22
        kernel%p2(3) = p23
        
        !set f2 functions
        kernel%f2(1)%f => f21
        kernel%f2(2)%f => f22
        kernel%f2(3)%f => f23 
        
        !set f3 functions
        kernel%f3(1)%f => f31
        kernel%f3(2)%f => f32        
        kernel%f3(3)%f => f33  
        
        !set df2 functions
        kernel%df2(1)%f => df21
        kernel%df2(2)%f => df22
        kernel%df2(3)%f => df23
        
        !set df3 functions
        kernel%df3(1)%f => df31
        kernel%df3(2)%f => df32  
        kernel%df3(3)%f => df33
               
        !set d2f2 functions
        kernel%d2f2(1)%f => d2f21
        kernel%d2f2(2)%f => d2f22
        kernel%d2f2(3)%f => d2f23
        
        !set d2f3 functions
        kernel%d2f3(1)%f => d2f31
        kernel%d2f3(2)%f => d2f32  
        kernel%d2f3(3)%f => d2f33 
        return     
    end subroutine init_rp_3_3
    
    !> initializes the kernel to reciprocal decay kernel with n = 2 and m = 2
    subroutine init_rp_2_2(kernel)
        use kernel_rp_2_2      
        implicit none
        type(kernel_1d), intent(out)  :: kernel
        kernel%kernel_type = RECIPROCAL_POWER_N2_M2_KERNEL
        kernel%M2   =  M2  !M2
        kernel%Npar =  Npar!Npar                             
        kernel%k    => k   !k
        kernel%dk   => dk  !dk
        kernel%d2k  => d2k !d2k
        
        !allocate memory
        allocate(kernel%p2  (M2),&
                 kernel%f2  (M2),&
                 kernel%f3  (M2),&
                 kernel%df2 (M2),&
                 kernel%df3 (M2),&
                 kernel%d2f2(M2),&
                 kernel%d2f3(M2))
                 
        !set p2 coefficients
        kernel%p2(1) = p21
        kernel%p2(2) = p22
        
        !set f2 functions
        kernel%f2(1)%f => f21
        kernel%f2(2)%f => f22
        
        !set f3 functions
        kernel%f3(1)%f => f31
        kernel%f3(2)%f => f32        
        
        !set df2 functions
        kernel%df2(1)%f => df21
        kernel%df2(2)%f => df22
        
        !set df3 functions
        kernel%df3(1)%f => df31
        kernel%df3(2)%f => df32 
        
        !set d2f2 functions
        kernel%d2f2(1)%f => d2f21
        kernel%d2f2(2)%f => d2f22
        
        !set d2f3 functions
        kernel%d2f3(1)%f => d2f31
        kernel%d2f3(2)%f => d2f32 
        return     
    end subroutine init_rp_2_2
    
    !> initializes the kernel to reciprocal decay kernel with n = 3 and m = 2
    subroutine init_rp_3_2(kernel)
        use kernel_rp_3_2       
        implicit none
        type(kernel_1d), intent(out)  :: kernel
        kernel%kernel_type = RECIPROCAL_POWER_N3_M2_KERNEL
        kernel%M2   =  M2  !M2
        kernel%Npar =  Npar!Npar                             
        kernel%k    => k   !k
        kernel%dk   => dk  !dk
        kernel%d2k  => d2k !d2k
        
        !allocate memory
        allocate(kernel%p2  (M2),&
                 kernel%f2  (M2),&
                 kernel%f3  (M2),&
                 kernel%df2 (M2),&
                 kernel%df3 (M2),&
                 kernel%d2f2(M2),&
                 kernel%d2f3(M2))
                 
        !set p2 coefficients
        kernel%p2(1) = p21
        kernel%p2(2) = p22
        kernel%p2(3) = p23
        
        !set f2 functions
        kernel%f2(1)%f => f21
        kernel%f2(2)%f => f22
        kernel%f2(3)%f => f23 
        
        !set f3 functions
        kernel%f3(1)%f => f31
        kernel%f3(2)%f => f32        
        kernel%f3(3)%f => f33  
        
        !set df2 functions
        kernel%df2(1)%f => df21
        kernel%df2(2)%f => df22
        kernel%df2(3)%f => df23
        
        !set df3 functions
        kernel%df3(1)%f => df31
        kernel%df3(2)%f => df32  
        kernel%df3(3)%f => df33
        
        !set d2f2 functions
        kernel%d2f2(1)%f => d2f21
        kernel%d2f2(2)%f => d2f22
        kernel%d2f2(3)%f => d2f23
        
        !set d2f3 functions
        kernel%d2f3(1)%f => d2f31
        kernel%d2f3(2)%f => d2f32  
        kernel%d2f3(3)%f => d2f33     
        return 
    end subroutine init_rp_3_2
    
    !> initializes the kernel to reciprocal decay kernel with n = 2 and m = 1
    subroutine init_rp_2_1(kernel)
        use kernel_rp_2_1       
        implicit none
        type(kernel_1d), intent(out)  :: kernel
        kernel%kernel_type = RECIPROCAL_POWER_N2_M1_KERNEL
        kernel%M2   =  M2  !M2
        kernel%Npar =  Npar!Npar                             
        kernel%k    => k   !k
        kernel%dk   => dk  !dk
        kernel%d2k  => d2k !d2k
        
        !allocate memory
        allocate(kernel%p2  (M2),&
                 kernel%f2  (M2),&
                 kernel%f3  (M2),&
                 kernel%df2 (M2),&
                 kernel%df3 (M2),&
                 kernel%d2f2(M2),&
                 kernel%d2f3(M2))
                 
        !set p2 coefficients
        kernel%p2(1) = p21
        kernel%p2(2) = p22
        
        !set f2 functions
        kernel%f2(1)%f => f21
        kernel%f2(2)%f => f22
        
        !set f3 functions
        kernel%f3(1)%f => f31
        kernel%f3(2)%f => f32        
        
        !set df2 functions
        kernel%df2(1)%f => df21
        kernel%df2(2)%f => df22
        
        !set df3 functions
        kernel%df3(1)%f => df31
        kernel%df3(2)%f => df32   
        
        !set d2f2 functions
        kernel%d2f2(1)%f => d2f21
        kernel%d2f2(2)%f => d2f22
        
        !set d2f3 functions
        kernel%d2f3(1)%f => d2f31
        kernel%d2f3(2)%f => d2f32  
        return  
    end subroutine init_rp_2_1
    
    !> initializes the kernel to reciprocal decay kernel with n = 3 and m = 1
    subroutine init_rp_3_1(kernel)
        use kernel_rp_3_1       
        implicit none
        type(kernel_1d), intent(out)  :: kernel
        kernel%kernel_type = RECIPROCAL_POWER_N3_M1_KERNEL
        kernel%M2   =  M2  !M2
        kernel%Npar =  Npar!Npar                             
        kernel%k    => k   !k
        kernel%dk   => dk  !dk
        kernel%d2k  => d2k !d2k
        
        !allocate memory
        allocate(kernel%p2  (M2),&
                 kernel%f2  (M2),&
                 kernel%f3  (M2),&
                 kernel%df2 (M2),&
                 kernel%df3 (M2),&
                 kernel%d2f2(M2),&
                 kernel%d2f3(M2))
                 
        !set p2 coefficients
        kernel%p2(1) = p21
        kernel%p2(2) = p22
        kernel%p2(3) = p23
        
        !set f2 functions
        kernel%f2(1)%f => f21
        kernel%f2(2)%f => f22
        kernel%f2(3)%f => f23 
        
        !set f3 functions
        kernel%f3(1)%f => f31
        kernel%f3(2)%f => f32        
        kernel%f3(3)%f => f33  
        
        !set df2 functions
        kernel%df2(1)%f => df21
        kernel%df2(2)%f => df22
        kernel%df2(3)%f => df23
        
        !set df3 functions
        kernel%df3(1)%f => df31
        kernel%df3(2)%f => df32  
        kernel%df3(3)%f => df33 
        
        !set d2f2 functions
        kernel%d2f2(1)%f => d2f21
        kernel%d2f2(2)%f => d2f22
        kernel%d2f2(3)%f => d2f23
        
        !set d2f3 functions
        kernel%d2f3(1)%f => d2f31
        kernel%d2f3(2)%f => d2f32  
        kernel%d2f3(3)%f => d2f33    
        return 
    end subroutine init_rp_3_1
    
    !> initializes the kernel to reciprocal decay kernel with n = 2 and m = 0
    subroutine init_rp_2_0(kernel)
        use kernel_rp_2_0        
        implicit none
        type(kernel_1d), intent(out)  :: kernel
        kernel%kernel_type = RECIPROCAL_POWER_N2_M0_KERNEL
        kernel%M2   =  M2  !M2
        kernel%Npar =  Npar!Npar                             
        kernel%k    => k   !k
        kernel%dk   => dk  !dk
        kernel%d2k  => d2k !d2k
        
        !allocate memory
        allocate(kernel%p2  (M2),&
                 kernel%f2  (M2),&
                 kernel%f3  (M2),&
                 kernel%df2 (M2),&
                 kernel%df3 (M2),&
                 kernel%d2f2(M2),&
                 kernel%d2f3(M2))
                 
        !set p2 coefficients
        kernel%p2(1) = p21
        kernel%p2(2) = p22
        
        !set f2 functions
        kernel%f2(1)%f => f21
        kernel%f2(2)%f => f22
        
        !set f3 functions
        kernel%f3(1)%f => f31
        kernel%f3(2)%f => f32        
        
        !set df2 functions
        kernel%df2(1)%f => df21
        kernel%df2(2)%f => df22
        
        !set df3 functions
        kernel%df3(1)%f => df31
        kernel%df3(2)%f => df32
        
        !set d2f2 functions
        kernel%d2f2(1)%f => d2f21
        kernel%d2f2(2)%f => d2f22
        
        !set d2f3 functions
        kernel%d2f3(1)%f => d2f31
        kernel%d2f3(2)%f => d2f32 
        return      
    end subroutine init_rp_2_0
    
    !> initializes the kernel to reciprocal decay kernel with n = 3 and m = 0
    subroutine init_rp_3_0(kernel)
        use kernel_rp_3_0       
        implicit none
        type(kernel_1d), intent(out)  :: kernel
        kernel%kernel_type = RECIPROCAL_POWER_N3_M0_KERNEL
        kernel%M2   =  M2  !M2
        kernel%Npar =  Npar!Npar                             
        kernel%k    => k   !k
        kernel%dk   => dk  !dk
        kernel%d2k  => d2k !d2k
        
        !allocate memory
        allocate(kernel%p2  (M2),&
                 kernel%f2  (M2),&
                 kernel%f3  (M2),&
                 kernel%df2 (M2),&
                 kernel%df3 (M2),&
                 kernel%d2f2(M2),&
                 kernel%d2f3(M2))
                 
        !set p2 coefficients
        kernel%p2(1) = p21
        kernel%p2(2) = p22
        kernel%p2(3) = p23
        
        !set f2 functions
        kernel%f2(1)%f => f21
        kernel%f2(2)%f => f22
        kernel%f2(3)%f => f23 
        
        !set f3 functions
        kernel%f3(1)%f => f31
        kernel%f3(2)%f => f32        
        kernel%f3(3)%f => f33  
        
        !set df2 functions
        kernel%df2(1)%f => df21
        kernel%df2(2)%f => df22
        kernel%df2(3)%f => df23
        
        !set df3 functions
        kernel%df3(1)%f => df31
        kernel%df3(2)%f => df32  
        kernel%df3(3)%f => df33
        
        !set d2f2 functions
        kernel%d2f2(1)%f => d2f21
        kernel%d2f2(2)%f => d2f22
        kernel%d2f2(3)%f => d2f23
        
        !set d2f3 functions
        kernel%d2f3(1)%f => d2f31
        kernel%d2f3(2)%f => d2f32  
        kernel%d2f3(3)%f => d2f33     
        return 
    end subroutine init_rp_3_0
    
    !> initializes the kernel to exponential decay kernel with n = 2
    subroutine init_exp_2(kernel)
        use kernel_exp_2    
        implicit none
        type(kernel_1d), intent(out)  :: kernel
        kernel%kernel_type = EXPONENTIAL_DECAY_N2_KERNEL
        kernel%M2   = M2   !M2
        kernel%Npar = Npar !Npar
        if(allocated(kernel%par)) deallocate(kernel%par)
        allocate(kernel%par(Npar))
        kernel%par = 1d0                                                      
        kernel%k    => k   !k
        kernel%dk   => dk  !dk
        kernel%d2k  => d2k !d2k
        
        !allocate memory
        allocate(kernel%p2  (M2),&
                 kernel%f2  (M2),&
                 kernel%f3  (M2),&
                 kernel%df2 (M2),&
                 kernel%df3 (M2),&
                 kernel%d2f2(M2),&
                 kernel%d2f3(M2))
                 
        !set p2 coefficients
        kernel%p2(1) = p21
        kernel%p2(2) = p22
        
        !set f2 functions
        kernel%f2(1)%f => f21
        kernel%f2(2)%f => f22
        
        !set f3 functions
        kernel%f3(1)%f => f31
        kernel%f3(2)%f => f32        
        
        !set df2 functions
        kernel%df2(1)%f => df21
        kernel%df2(2)%f => df22
        
        !set df3 functions
        kernel%df3(1)%f => df31
        kernel%df3(2)%f => df32
        
        !set d2f2 functions
        kernel%d2f2(1)%f => d2f21
        kernel%d2f2(2)%f => d2f22
        
        !set d2f3 functions
        kernel%d2f3(1)%f => d2f31
        kernel%d2f3(2)%f => d2f32
        return      
    end subroutine init_exp_2
    
    !> initializes the kernel to exponential decay kernel with n = 3
    subroutine init_exp_3(kernel)
        use kernel_exp_3   
        implicit none
        type(kernel_1d), intent(out)  :: kernel
        kernel%kernel_type = EXPONENTIAL_DECAY_N3_KERNEL
        kernel%M2   = M2   !M2
        kernel%Npar = Npar !Npar
        if(allocated(kernel%par)) deallocate(kernel%par)
        allocate(kernel%par(Npar))
        kernel%par = 1d0                           
        kernel%k    => k   !k
        kernel%dk   => dk  !dk
        kernel%d2k  => d2k !d2k
        
        !allocate memory
        allocate(kernel%p2  (M2),&
                 kernel%f2  (M2),&
                 kernel%f3  (M2),&
                 kernel%df2 (M2),&
                 kernel%df3 (M2),&
                 kernel%d2f2(M2),&
                 kernel%d2f3(M2))
                 
        !set p2 coefficients
        kernel%p2(1) = p21
        kernel%p2(2) = p22
        kernel%p2(3) = p23
        kernel%p2(4) = p24
        kernel%p2(5) = p25
        
        !set f2 functions
        kernel%f2(1)%f => f21
        kernel%f2(2)%f => f22
        kernel%f2(3)%f => f23
        kernel%f2(4)%f => f24
        kernel%f2(5)%f => f25
        
        !set f3 functions
        kernel%f3(1)%f => f31
        kernel%f3(2)%f => f32    
        kernel%f3(3)%f => f33
        kernel%f3(4)%f => f34       
        kernel%f3(5)%f => f35
            
        !set df2 functions
        kernel%df2(1)%f => df21
        kernel%df2(2)%f => df22
        kernel%df2(3)%f => df23
        kernel%df2(4)%f => df24
        kernel%df2(5)%f => df25
        
        !set df3 functions
        kernel%df3(1)%f => df31
        kernel%df3(2)%f => df32    
        kernel%df3(3)%f => df33
        kernel%df3(4)%f => df34    
        kernel%df3(5)%f => df35
            
        !set d2f2 functions
        kernel%d2f2(1)%f => d2f21
        kernel%d2f2(2)%f => d2f22
        kernel%d2f2(3)%f => d2f23
        kernel%d2f2(4)%f => d2f24
        kernel%d2f2(5)%f => d2f25
        
        !set d2f3 functions
        kernel%d2f3(1)%f => d2f31
        kernel%d2f3(2)%f => d2f32    
        kernel%d2f3(3)%f => d2f33
        kernel%d2f3(4)%f => d2f34    
        kernel%d2f3(5)%f => d2f35
        return 
    end subroutine init_exp_3
    
    !> initializes the kernel to Taylor spline kernel with n = 2
    subroutine init_ts_2(kernel)
        use kernel_ts_2    
        implicit none
        type(kernel_1d), intent(out)  :: kernel
        kernel%kernel_type = TAYLOR_SPLINE_N2_KERNEL
        kernel%M2   =  M2  !M2
        kernel%Npar =  Npar!Npar                             
        kernel%k    => k   !k
        kernel%dk   => dk  !dk
        kernel%d2k  => d2k !d2k
        
        !allocate memory
        allocate(kernel%p2  (M2),&
                 kernel%f2  (M2),&
                 kernel%f3  (M2),&
                 kernel%df2 (M2),&
                 kernel%df3 (M2),&
                 kernel%d2f2(M2),&
                 kernel%d2f3(M2))
                 
        !set p2 coefficients
        kernel%p2(1) = p21
        kernel%p2(2) = p22
        kernel%p2(3) = p23
        
        !set f2 functions
        kernel%f2(1)%f => f21
        kernel%f2(2)%f => f22
        kernel%f2(3)%f => f23
        
        !set f3 functions
        kernel%f3(1)%f => f31
        kernel%f3(2)%f => f32        
        kernel%f3(3)%f => f33    
        
        !set df2 functions
        kernel%df2(1)%f => df21
        kernel%df2(2)%f => df22
        kernel%df2(3)%f => df23 
        
        !set df3 functions
        kernel%df3(1)%f => df31
        kernel%df3(2)%f => df32
        kernel%df3(3)%f => df33
        
        !set d2f2 functions
        kernel%d2f2(1)%f => d2f21
        kernel%d2f2(2)%f => d2f22
        kernel%d2f2(3)%f => d2f23 
        
        !set d2f3 functions
        kernel%d2f3(1)%f => d2f31
        kernel%d2f3(2)%f => d2f32
        kernel%d2f3(3)%f => d2f33     
        return 
    end subroutine init_ts_2
    
    !> initializes the kernel to Taylor spline kernel with n = 3
    subroutine init_ts_3(kernel)
        use kernel_ts_3   
        implicit none
        type(kernel_1d), intent(out)  :: kernel
        kernel%kernel_type = TAYLOR_SPLINE_N3_KERNEL
        kernel%M2   =  M2  !M2
        kernel%Npar =  Npar!Npar                             
        kernel%k    => k   !k
        kernel%dk   => dk  !dk
        kernel%d2k  => d2k !d2k
        
        !allocate memory
        allocate(kernel%p2  (M2),&
                 kernel%f2  (M2),&
                 kernel%f3  (M2),&
                 kernel%df2 (M2),&
                 kernel%df3 (M2),&
                 kernel%d2f2(M2),&
                 kernel%d2f3(M2))
                 
        !set p2 coefficients
        kernel%p2(1) = p21
        kernel%p2(2) = p22
        kernel%p2(3) = p23
        kernel%p2(4) = p24
        kernel%p2(5) = p25
        
        !set f2 functions
        kernel%f2(1)%f => f21
        kernel%f2(2)%f => f22
        kernel%f2(3)%f => f23
        kernel%f2(4)%f => f24
        kernel%f2(5)%f => f25
        
        !set f3 functions
        kernel%f3(1)%f => f31
        kernel%f3(2)%f => f32        
        kernel%f3(3)%f => f33 
        kernel%f3(4)%f => f34        
        kernel%f3(5)%f => f35    
        
        !set df2 functions
        kernel%df2(1)%f => df21
        kernel%df2(2)%f => df22
        kernel%df2(3)%f => df23
        kernel%df2(4)%f => df24
        kernel%df2(5)%f => df25  
        
        !set df3 functions
        kernel%df3(1)%f => df31
        kernel%df3(2)%f => df32
        kernel%df3(3)%f => df33   
        kernel%df3(4)%f => df34
        kernel%df3(5)%f => df35 
        
        !set d2f2 functions
        kernel%d2f2(1)%f => d2f21
        kernel%d2f2(2)%f => d2f22
        kernel%d2f2(3)%f => d2f23
        kernel%d2f2(4)%f => d2f24
        kernel%d2f2(5)%f => d2f25  
        
        !set d2f3 functions
        kernel%d2f3(1)%f => d2f31
        kernel%d2f3(2)%f => d2f32
        kernel%d2f3(3)%f => d2f33   
        kernel%d2f3(4)%f => d2f34
        kernel%d2f3(5)%f => d2f35
        return 
    end subroutine init_ts_3   
    
    !> initializes the kernel to Laplacian kernel
    subroutine init_laplacian(kernel)
        use kernel_laplacian  
        implicit none
        type(kernel_1d), intent(out)  :: kernel
        kernel%kernel_type = LAPLACIAN_KERNEL
        kernel%M2   = M2   !M2
        kernel%Npar = Npar !Npar
        if(allocated(kernel%par)) deallocate(kernel%par)
        allocate(kernel%par(Npar))
        kernel%par = 1d0                                                      
        kernel%k    => k   !k
        kernel%dk   => dk  !dk
        kernel%d2k  => d2k !d2k
        
        !allocate memory
        allocate(kernel%p2  (M2),&
                 kernel%f2  (M2),&
                 kernel%f3  (M2),&
                 kernel%df2 (M2),&
                 kernel%df3 (M2),&
                 kernel%d2f2(M2),&
                 kernel%d2f3(M2))
                 
        !set p2 coefficients
        kernel%p2(1) = p21
        
        !set f2 functions
        kernel%f2(1)%f => f21
        
        !set f3 functions
        kernel%f3(1)%f => f31        
        
        !set df2 functions
        kernel%df2(1)%f => df21
        
        !set df3 functions
        kernel%df3(1)%f => df31 
        
        !set d2f2 functions
        kernel%d2f2(1)%f => d2f21
        
        !set d2f3 functions
        kernel%d2f3(1)%f => d2f31
        return  
    end subroutine init_laplacian
    
    !> initializes the kernel to Bernoulli kernel with n = 2
    subroutine init_bernoulli_2(kernel)
        use kernel_bernoulli_2
        implicit none
        type(kernel_1d), intent(out)  :: kernel
        kernel%kernel_type = BERNOULLI_N2_KERNEL
        kernel%M2   =  M2  !M2
        kernel%Npar =  Npar!Npar                             
        kernel%k    => k   !k
        kernel%dk   => dk  !dk
        kernel%d2k  => d2k !d2k
        
        !allocate memory
        allocate(kernel%p2  (M2),&
                 kernel%f2  (M2),&
                 kernel%f3  (M2),&
                 kernel%df2 (M2),&
                 kernel%df3 (M2),&
                 kernel%d2f2(M2),&
                 kernel%d2f3(M2))
                 
        !set p2 coefficients
        kernel%p2(1)  = p21
        kernel%p2(2)  = p22
        kernel%p2(3)  = p23
        kernel%p2(4)  = p24 
        kernel%p2(5)  = p25 
        kernel%p2(6)  = p26 
        kernel%p2(7)  = p27 
        kernel%p2(8)  = p28 
        kernel%p2(9)  = p29 
        kernel%p2(10) = p210 
        kernel%p2(11) = p211 
        kernel%p2(12) = p212 
        
        !set f2 functions
        kernel%f2(1)%f  => f21 
        kernel%f2(2)%f  => f22
        kernel%f2(3)%f  => f23 
        kernel%f2(4)%f  => f24 
        kernel%f2(5)%f  => f25 
        kernel%f2(6)%f  => f26 
        kernel%f2(7)%f  => f27 
        kernel%f2(8)%f  => f28 
        kernel%f2(9)%f  => f29 
        kernel%f2(10)%f => f210 
        kernel%f2(11)%f => f211 
        kernel%f2(12)%f => f212 
        
        !set f3 functions
        kernel%f3(1)%f  => f31  
        kernel%f3(2)%f  => f32    
        kernel%f3(3)%f  => f33    
        kernel%f3(4)%f  => f34    
        kernel%f3(5)%f  => f35    
        kernel%f3(6)%f  => f36    
        kernel%f3(7)%f  => f37    
        kernel%f3(8)%f  => f38    
        kernel%f3(9)%f  => f39    
        kernel%f3(10)%f => f310    
        kernel%f3(11)%f => f311    
        kernel%f3(12)%f => f312             
        
        !set df2 functions
        kernel%df2(1)%f  => df21
        kernel%df2(2)%f  => df22 
        kernel%df2(3)%f  => df23 
        kernel%df2(4)%f  => df24 
        kernel%df2(5)%f  => df25 
        kernel%df2(6)%f  => df26 
        kernel%df2(7)%f  => df27 
        kernel%df2(8)%f  => df28 
        kernel%df2(9)%f  => df29 
        kernel%df2(10)%f => df210 
        kernel%df2(11)%f => df211 
        kernel%df2(12)%f => df212 
        
        !set df3 functions
        kernel%df3(1)%f  => df31  
        kernel%df3(2)%f  => df32 
        kernel%df3(3)%f  => df33 
        kernel%df3(4)%f  => df34 
        kernel%df3(5)%f  => df35 
        kernel%df3(6)%f  => df36 
        kernel%df3(7)%f  => df37 
        kernel%df3(8)%f  => df38 
        kernel%df3(9)%f  => df39 
        kernel%df3(10)%f => df310 
        kernel%df3(11)%f => df311 
        kernel%df3(12)%f => df312 
        
        !set d2f2 functions
        kernel%d2f2(1)%f  => d2f21
        kernel%d2f2(2)%f  => d2f22 
        kernel%d2f2(3)%f  => d2f23 
        kernel%d2f2(4)%f  => d2f24 
        kernel%d2f2(5)%f  => d2f25 
        kernel%d2f2(6)%f  => d2f26 
        kernel%d2f2(7)%f  => d2f27 
        kernel%d2f2(8)%f  => d2f28 
        kernel%d2f2(9)%f  => d2f29 
        kernel%d2f2(10)%f => d2f210 
        kernel%d2f2(11)%f => d2f211 
        kernel%d2f2(12)%f => d2f212 
        
        !set d2f3 functions
        kernel%d2f3(1)%f  => d2f31  
        kernel%d2f3(2)%f  => d2f32 
        kernel%d2f3(3)%f  => d2f33 
        kernel%d2f3(4)%f  => d2f34 
        kernel%d2f3(5)%f  => d2f35 
        kernel%d2f3(6)%f  => d2f36 
        kernel%d2f3(7)%f  => d2f37 
        kernel%d2f3(8)%f  => d2f38 
        kernel%d2f3(9)%f  => d2f39 
        kernel%d2f3(10)%f => d2f310 
        kernel%d2f3(11)%f => d2f311 
        kernel%d2f3(12)%f => d2f312 
        return 
    end subroutine init_bernoulli_2
        
end module reproducing_kernels


!-----------------------------------------------------------------------
!> @author
!> Oliver T. Unke, University of Basel
!
!> @brief
!> main RKHS module (this is the only module that needs to be directly used)
!
!> @details
!> This module contains all necessary type definitions and functions in order
!! to setup an arbitrary kernel for interpolating multi-dimensional functions.
!! The method was first described in:
!! T.-S. Ho and H. Rabitz. "A general method for constructing multidimensional molecular potential energy surfaces from abinitio calculations." The Journal of chemical physics 104.7 (1996): 2584-2597.
!! The fast algorithm was first described in:
!! T. Hollebeek, T.-S. Ho, and H. Rabitz. "A fast algorithm for evaluating multidimensional potential energy surfaces." The Journal of chemical physics 106.17 (1997): 7223-7227.
!! And a useful review article can be found in:
!! T. Hollebeek, T.-S. Ho, and H. Rabitz. "Constructing multidimensional molecular potential energy surfaces from ab initio data." Annual review of physical chemistry 50.1 (1999): 537-570.
!-----------------------------------------------------------------------
module RKHS
use file_io
use reproducing_kernels
use LinearAlgebra
implicit none
private
public :: kernel, &
          debug_kernel, &
          RECIPROCAL_POWER_N2_M6_KERNEL , &
          RECIPROCAL_POWER_N3_M6_KERNEL , &
          RECIPROCAL_POWER_N2_M5_KERNEL , &
          RECIPROCAL_POWER_N3_M5_KERNEL , &
          RECIPROCAL_POWER_N2_M4_KERNEL , &
          RECIPROCAL_POWER_N3_M4_KERNEL , &
          RECIPROCAL_POWER_N2_M3_KERNEL , &
          RECIPROCAL_POWER_N3_M3_KERNEL , &
          RECIPROCAL_POWER_N2_M2_KERNEL , &
          RECIPROCAL_POWER_N3_M2_KERNEL , &
          RECIPROCAL_POWER_N2_M1_KERNEL , &
          RECIPROCAL_POWER_N3_M1_KERNEL , &
          RECIPROCAL_POWER_N2_M0_KERNEL , &
          RECIPROCAL_POWER_N3_M0_KERNEL , &
          EXPONENTIAL_DECAY_N2_KERNEL   , &
          EXPONENTIAL_DECAY_N3_KERNEL   , &
          TAYLOR_SPLINE_N2_KERNEL       , &
          TAYLOR_SPLINE_N3_KERNEL       , &
          LAPLACIAN_KERNEL              , &
          BERNOULLI_N2_KERNEL           

!> contained in kernel type, stores grid points for each dimension
type grid
    real(kind(0d0)), dimension(:),   allocatable :: x
end type grid

!> wrapper for matrices, used to store a kernel matrix in tensor product form 
type kernel_matrix
    real(kind(0d0)), dimension(:,:), allocatable :: M
end type kernel_matrix

!> multi-dimensional kernel type (contains 1-d kernels and lookup tables, coefficients, etc.)
type kernel
    !> number of dimensions/1-dimensional kernel functions
    integer                                      :: Ndim = 0
    !> total number of reference points 
    integer                                      :: Npoint = 0
    !> array of 1-dimensional kernel functions
    type(kernel_1d), dimension(:),   allocatable :: k1d 
    !> stores the values of the grid points of each 1-dimensional grid 
    type(grid),      dimension(:),   allocatable :: grid
    !> stores the values of the Npoint reference points 
    real(kind(0d0)), dimension(:),   allocatable :: values 
    !> stores the inverse of the kernel matrix (needed for calculating error bound)
    real(kind(0d0)), dimension(:,:), allocatable :: invQ     
    !> stores which values of the reference points are not missing
    logical,         dimension(:),   allocatable :: valueIsPresent    
    !> stores the kernel coefficients (obtained by matrix inversion)
    real(kind(0d0)), dimension(:),   allocatable :: alpha  
    !> stores the lookup table of precomputed sums
    real(kind(0d0)), dimension(:,:), allocatable :: sigma  
    !> stores intermediate results in fast evaluation 
    !(defined here instead of in function to avoid frequent memory allocation)
    real(kind(0d0)), dimension(:),   allocatable :: gamma_old, gamma_new  
    
    contains
        procedure :: save_to_file                => save_kernel_to_file
        procedure :: load_from_file              => load_kernel_from_file
        procedure :: free                        => deallocate_kernel
        procedure :: read_grid                   => read_grid_data
        procedure :: write_grid                  => write_grid_data
        procedure :: evaluate_fast               => evaluate_kernel_fast
        procedure :: calculate_error_bound       => calculate_error_bound
        procedure :: evaluate_slow               => evaluate_kernel_slow
        procedure :: calculate_coefficients_fast => calculate_coefficients_fast
        procedure :: calculate_coefficients_slow => calculate_coefficients_slow
        procedure :: calculate_sums              => calculate_sums
        procedure :: get_alpha_idx               => get_idx_from_indices
        procedure :: get_sig_idx                 => get_sig_idx_from_m_and_k
        procedure :: get_lookup_idx              => get_sig_idx_from_indices
end type kernel

contains 
!-----------------------------------------------------------------------
!> @author
!> Oliver T. Unke, University of Basel
!
!> @brief
!> Throws error messages if problems are encountered
!-----------------------------------------------------------------------
subroutine throw_error(proc,errormsg)
    implicit none
    character(len=*), intent(in) :: proc, errormsg
    
    write(*,*) "ERROR in module RKHS: "//proc//": "//errormsg
    stop
end subroutine throw_error

!-----------------------------------------------------------------------
!> @author
!> Oliver T. Unke, University of Basel
!
!> @brief
!> Saves kernel in a binary file, so that it can be reused later
!-----------------------------------------------------------------------
subroutine save_kernel_to_file(Q,datafile)
    implicit none
    character(len=*), intent(in) :: datafile
    class(kernel)                :: Q
    integer                      :: i, ios
    
    open(io_unit, file=datafile, form="unformatted", access="stream", &
                            action="write",     status="replace", iostat = ios)
    if(ios /= 0) then
        call throw_error("save_kernel_to_file",'Could not open "'//datafile//'".')
    end if
    write(io_unit) Q%Ndim                        
    write(io_unit) Q%Npoint
    if(allocated(Q%k1d)) then
        write(io_unit) .true.
        do i = 1,size(Q%k1d)
            write(io_unit) Q%k1d(i)%kernel_type
            if(Q%k1d(i)%Npar > 0) then
                write(io_unit) Q%k1d(i)%par    
            end if
        end do 
    else
        write(io_unit) .false.
    end if
    
    if(allocated(Q%grid)) then
        write(io_unit) .true.
        do i = 1,size(Q%grid)
            if(allocated(Q%grid(i)%x)) then
                write(io_unit) .true. 
                write(io_unit) size(Q%grid(i)%x)
                write(io_unit) Q%grid(i)%x
            else
                write(io_unit) .false.
            end if
        end do 
    else
        write(io_unit) .false.
    end if   
    
    if(allocated(Q%values)) then
        write(io_unit) .true.
        write(io_unit) Q%values 
    else
        write(io_unit) .false.
    end if
    
    if(allocated(Q%valueIsPresent)) then
        write(io_unit) .true.
        write(io_unit) Q%valueIsPresent
    else
        write(io_unit) .false.
    end if
    
    if(allocated(Q%alpha)) then
        write(io_unit) .true.
        write(io_unit) Q%alpha 
    else
        write(io_unit) .false.
    end if
    
    if(allocated(Q%invQ)) then
        write(io_unit) .true.
        write(io_unit) Q%invQ
    else
        write(io_unit) .false.
    end if
        
    if(allocated(Q%sigma)) then
        write(io_unit) .true.
        write(io_unit) Q%sigma
    else
        write(io_unit) .false.
    end if
                        
    close(io_unit)
    return 
end subroutine save_kernel_to_file

!-----------------------------------------------------------------------
!> @author
!> Oliver T. Unke, University of Basel
!
!> @brief
!> Loads kernel from a binary file
!-----------------------------------------------------------------------
subroutine load_kernel_from_file(Q,datafile)
    implicit none
    character(len=*), intent(in) :: datafile
    class(kernel)                :: Q
    integer                      :: i, ios, grid_size, kernel_type
    logical                      :: isWritten
    
    open(io_unit, file=datafile, form="unformatted", access="stream", &
                            action="read", status="old", iostat = ios)
    if(ios /= 0) call throw_error("load_kernel_from_file",'Could not open "'//datafile//'".')
    
    call Q%free() !make sure that Q is a completely fresh kernel
    
    read(30, iostat = ios) Q%Ndim 
    if(ios /= 0) call throw_error("load_kernel_from_file",'Could not read "'//datafile//'".') 
                          
    read(30, iostat = ios) Q%Npoint
    if(ios /= 0) call throw_error("load_kernel_from_file",'Could not read "'//datafile//'".')   
    
    read(30, iostat = ios) isWritten
    if(ios /= 0) call throw_error("load_kernel_from_file",'Could not read "'//datafile//'".') 
    if(isWritten) then    
        allocate(Q%k1d(Q%Ndim))
        do i = 1,size(Q%k1d)
            read(30, iostat = ios) kernel_type
            if(ios /= 0) call throw_error("load_kernel_from_file",'Could not read "'//datafile//'".')   
            call Q%k1d(i)%init(kernel_type)
            if(Q%k1d(i)%Npar > 0) then
                read(30, iostat = ios) Q%k1d(i)%par    
                if(ios /= 0) call throw_error("load_kernel_from_file",'Could not read "'//datafile//'".')   
            end if
        end do 
    end if
    
    read(30, iostat = ios) isWritten
    if(ios /= 0) call throw_error("load_kernel_from_file",'Could not read "'//datafile//'".') 
    if(isWritten) then 
        allocate(Q%grid(Q%Ndim))
        do i = 1,size(Q%grid)
            read(30, iostat = ios) isWritten
            if(ios /= 0) call throw_error("load_kernel_from_file",'Could not read "'//datafile//'".')
            if(isWritten) then
                read(30, iostat = ios) grid_size
                if(ios /= 0) call throw_error("load_kernel_from_file",'Could not read "'//datafile//'".')   
                allocate(Q%grid(i)%x(grid_size))
                read(30, iostat = ios) Q%grid(i)%x
                if(ios /= 0) call throw_error("load_kernel_from_file",'Could not read "'//datafile//'".')   
            end if
        end do 
    end if
    
    read(30, iostat = ios) isWritten
    if(ios /= 0) call throw_error("load_kernel_from_file",'Could not read "'//datafile//'".')    
    if(isWritten) then      
        allocate(Q%values(Q%Npoint))
        read(30, iostat = ios) Q%values 
        if(ios /= 0) call throw_error("load_kernel_from_file",'Could not read "'//datafile//'".') 
    end if  
    
    read(30, iostat = ios) isWritten
    if(ios /= 0) call throw_error("load_kernel_from_file",'Could not read "'//datafile//'".')   
    if(isWritten) then   
        allocate(Q%valueIsPresent(Q%Npoint))
        read(30, iostat = ios) Q%valueIsPresent
        if(ios /= 0) call throw_error("load_kernel_from_file",'Could not read "'//datafile//'".')
    end if  
    
    read(30, iostat = ios) isWritten
    if(ios /= 0) call throw_error("load_kernel_from_file",'Could not read "'//datafile//'".')  
    if(isWritten) then    
        allocate(Q%alpha(Q%Npoint))
        read(30, iostat = ios) Q%alpha
        if(ios /= 0) call throw_error("load_kernel_from_file",'Could not read "'//datafile//'".')  
    end if 
    
    read(30, iostat = ios) isWritten
    if(ios /= 0) call throw_error("load_kernel_from_file",'Could not read "'//datafile//'".') 
    if(isWritten) then     
        allocate(Q%invQ(Q%Npoint,Q%Npoint))
        read(30, iostat = ios) Q%invQ
        if(ios /= 0) call throw_error("load_kernel_from_file",'Could not read "'//datafile//'".')
    end if 
    
    read(30, iostat = ios) isWritten
    if(ios /= 0) call throw_error("load_kernel_from_file",'Could not read "'//datafile//'".') 
    if(isWritten) then 
        call allocate_sigma_gamma(Q)
        read(30, iostat = ios) Q%sigma
        if(ios /= 0) call throw_error("load_kernel_from_file",'Could not read "'//datafile//'".')
    end if   
            
    close(io_unit)
    return 
end subroutine load_kernel_from_file

!-----------------------------------------------------------------------
!> @author
!> Oliver T. Unke, University of Basel
!
!> @brief
!> Frees kernel resources (deallocates memory)
!-----------------------------------------------------------------------
subroutine deallocate_kernel(Q)
    implicit none
    class(kernel) :: Q
    integer       :: i
    
    if(allocated(Q%k1d)) then
        do i = 1,size(Q%k1d)
            if(allocated(Q%k1d(i)%p2))   deallocate(Q%k1d(i)%p2)
            if(allocated(Q%k1d(i)%f2))   deallocate(Q%k1d(i)%f2)
            if(allocated(Q%k1d(i)%f3))   deallocate(Q%k1d(i)%f3)
            if(allocated(Q%k1d(i)%df2))  deallocate(Q%k1d(i)%df2)
            if(allocated(Q%k1d(i)%df3))  deallocate(Q%k1d(i)%df3)
            if(allocated(Q%k1d(i)%d2f2)) deallocate(Q%k1d(i)%d2f2)
            if(allocated(Q%k1d(i)%d2f3)) deallocate(Q%k1d(i)%d2f3)
            if(allocated(Q%k1d(i)%par))  deallocate(Q%k1d(i)%par)
        end do
        deallocate(Q%k1d)
    end if    
    if(allocated(Q%grid)) then
        do i = 1,size(Q%grid)
            if(allocated(Q%grid(i)%x)) deallocate(Q%grid(i)%x)
        end do
        deallocate(Q%grid)
    end if
    if(allocated(Q%values))         deallocate(Q%values)
    if(allocated(Q%valueIsPresent)) deallocate(Q%valueIsPresent)
    if(allocated(Q%alpha))          deallocate(Q%alpha)
    if(allocated(Q%invQ))           deallocate(Q%invQ)
    if(allocated(Q%sigma))          deallocate(Q%sigma)
    if(allocated(Q%gamma_old))      deallocate(Q%gamma_old)
    if(allocated(Q%gamma_new))      deallocate(Q%gamma_new)
    return 
end subroutine deallocate_kernel

!-----------------------------------------------------------------------
!> @author
!> Oliver T. Unke, University of Basel
!
!> @brief
!> Computes the index in the lookup table of precomputed sums for m and k
!
!> @details
!> Given a combination of values for m and k, the corresponding index in
!! the lookup table of precomputed sums is returned. This is necessary,
!! because the lookup table is a highly dimensional object (of unknown dimensions)
!! and is stored in a flattened 1-dimensional array
!-----------------------------------------------------------------------
function get_sig_idx_from_m_and_k(Q,Ndim,m,k) result(idx)
    implicit none
    class(kernel), intent(in)            :: Q
    integer, intent(in)                  :: Ndim !idx in how many dimensions
    integer, dimension(Ndim), intent(in) :: m,k
    integer                              :: idx, d, prod
    !NOTE: m here is a vector that contains either 0 or 1, it determines
    !whether a f2 or a f3 function is used. A 0 means the f2 function is used
    !and a 1 means the f3 function is used
    
    prod = 1
    idx  = 2*k(1)-m(1)
    do d = 2,Ndim
        prod = prod*(2*Q%k1d(d-1)%M2)
        idx = idx + (2*k(d)-m(d)-1)*prod
    end do
end function get_sig_idx_from_m_and_k

!-----------------------------------------------------------------------
!> @author
!> Oliver T. Unke, University of Basel
!
!> @brief
!> Computes the index for retrieving the correct lookup table for a certain index combination
!
!> @details
!> Given a combination of indices, the index for retrieving the correct lookup table
!! (containing all possible combinations of m and k) is returned. This is necessary
!! because the lookup table is stored in a flattened 1-dimensional array.
!-----------------------------------------------------------------------
function get_sig_idx_from_indices(Q,indices) result(idx)
    implicit none
    class(kernel), intent(in)              :: Q
    integer, dimension(Q%Ndim), intent(in) :: indices
    integer                                :: idx, d, prod
       
    prod = 1
    idx = indices(Q%Ndim) + 1
    do d = Q%Ndim-1,1,-1
        prod = prod*(size(Q%grid(d+1)%x)+1)
        idx = idx + indices(d) * prod
    end do  
end function get_sig_idx_from_indices

!-----------------------------------------------------------------------
!> @author
!> Oliver T. Unke, University of Basel
!
!> @brief
!> Computes the index for retrieving the correct alpha coefficient
!
!> @details
!> Given a combination of indices, the index for retrieving the correct alpha coefficient
!! is returned. This is necessary because the coefficients are stored in a flattened 1-dimensional array.
!-----------------------------------------------------------------------
function get_idx_from_indices(Q,indices) result(idx)
    implicit none
    class(kernel), intent(in)              :: Q
    integer, dimension(Q%Ndim), intent(in) :: indices
    integer                                :: idx, d, prod
       
    prod = 1
    idx = indices(Q%Ndim)
    do d = Q%Ndim-1,1,-1
        prod = prod*size(Q%grid(d+1)%x)
        idx = idx + (indices(d)-1) * prod
    end do  
end function get_idx_from_indices

!-----------------------------------------------------------------------
!> @author
!> Oliver T. Unke, University of Basel
!
!> @brief
!> Reads grid data from a .csv file and automatically allocates all necessary memory
!
!> @details
!> Given a combination of indices, the index for retrieving the correct alpha coefficient
!! is returned. This is necessary because the coefficients are stored in a flattened 1-dimensional array.
!! IMPORTANT:
!! The file must have the format x1,x2,x3,...,xn,f where the xi are the individual
!! coordinates and f the reference values. The file needs to be sorted 
!! in this order: x1>x2>x3>...>xn
! For example:
! x1_1, x2_1, x3_1 , ..., xn_1, f111...1
! x1_1, x2_1, x3_1 , ..., xn_2, f111...2
! x1_1, x2_1, x3_1 , ..., xn_3, f111...3
!   . ,   . ,   .  ,  . ,   . ,     .
!   . ,   . ,   .  ,  . ,   . ,     .
!   . ,   . ,   .  ,  . ,   . ,     .
! x1_n1,x2_n2,x3_n3,... ,xn_nn, fn1n2n3...nn
!-----------------------------------------------------------------------
subroutine read_grid_data(Q,datafile)
    implicit none
    character(len=*), intent(in)                 :: datafile
    class(kernel)                                :: Q
    character(len=1)                             :: dummy
    integer, dimension(:), allocatable           :: dimcount !counts how many points are in each dimension
    real(kind(0d0)), dimension(:),   allocatable :: currvalues
    real(kind(0d0)), dimension(:,:), allocatable :: gridvalues
    integer                                      :: alloc_status, ios, i, j
    
    !deallocate the kernel first
    call Q%free()
    
    !open grid datafile
    open(io_unit,file=datafile,status="old",action="read",iostat = ios)
    if(ios /= 0) call throw_error("read_grid_data", 'File "'//datafile//'" could not be opened.')
    
    !determine the number of dimensions/coordinates by counting the amount of ',' in the grid datafile
    i = 0
    do while(.true.)
        read(30,'(A1)',advance='no',iostat = ios) dummy
        if(ios /= 0) exit
        if(dummy == ',') i = i + 1
    end do
    rewind(30) !rewind file
    Q%Ndim = i
    
    !determine the total number of gridpoints
    i = 0
    do while(.true.)
        read(30,*,iostat = ios) dummy
        if(ios /= 0) exit
        i = i + 1
    end do
    rewind(30)
    Q%Npoint = i
     
    !allocate the required memory to the members of Q
    allocate(Q%k1d(Q%Ndim), Q%grid(Q%Ndim),Q%values(Q%Npoint),Q%alpha(Q%Npoint),&
             Q%valueIsPresent(Q%Npoint),stat=alloc_status)
    if(alloc_status /= 0) call throw_error("read_grid_data","could not allocate memory.")
       
    !allocate memory to variables used to find out grid dimensions
    allocate(dimcount(Q%Ndim), gridvalues(Q%Ndim,Q%Npoint), currvalues(Q%Ndim),stat=alloc_status)
    if(alloc_status /= 0) call throw_error("read_grid_data","could not allocate memory.")
    dimcount = 1
    !read the first line of gridvalues
    read(30,*,iostat=ios) gridvalues(:,1), Q%values(1)
    Q%valueIsPresent(1) = .not.isNaN(Q%values(1)) !detects holes
      
    if(ios /= 0) call throw_error("read_grid_data", 'File "'//datafile//'" could not be read properly.')
    !read the remaining grid values
    do i = 2,Q%Npoint
        read(30,*,iostat=ios) currvalues(:), Q%values(i)
        Q%valueIsPresent(i) = .not.isNaN(Q%values(i)) !detects holes
        if(ios /= 0) call throw_error("read_grid_data", 'File "'//datafile//'" could not be read properly.')
        do j = 1,Q%Ndim
            if(currvalues(j) > gridvalues(j,dimcount(j))) then
                dimcount(j) = dimcount(j) + 1
                gridvalues(j,dimcount(j)) = currvalues(j)
            end if
        end do
    end do    

    !store the individual coordinate grids in the kernel
    do i = 1,Q%Ndim
        allocate(Q%grid(i)%x(dimcount(i)))
        Q%grid(i)%x(:) = gridvalues(i,1:dimcount(i))
    end do
    
    !deallocate unused memory again
    deallocate(dimcount, gridvalues, currvalues)
    return  
end subroutine read_grid_data

!-----------------------------------------------------------------------
!> @author
!> Oliver T. Unke, University of Basel
!
!> @brief
!> Writes grid data as .csv file
!
!> @details
!> Writes the complete information from which the kernel was constructed into
!! a .csv file. This is useful to reconstruct the grid in case only the binary format
!! of the kernel is available.
!! The file has the format x1,x2,x3,...,xn,f where the xi are the individual
!! coordinates and f the reference values. The file is sorted in this order: x1>x2>x3>...>xn
! For example:
! x1_1, x2_1, x3_1 , ..., xn_1, f111...1
! x1_1, x2_1, x3_1 , ..., xn_2, f111...2
! x1_1, x2_1, x3_1 , ..., xn_3, f111...3
!   . ,   . ,   .  ,  . ,   . ,     .
!   . ,   . ,   .  ,  . ,   . ,     .
!   . ,   . ,   .  ,  . ,   . ,     .
! x1_n1,x2_n2,x3_n3,... ,xn_nn, fn1n2n3...nn
!-----------------------------------------------------------------------
subroutine write_grid_data(Q,datafile)
    implicit none
    character(len=*), intent(in) :: datafile
    class(kernel)                :: Q
    integer, dimension(Q%Ndim)   :: dimindx !for looping over the dimensions
    integer                      :: ios, counter, d
    
    !open grid datafile
    open(io_unit,file=datafile,status="replace",action="write",iostat = ios)
    if(ios /= 0) call throw_error("write_grid_data", 'File "'//datafile//'" could not be opened.')
    !write grid data to file    
    dimindx = 1
    counter = 0
    do while(dimindx(1) <= size(Q%grid(1)%x))        
        counter = counter + 1
        do d = 1,Q%Ndim
            write(io_unit,'(ES23.16,A1)',advance='no') Q%grid(d)%x(dimindx(d)),','
        end do        
        if(Q%valueIsPresent(counter)) then
            write(io_unit,'(ES23.16)') Q%values(counter)
        else
            write(io_unit,'(A23)') "NaN"
        end if
                
        !increase dimindx
        dimindx(Q%Ndim) = dimindx(Q%Ndim) + 1
        do d = Q%Ndim,2,-1
            if(dimindx(d) > size(Q%grid(d)%x)) then
                dimindx(d)   = 1
                dimindx(d-1) = dimindx(d-1) + 1
            end if
        end do
    end do        
    close(io_unit)
    return  
end subroutine write_grid_data

!-----------------------------------------------------------------------
!> @author
!> Oliver T. Unke, University of Basel
!
!> @brief
!> Calculates the kernel coefficients using the naive way (direct inversion of full matrix)
!
!> @details
!> Note that this method to calculate the kernel coefficients should never be used and instead, calculate_coefficients_fast 
!! should be called. The only exception is if the matrix is ill-conditioned, in which case the Tikhonov regularization 
!! procedure needs to be used, which is only possible with the slow method. Note that this subroutine also allows the
!! computation of the inverse of the kernel matrix by providing the additional argument calcInverse = .true.
!! This is needed for calculating error bounds if so desired (for large matrices, the calculation of the inverse
!! matrix is exceedingly slow).
!-----------------------------------------------------------------------
subroutine calculate_coefficients_slow(Q,alpha,calcInverse)
    implicit none
    class(kernel)                         :: Q
    real(kind(0d0)), dimension(&
        count(Q%valueIsPresent),&
        count(Q%valueIsPresent))          :: KM, invKM
    real(kind(0d0)), dimension(&
        count(Q%valueIsPresent),&
        Q%Ndim)                           :: fullgrid !matrix that stores the complete grid
    real(kind(0d0)), dimension(&
        count(Q%valueIsPresent))          :: coefficients, functionvalues
    real(kind(0d0)), optional, intent(in) :: alpha !regularization parameter
    integer, dimension(Q%Ndim)            :: dimindx !only needed for the traditional implementation
    integer                               :: i,j,k,counter,counter2,Npoint, alloc_status
    logical, optional, intent(in)         :: calcInverse !should the inverse be computed?
    
    !Npoint here stores the number of existing points
    Npoint = count(Q%valueIsPresent)
    
    !build the full grid
    counter  = 0 !keeps track of where we are in the total grid
    counter2 = 0 !keeps track of at which existing point we are
    dimindx = 1 !set dimindx = 1
    do while(dimindx(1) <= size(Q%grid(1)%x))        
        counter = counter + 1
        if(Q%valueIsPresent(counter)) then
            counter2 = counter2 + 1
            do i = 1,Q%Ndim
                fullgrid(counter2,i) = Q%grid(i)%x(dimindx(i))
            end do
        end if
                
        !increase dimindx
        dimindx(Q%Ndim) = dimindx(Q%Ndim) + 1
        do i = Q%Ndim,2,-1
            if(dimindx(i) > size(Q%grid(i)%x)) then
                dimindx(i)   = 1
                dimindx(i-1) = dimindx(i-1) + 1
            end if
        end do
    end do    
    
    !build full kernel matrix
    do i = 1,Npoint
        do j = 1,i
            KM(i,j) = 1d0
            !loop over all dimensions
            do k = 1,Q%Ndim
                KM(i,j) = KM(i,j)*Q%k1d(k)%k(fullgrid(i,k),fullgrid(j,k),Q%k1d(k)%par)
            end do
        end do
    end do
    
    !add the regularization parameter
    if(present(alpha)) then
        do i = 1,Npoint
            KM(i,i) = KM(i,i) + alpha
        end do
    end if
    
    !do cholesky decomposition and solve for coefficients
    functionvalues = pack(Q%values,Q%valueIsPresent) !only existing functionvalues
    call cholesky_decomposition(KM)
    call cholesky_solve(KM,coefficients,functionvalues)
    
    !store the coefficients 
    counter = 0
    do i = 1,Q%Npoint
        if(Q%valueIsPresent(i)) then
            counter = counter + 1
            Q%alpha(i) = coefficients(counter)
        else
            Q%alpha(i) = 0d0
        end if
    end do
    
    !also calculate the inverse of the kernel matrix
    !and store it (for calculating error bounds)
    if(.not.present(calcInverse)) return
    if(.not.calcInverse) return
        
    if(.not.allocated(Q%invQ)) then 
        allocate(Q%invQ(Q%Npoint,Q%Npoint),stat=alloc_status)
        if(alloc_status /= 0) call throw_error("calculate_coefficients_slow","could not allocate memory.")
        call cholesky_inverse(KM,invKM)
        counter = 0    
        do i = 1,Q%Npoint
            if(Q%valueIsPresent(i)) then
                counter = counter + 1
                counter2 = 0
                do j = 1,i
                   if(Q%valueIsPresent(j)) then
                       counter2 = counter2 + 1
                       Q%invQ(i,j) = invKM(counter,counter2)
                       Q%invQ(j,i) = Q%invQ(i,j)
                   else
                       Q%invQ(i,j) = 0d0
                   end if 
                end do
            else
                Q%invQ(i,:) = 0d0
            end if
        end do
    end if      
    return       
end subroutine calculate_coefficients_slow

!-----------------------------------------------------------------------
!> @author
!> Oliver T. Unke, University of Basel
!
!> @brief
!> Solves a system of linear equations in Tensor product form efficiently
!
!> @details
!> Solution contains the desired solution to the linear equations on input
!! and the calculated coefficients on output.
!! The algorithm is adapted from:
!! P. Fernandes and B. Plateau. "Triangular solution of linear systems in tensor product format." 
!! ACM SIGMETRICS Performance Evaluation Review 28.4 (2001): 30-32.
!-----------------------------------------------------------------------
subroutine fast_tensor_solve(Q,KM,solution)
    implicit none
    class(kernel)                                           :: Q
    type(kernel_matrix), dimension(Q%Ndim), intent(in)      :: KM
    real(kind(0d0)), dimension(size(Q%values)), intent(out) :: solution
    real(kind(0d0)), dimension(size(Q%values))              :: zin, zout    
    integer                                                 :: i,j,k,l,base,indx,nleft,nright
       
    !forward substitution
    !initialize nleft and nright (needed for the algorithm)
    nleft  = 1
    nright = 1
    do i = 1,Q%Ndim
        nright = nright * size(Q%grid(i)%x)
    end do
    
    !loop over all dimensions
    do i = 1,Q%Ndim      
        !calculate new nright
        nright = nright/size(Q%grid(i)%x)
        
        base = 0
        do k = 1,nleft
            do j = 1,nright
                indx = base + j
                do l = 1,size(Q%grid(i)%x)
                    zin(l) = solution(indx)
                    indx = indx + nright
                end do   
                call forward_substitution(KM(i)%M,zout(1:size(Q%grid(i)%x)),zin(1:size(Q%grid(i)%x)))
                indx = base + j
                do l = 1,size(Q%grid(i)%x)
                    solution(indx) = zout(l)
                    indx = indx + nright
                end do 
            end do
            base = base + (nright * size(Q%grid(i)%x))
        end do        
        !calculate new nleft
        nleft = nleft*size(Q%grid(i)%x)
    end do
    
    !backward substitution
    !initialize nleft and nright (needed for the algorithm)
    nleft  = 1
    nright = 1
    do i = 1,Q%Ndim
        nright = nright * size(Q%grid(i)%x)
    end do
    
    !loop over all dimensions
    do i = 1,Q%Ndim      
        !calculate new nright
        nright = nright/size(Q%grid(i)%x)
        
        base = 0
        do k = 1,nleft
            do j = 1,nright
                indx = base + j
                do l = 1,size(Q%grid(i)%x)
                    zin(l) = solution(indx)
                    indx = indx + nright
                end do   
                call backward_substitution(transpose(KM(i)%M),zout(1:size(Q%grid(i)%x)),zin(1:size(Q%grid(i)%x))) 
                indx = base + j
                do l = 1,size(Q%grid(i)%x)
                    solution(indx) = zout(l)
                    indx = indx + nright
                end do 
            end do
            base = base + (nright * size(Q%grid(i)%x))
        end do        
        !calculate new nleft
        nleft = nleft*size(Q%grid(i)%x)
    end do     
    return 
end subroutine fast_tensor_solve

!-----------------------------------------------------------------------
!> @author
!> Oliver T. Unke, University of Basel
!
!> @brief
!> Calculates the kernel coefficients using the fast way
!
!> @details
!> This subroutine calculates the kernel coefficients using the fact that we have a linear system in
!! tensor product form, therefore making the solution much easier. Note that this method to calculate the
!! kernel coefficients should always be prefered over calculate_coefficients_slow, unless the matrix is
!! ill conditioned. In that case the Tikhonov regularization procedure needs to be used, which is only
!! possible with the slow method. Note that it is not possible to calculate the error bound if the fast
!! method is used (because a calculation of the full inverse matrix is avoided - hence it's fast). 
!! If error bounds are desired, the slow method needs to be invoked. 
!! Holes in the grid are handled automatically.
!-----------------------------------------------------------------------
subroutine calculate_coefficients_fast(Q)
    implicit none
    class(kernel)                                    :: Q
    type(kernel)                                     :: KQ !temporary kernel for interpolations               
    type(kernel_matrix), dimension(Q%Ndim)           :: KM
    real(kind(0d0)), dimension(&
        count(.not.Q%valueIsPresent),&
        count(.not.Q%valueIsPresent))                :: cholQ
    real(kind(0d0)), dimension(&
        Q%Npoint,&
        count(.not.Q%valueIsPresent))                :: invQ !only a subslice of the full inverse matrix
    real(kind(0d0)), dimension(&
        count(.not.Q%valueIsPresent))                :: deltai
    integer, dimension(count(.not.Q%valueIsPresent)) :: indxNotPresent !keeps track of which values are missing
    real(kind(0d0)), dimension(Q%Npoint)             :: invQslice !slice of the full inverse matrix
    real(kind(0d0))                                  :: cdiff !for correcting coefficients when holes are present
    integer                                          :: i,j,k,counter,counteri,counterj,indxi,indxj,sizeDim
    integer                                          :: alloc_status
                 
    !build the kernel matrices for all dimensions and make their cholesky decomposition
    do i = 1,Q%Ndim
        !initialize
        if(allocated(KM(i)%M)) deallocate(KM(i)%M)
        allocate(KM(i)%M(size(Q%grid(i)%x),size(Q%grid(i)%x)),stat=alloc_status)
        if(alloc_status /= 0) call throw_error("calculate_coefficients_fast","could not allocate enough memory.") 
        
        !build the kernel matrices (only the lower diagonal part)
        do j = 1,size(Q%grid(i)%x)
            do k = 1,j
                KM(i)%M(j,k) = Q%k1d(i)%k(Q%grid(i)%x(j),Q%grid(i)%x(k),Q%k1d(i)%par)
            end do
        end do
        
        !do the Cholesky decomposition of the kernel matrices
        call cholesky_decomposition(KM(i)%M)
    end do
    
    !calculate coefficients
    if(.not.any(.not.Q%valueIsPresent)) then !straightforward if no holes present
        Q%alpha = Q%values !initialize alpha to solutions
        call fast_tensor_solve(Q,KM,Q%alpha) !solve for alpha
    else !a bit more sophisticated if holes are present
        !set up the temporary 1d kernel for constructing guess interpolations
        call KQ%free() !make sure nothing is allocated
        KQ%Ndim   = 1
        KQ%Npoint = size(Q%grid(Q%Ndim)%x)
        allocate(KQ%k1d(KQ%Ndim), KQ%grid(KQ%Ndim), KQ%values(KQ%Npoint), KQ%alpha(KQ%Npoint),&
                 KQ%valueIsPresent(KQ%Npoint), stat=alloc_status)
        if(alloc_status /= 0) call throw_error("calculate_coefficients_fast","could not allocate memory.")
        KQ%k1d(1) = Q%k1d(Q%Ndim)
        allocate(KQ%grid(1)%x(KQ%Npoint), stat=alloc_status)
        if(alloc_status /= 0) call throw_error("calculate_coefficients_fast","could not allocate memory.")
        KQ%grid(1)%x = Q%grid(Q%Ndim)%x
        
        !construct guesses for missing function values
        do i = 1,Q%Npoint,KQ%Npoint   
            !set up 1-d kernel for this cut
            KQ%valueIsPresent = Q%valueIsPresent(i:i+KQ%Npoint-1)
            if(count(KQ%valueIsPresent) > 0) then !if there is at least 1 point in the grid
                KQ%values         = Q%values(i:i+KQ%Npoint-1)    
                call KQ%calculate_coefficients_slow()
                call KQ%calculate_sums()
                do j = 1,KQ%Npoint
                    if(.not.KQ%valueIsPresent(j)) then
                        call KQ%evaluate_fast((/KQ%grid(1)%x(j)/),Q%values(i+j-1))
                    end if
                end do
            else !no point at all present -> no available information!
                Q%values(i:i+KQ%Npoint-1) = 0d0 
            end if
        end do
                
        !first calculate solution with "arbitrary" function values
        Q%alpha = Q%values !initialize alpha to solutions
        call fast_tensor_solve(Q,KM,Q%alpha)              
                 
        !build small submatrix
        counteri = 0
        !first find the indices of missing points
        do i = 1,Q%Npoint
            if(.not.Q%valueIsPresent(i)) then
                counteri = counteri + 1
                indxNotPresent(counteri) = i
            end if
        end do
        
        !construct submatrix of inverse (for all missing points)
        do i = 1,size(indxNotPresent)
            !compute a slice of the full inverse matrix
            invQ(:,i)                 = 0d0
            invQ(indxNotPresent(i),i) = 1d0
            call fast_tensor_solve(Q,KM,invQ(:,i)) 
            do j = 1,i !only the lower triangular part is needed
                cholQ(i,j) = invQ(indxNotPresent(j),i) !store part of the full slice          
            end do
        end do
                            
        !do cholesky decomposition of submatrix and solve for the unknown deltai
        call cholesky_decomposition(cholQ)
        call cholesky_solve(cholQ,deltai,pack(Q%alpha,.not.Q%valueIsPresent))
        
        !now that deltai is known, we can correct the coefficients
        counteri = 0
        do i = 1,Q%Npoint
            if(Q%valueIsPresent(i)) then
                counteri = counteri + 1
                cdiff = 0d0
                do j = 1,size(indxNotPresent)
                    cdiff = cdiff + deltai(j)*invQ(i,j)        
                end do                  
                !correct the coefficients
                Q%alpha(i) =  Q%alpha(i) - cdiff  
            else !points that are not present MUST have a coefficient of 0
                Q%alpha(i) = 0d0
            end if  
        end do     
        call KQ%free() !free temporary 1d kernel            
    end if 
     
    !deallocate memory again
    do i = 1,Q%Ndim
        if(allocated(KM(i)%M))    deallocate(KM(i)%M)
    end do  
    return
end subroutine calculate_coefficients_fast

!-----------------------------------------------------------------------
!> @author
!> Oliver T. Unke, University of Basel
!
!> @brief
!> Computes the presummed values needed for the fast kernel evaluation for a single point
!
!> @details
!> This subroutine calculates all sums for a given combination of m and k
!! needed for the fast evaluation of a kernel. Note that the sums are
!! calculated using a recurrence relation and the complexity of computing
!! all sums therefore is O(Npoint) instead of O(Npoint^2) (naive implementaion).
!-----------------------------------------------------------------------
subroutine calculate_single_sum(Q,sigma,m,k)
    implicit none
    type(kernel)                                      :: Q
    real(kind(0d0)), dimension(Q%Npoint), intent(out) :: sigma
    integer, dimension(Q%Ndim), intent(in)            :: m,k
    type(f_ptr), dimension(Q%Ndim)                    :: fmk !array of function pointers
    real(kind(0d0)), dimension(Q%Ndim)                :: p2k !array of p values
    real(kind(0d0))                                   :: prod, Sjk  !intermediate product and intermediate sum
    integer, dimension(Q%Ndim)                        :: z, z_start, z_end, z_incr !keeps track of z values
    integer, dimension(Q%Ndim)                        :: zvec, ivec !index vectors
    integer, dimension(Q%Ndim)                        :: loop_end, loop_idx
    integer                                           :: d,idx,j,pos1,pos2,pos3,counter
    logical                                           :: sum_exists
    
    sigma = 0d0 !set everything to 0
       
    !build the array of function pointers and p values for this combination of m and k
    !and also initialize loop_start and loop_end  for each dimension (depends on m),
    !as well as decide how to loop over the z values (also depends on m)
    !NOTE: m = 0 encodes f2 functions, m=1 encodes f3 functions
    do d = 1,Q%Ndim
        p2k(d) = Q%k1d(d)%p2(k(d)) !p2k is independent of m
        if(m(d) == 0) then !this is like m=2 in the paper (f2k functions)
            fmk(d)%f => Q%k1d(d)%f2(k(d))%f !assign to f2k function
            z_start(d) = 0 !loop should start with z(d) = 0
            z_end(d)   = size(Q%grid(d)%x) !loop should end with z(d) = Nd
            z_incr(d)  = 1  
        else               !this is like m=3 in the paper (f3k functions)
            fmk(d)%f => Q%k1d(d)%f3(k(d))%f !assign to f3k function
            z_start(d) = size(Q%grid(d)%x) !loop should start with z(d) = Nd
            z_end(d)   = 0 !loop should end with z(d) = 0
            z_incr(d)  = -1                        
        end if  
        loop_end(d) = size(Q%grid(d)%x) !loop end is always N
    end do  
    
    !initialize z values
    z = z_start
    
    !loop over all possible grid points
    do while(z(1) /= (z_end(1) + z_incr(1))) 
        idx        = Q%get_lookup_idx(z)
        sigma(idx) = 0d0 !initialize sum to 0
        
        !assign appropriate loop_idx, depending on z and m
        do d = 1,Q%Ndim
            if(m(d) == 0) then !this is like m=2 in the paper (f2k functions)
                loop_idx(d) = z(d)
            else !this is like m=3 in the paper (f3k functions)
                loop_idx(d) = z(d)+1
            end if    
        end do
        
        !initialize to the value at the current point (if meaningful, else it stays at 0)
        if(.not.any(loop_idx > loop_end) .and. .not.any(loop_idx < 1)) then
            prod = 1d0 !initialize product
            do d = 1,Q%Ndim
                prod = prod * p2k(d)*fmk(d)%f(Q%grid(d)%x(loop_idx(d)),Q%k1d(d)%par)
            end do
            sigma(idx) = sigma(idx) + Q%alpha(Q%get_alpha_idx(loop_idx)) * prod !add to sum
        end if
        
        !loop over all possible values of j
        do j = 1,Q%Ndim
            !initialize intermediate sum to 0
            Sjk = 0d0
            
            !loop over the set of vectors containing j ones and Ndim-j zeros
            counter = 0
            pos1 = 1
            pos2 = pos1+j-1
            pos3 = pos2+1
            ivec = 0
            ivec(pos1:pos2) = 1
            do 
                !build the z vector (depends on what values of m there are)
                sum_exists = .true.
                do d = 1,Q%Ndim
                    if(m(d) == 0) then !m=2, ivec gets subtracted
                        zvec(d) = z(d) - ivec(d)
                        if(zvec(d) < 0) then
                            sum_exists = .false.
                            exit
                        end if
                    else !m=3, ivec gets added
                        zvec(d) = z(d) + ivec(d)
                        if(zvec(d) > size(Q%grid(d)%x)) then
                            sum_exists = .false.
                            exit
                        end if
                    end if
                end do
                !add to Sjk (if that value even exists)   
                if(sum_exists) then
                    Sjk = Sjk + sigma(Q%get_lookup_idx(zvec))
                end if

                !everything below is just to make sure that every possible
                !ivec is exactly used once!
                if(.not.any(ivec(Q%Ndim-j+1:Q%Ndim) /= 1)) exit
                
                if(pos3 <= Q%Ndim-counter) then 
                    ivec(pos3)   = 1
                    ivec(pos3-1) = 0
                    pos3         = pos3 + 1
                else
                    counter = counter + 1
                    pos2       = pos2-1
                    ivec(pos2) = 0
                    pos3       = pos2+1
                    
                    if(pos3 <= Q%Ndim-counter) then
                        ivec(pos3)   = 1
                        ivec(pos3-1) = 0
                        pos3         = pos3 + 1
                    end if     
                end if
                                
                if(pos2 <= pos1) then
                    counter = 0
                    pos1 = pos1+1
                    pos2 = pos1+j-1
                    pos3 = pos2+1
                    ivec = 0
                    ivec(pos1:pos2) = 1
                end if
            end do
                        
            if(mod(j,2) == 1) then
                sigma(idx) = sigma(idx) + Sjk 
            else
                sigma(idx) = sigma(idx) - Sjk 
            end if        
        end do
        
        !increase z successively (loops over all possible grid points this way)
        z(Q%Ndim) = z(Q%Ndim) + z_incr(Q%Ndim)
        do d = Q%Ndim,2,-1
            if(z(d) == (z_end(d) + z_incr(d))) then
                z(d)   = z_start(d)
                z(d-1) = z(d-1) + z_incr(d-1)
            end if
        end do
    end do    
    return
end subroutine calculate_single_sum

!-----------------------------------------------------------------------
!> @author
!> Oliver T. Unke, University of Basel
!
!> @brief
!> Allocates memory for the lookup table
!
!> @details
!> Since the amount of memory for storing the lookup table can get quite
!! large when an extremely high-dimensional grid with many points is used,
!! this subroutine will print the approximate RAM requirement in case
!! the allocation fails.
!-----------------------------------------------------------------------
subroutine allocate_sigma_gamma(Q)
    type(kernel)        :: Q
    integer             :: alloc_status
    character(len=1024) :: alloc_err_msg
    integer             :: prod,d
    
    if(allocated(Q%sigma))     deallocate(Q%sigma)
    if(allocated(Q%gamma_old)) deallocate(Q%gamma_old)
    if(allocated(Q%gamma_new)) deallocate(Q%gamma_new)
    
    !determine how much memory to allocate to sigma
    prod = 1
    do d = 1,Q%Ndim
        prod = prod * (size(Q%grid(d)%x)+1)
    end do
    
    !allocate memory for the storage of the sums
    allocate(Q%sigma(2**Q%Ndim*product(Q%k1d(:)%M2),prod),&
             Q%gamma_old(2**Q%Ndim*product(Q%k1d(:)%M2)),&
             Q%gamma_new(2**Q%Ndim*product(Q%k1d(:)%M2)), stat=alloc_status)
    if(alloc_status /= 0) then
        write(alloc_err_msg,'(I0)') sizeof(0d0)*2**Q%Ndim*product(Q%k1d(:)%M2)*(prod+2)
        call throw_error("allocate_sigma_gamma","Could not allocate enough memory for storing precomputed sums. "//&
                                          "Estimated RAM requirement: "//trim(alloc_err_msg)//" bytes.")
    end if
    return
end subroutine allocate_sigma_gamma

!-----------------------------------------------------------------------
!> @author
!> Oliver T. Unke, University of Basel
!
!> @brief
!> Computes the presummed values needed for the fast kernel evaluation
!
!> @details
!> This subroutine calculates all sums needed for the fast evaluation of a kernel. 
!! Note that the sums are calculated using a recurrence relation and the complexity of computing
!! all sums therefore is O(Npoint) instead of O(Npoint^2) (naive implementaion).
!-----------------------------------------------------------------------
subroutine calculate_sums(Q)
    implicit none
    class(kernel)              :: Q
    integer, dimension(Q%Ndim) :: m,k
    integer                    :: d, idx
    
    call allocate_sigma_gamma(Q)
    
    !loop over all possible combinations of m and k
    m = 0
    k = 1
    do while(k(1) <= Q%k1d(1)%M2)
        m = 0
        do while(m(1) <= 1)
            !find appropriate index for this combination of m and k
            idx = Q%get_sig_idx(Q%Ndim,m,k) 
                                    
            !do actual summation for this combination of m and k
            call calculate_single_sum(Q,Q%sigma(idx,:),m,k)  
            
            !increase m
            m(Q%Ndim) = m(Q%Ndim) + 1
            do d = Q%Ndim,2,-1
                if(m(d) > 1) then
                    m(d)   = 0
                    m(d-1) = m(d-1) + 1
                end if
            end do
        end do
        !increase k
        k(Q%Ndim) = k(Q%Ndim) + 1
        do d = Q%Ndim,2,-1
            if(k(d) > Q%k1d(d)%M2) then
                k(d)   = 1
                k(d-1) = k(d-1) + 1
            end if
        end do
    end do
    return 
end subroutine calculate_sums

!-----------------------------------------------------------------------
!> @author
!> Oliver T. Unke, University of Basel
!
!> @brief
!> Binary search algorithm to find correct indices for lookup tables
!
!> @details
!> Given an ascendingly sorted array of values and a search value, this function returns the
!! index i of the element xi in array which satisfies xi <= value < xi+1.
!! In the case that the smallest element is smaller than value, the function returns 0.
!-----------------------------------------------------------------------
integer function binary_search(array,value)
    implicit none
    real(kind(0d0)), dimension(:), intent(in) :: array
    real(kind(0d0)), intent(in)               :: value
    integer                                   :: l,r,m
    
    !if value is smaller or larger than the contents of the array, we can immediately leave
    if (array(1) > value) then
        binary_search = 0
        return
    else if (array(size(array)) < value) then
        binary_search = size(array)
        return
    end if
    
    l = 1
    r = size(array)
    do while(l <= r)
        m = (l+r)/2
        if(array(m) > value) then
            r = m-1
        else
            l = m+1
        end if
    end do
    binary_search = (l+r)/2
end function binary_search

!-----------------------------------------------------------------------
!> @author
!> Oliver T. Unke, University of Basel
!
!> @brief
!> Evaluates the kernel using the fast algorithm and returns the function value (and derivatives).
!
!> @details
!> This method should always be preferred over evaluate_kernel_slow, unless the memory requirement
!! for the lookup table is larger than the available RAM. In that case, the only way
!! to evaluate the kernel is to use the slow algorithm. For very small grids, the slow evaluation
!! might actually be slightly faster (this should be benchmarked in performance critical applications).
!-----------------------------------------------------------------------
subroutine evaluate_kernel_fast(Q,x,f,df,d2f,df_mask,d2f_mask)
    implicit none
    class(kernel)                                                      :: Q
    real(kind(0d0)), dimension(:), intent(in)                          :: x !point at which to evaluate
    real(kind(0d0)), intent(out)                                       :: f !function value at the evaluation point
    real(kind(0d0)), dimension(size(x)), intent(out), optional         :: df !partial derivatives at that point
    real(kind(0d0)), dimension(size(x),size(x)), intent(out), optional :: d2f !Hessian at that point
    logical,         dimension(size(x)), intent(in),  optional         :: df_mask !for calculating only some derivatives
    logical,         dimension(size(x),size(x)), intent(in),  optional :: d2f_mask !for calculating only some entries
    integer, dimension(Q%Ndim)                                         :: m,k
    integer, dimension(Q%Ndim)                                         :: z !index array
    integer                                                            :: deriv,deriv2,d,i,kp1,new_idx,old_idx
    
    !find appropriate indices to initialize gamma from lookup table
    do d = 1,Q%Ndim
        z(d) = binary_search(Q%grid(d)%x,x(d))  
    end do
    
    !initialize to sigma from lookup table
    Q%gamma_new = Q%sigma(:,Q%get_lookup_idx(z))
    Q%gamma_old = Q%gamma_new
    
    !iteratively update gamma values  
    do d = Q%Ndim-1,1,-1 !iteratively update gamma values until arriving at dimension 1 
        !loop over all possible k and m values
        k = 1
        do while(k(1) <= Q%k1d(1)%M2)
            m = 0
            do while(m(1) <= 1)
                new_idx = Q%get_sig_idx(d,m(1:d),k(1:d)) !get the appropriate index for the gamma_new array
                Q%gamma_new(new_idx) = 0d0 !initialize to 0
                
                !loop over possible k values in dimension d+1
                do kp1 = 1,Q%k1d(d+1)%M2
                    k(d+1) = kp1
                       
                    !first, calculate the f3k part (uses gamma2k)
                    m(d+1) = 0
                    old_idx = Q%get_sig_idx(d+1,m(1:d+1),k(1:d+1))
                    Q%gamma_new(new_idx) = Q%gamma_new(new_idx) &
                                         + Q%gamma_old(old_idx) &
                                         * Q%k1d(d+1)%f3(kp1)%f(x(d+1),Q%k1d(d+1)%par)
                    
                    !second, calculate the f2k part (uses gamma3k)
                    m(d+1) = 1
                    old_idx = Q%get_sig_idx(d+1,m(1:d+1),k(1:d+1))
                    Q%gamma_new(new_idx) = Q%gamma_new(new_idx) &
                                         + Q%gamma_old(old_idx) &
                                         * Q%k1d(d+1)%f2(kp1)%f(x(d+1),Q%k1d(d+1)%par)                              
                end do
                                      
                !increment m
                m(d) = m(d) + 1
                do i = d,2,-1
                    if(m(i) > 1) then
                        m(i)   = 0
                        m(i-1) = m(i-1) + 1
                    end if
                end do
            end do
            !increment k
            k(d) = k(d) + 1
            do i = d,2,-1
                if(k(i) > Q%k1d(i)%M2) then
                    k(i)   = 1
                    k(i-1) = k(i-1) + 1
                end if
            end do
        end do  
        
        !update old gamma for next iteration, but only the values that are actually used
        kp1 = 2**d*product(Q%k1d(1:d)%M2) 
        Q%gamma_old(1:kp1) = Q%gamma_new(1:kp1)
    end do
            
    !Now, finally, the gamma value is updated and the energy can be computed
    f = 0d0 !initialize to 0
    !loop over possible k values in dimension 1
    do kp1 = 1,Q%k1d(1)%M2
        !first, calculate the f3k part (uses gamma2k)        
        new_idx = Q%get_sig_idx(1,(/0/),(/kp1/))
        f = f + Q%gamma_new(new_idx) &
              * Q%k1d(1)%f3(kp1)%f(x(1),Q%k1d(1)%par)
        
        !second, calculate the f2k part (uses gamma3k)
        new_idx = Q%get_sig_idx(1,(/1/),(/kp1/))
        f = f + Q%gamma_new(new_idx) &
              * Q%k1d(1)%f2(kp1)%f(x(1),Q%k1d(1)%par)                                
    end do
    
    !also compute partial derivatives if needed    
    if(present(df)) then
        !loop over all possible partial derivatives
        do deriv = 1,Q%Ndim  
            !skip values if a mask is provided
            if(present(df_mask)) then
                if(.not.df_mask(deriv)) cycle
            end if                 
            Q%gamma_old = Q%sigma(:,Q%get_lookup_idx(z)) !initialize to sigma from lookup table
            !iteratively update gamma values
            do d = Q%Ndim-1,1,-1 !iteratively update gamma values until arriving at dimension 1      
                !loop over all possible k and m values
                k = 1
                do while(k(1) <= Q%k1d(1)%M2)
                    m = 0
                    do while(m(1) <= 1)
                        new_idx = Q%get_sig_idx(d,m(1:d),k(1:d)) !get the appropriate index for the gamma_new array
                        Q%gamma_new(new_idx) = 0d0 !initialize to 0
                        
                        !loop over possible k values in dimension d+1
                        do kp1 = 1,Q%k1d(d+1)%M2
                            k(d+1) = kp1
                            
                            if(deriv == d+1) then
                                !first, calculate the f3k part (uses gamma2k)
                                m(d+1) = 0
                                old_idx = Q%get_sig_idx(d+1,m(1:d+1),k(1:d+1))
                                Q%gamma_new(new_idx) = Q%gamma_new(new_idx) &
                                                     + Q%gamma_old(old_idx) &
                                                     * Q%k1d(d+1)%df3(kp1)%f(x(d+1),Q%k1d(d+1)%par)
                                
                                !second, calculate the f2k part (uses gamma3k)
                                m(d+1) = 1
                                old_idx = Q%get_sig_idx(d+1,m(1:d+1),k(1:d+1))
                                Q%gamma_new(new_idx) = Q%gamma_new(new_idx) &
                                                     + Q%gamma_old(old_idx) &
                                                     * Q%k1d(d+1)%df2(kp1)%f(x(d+1),Q%k1d(d+1)%par)   
                            else
                                !first, calculate the f3k part (uses gamma2k)
                                m(d+1) = 0
                                old_idx = Q%get_sig_idx(d+1,m(1:d+1),k(1:d+1))
                                Q%gamma_new(new_idx) = Q%gamma_new(new_idx) &
                                                     + Q%gamma_old(old_idx) &
                                                     * Q%k1d(d+1)%f3(kp1)%f(x(d+1),Q%k1d(d+1)%par)
                                
                                !second, calculate the f2k part (uses gamma3k)
                                m(d+1) = 1
                                old_idx = Q%get_sig_idx(d+1,m(1:d+1),k(1:d+1))
                                Q%gamma_new(new_idx) = Q%gamma_new(new_idx) &
                                                     + Q%gamma_old(old_idx) &
                                                     * Q%k1d(d+1)%f2(kp1)%f(x(d+1),Q%k1d(d+1)%par)   
                            end if                                                      
                        end do
                                              
                        !increment m
                        m(d) = m(d) + 1
                        do i = d,2,-1
                            if(m(i) > 1) then
                                m(i)   = 0
                                m(i-1) = m(i-1) + 1
                            end if
                        end do
                    end do
                    !increment k
                    k(d) = k(d) + 1
                    do i = d,2,-1
                        if(k(i) > Q%k1d(i)%M2) then
                            k(i)   = 1
                            k(i-1) = k(i-1) + 1
                        end if
                    end do
                end do  
                
                !update old gamma for next iteration, but only the values that are actually used
                kp1 = 2**d*product(Q%k1d(1:d)%M2) 
                Q%gamma_old(1:kp1) = Q%gamma_new(1:kp1)
            end do
                
            !Now, finally, the gamma value is updated and the energy can be computed
            df(deriv) = 0d0 !initialize to 0
            !loop over possible k values in dimension 1
            do kp1 = 1,Q%k1d(1)%M2
                if(deriv == 1) then
                    !first, calculate the f3k part (uses gamma2k)        
                    new_idx = Q%get_sig_idx(1,(/0/),(/kp1/))
                    df(deriv) = df(deriv) + Q%gamma_new(new_idx) &
                                          * Q%k1d(1)%df3(kp1)%f(x(1),Q%k1d(1)%par)
                    
                    !second, calculate the f2k part (uses gamma3k)
                    new_idx = Q%get_sig_idx(1,(/1/),(/kp1/))
                    df(deriv) = df(deriv) + Q%gamma_new(new_idx) &
                                          * Q%k1d(1)%df2(kp1)%f(x(1),Q%k1d(1)%par)                  
                else
                    !first, calculate the f3k part (uses gamma2k)        
                    new_idx = Q%get_sig_idx(1,(/0/),(/kp1/))
                    df(deriv) = df(deriv) + Q%gamma_new(new_idx) &
                                          * Q%k1d(1)%f3(kp1)%f(x(1),Q%k1d(1)%par)
                    
                    !second, calculate the f2k part (uses gamma3k)
                    new_idx = Q%get_sig_idx(1,(/1/),(/kp1/))
                    df(deriv) = df(deriv) + Q%gamma_new(new_idx) &
                                          * Q%k1d(1)%f2(kp1)%f(x(1),Q%k1d(1)%par)  
                end if                              
            end do
        end do
    end if 
    
    !also compute Hessian if needed    
    if(present(d2f)) then
        !loop over all possible entries of the Hessian
        do deriv = 1,Q%Ndim  
            do deriv2 = 1,Q%Ndim  
                !skip values if a mask is provided
                if(present(d2f_mask)) then
                    if(.not.d2f_mask(deriv,deriv2)) cycle
                end if     
                !because the Hessian is symmetric, we don't need to recalculate values
                !that were already calculated
                if(deriv2 > deriv) then
                    if(present(d2f_mask)) then
                        if(d2f_mask(deriv2,deriv)) then !we can only do this if the value was calculated before
                            d2f(deriv,deriv2) = d2f(deriv2,deriv)
                            cycle
                        end if
                    else
                        d2f(deriv,deriv2) = d2f(deriv2,deriv)
                        cycle
                    end if
                end if
                                          
                Q%gamma_old = Q%sigma(:,Q%get_lookup_idx(z)) !initialize to sigma from lookup table
                !iteratively update gamma values
                do d = Q%Ndim-1,1,-1 !iteratively update gamma values until arriving at dimension 1      
                    !loop over all possible k and m values
                    k = 1
                    do while(k(1) <= Q%k1d(1)%M2)
                        m = 0
                        do while(m(1) <= 1)
                            new_idx = Q%get_sig_idx(d,m(1:d),k(1:d)) !get the appropriate index for the gamma_new array
                            Q%gamma_new(new_idx) = 0d0 !initialize to 0
                            
                            !loop over possible k values in dimension d+1
                            do kp1 = 1,Q%k1d(d+1)%M2
                                k(d+1) = kp1                               
                                if(deriv == d+1 .and. deriv2 == d+1) then !second derivative needed
                                    !first, calculate the f3k part (uses gamma2k)
                                    m(d+1) = 0
                                    old_idx = Q%get_sig_idx(d+1,m(1:d+1),k(1:d+1))
                                    Q%gamma_new(new_idx) = Q%gamma_new(new_idx) &
                                                         + Q%gamma_old(old_idx) &
                                                         * Q%k1d(d+1)%d2f3(kp1)%f(x(d+1),Q%k1d(d+1)%par)
                                    
                                    !second, calculate the f2k part (uses gamma3k)
                                    m(d+1) = 1
                                    old_idx = Q%get_sig_idx(d+1,m(1:d+1),k(1:d+1))
                                    Q%gamma_new(new_idx) = Q%gamma_new(new_idx) &
                                                         + Q%gamma_old(old_idx) &
                                                         * Q%k1d(d+1)%d2f2(kp1)%f(x(d+1),Q%k1d(d+1)%par)   
                                else if(deriv == d+1 .or. deriv2 == d+1) then !first derivative needed
                                    !first, calculate the f3k part (uses gamma2k)
                                    m(d+1) = 0
                                    old_idx = Q%get_sig_idx(d+1,m(1:d+1),k(1:d+1))
                                    Q%gamma_new(new_idx) = Q%gamma_new(new_idx) &
                                                         + Q%gamma_old(old_idx) &
                                                         * Q%k1d(d+1)%df3(kp1)%f(x(d+1),Q%k1d(d+1)%par)
                                    
                                    !second, calculate the f2k part (uses gamma3k)
                                    m(d+1) = 1
                                    old_idx = Q%get_sig_idx(d+1,m(1:d+1),k(1:d+1))
                                    Q%gamma_new(new_idx) = Q%gamma_new(new_idx) &
                                                         + Q%gamma_old(old_idx) &
                                                         * Q%k1d(d+1)%df2(kp1)%f(x(d+1),Q%k1d(d+1)%par)   
                                else
                                    !first, calculate the f3k part (uses gamma2k)
                                    m(d+1) = 0
                                    old_idx = Q%get_sig_idx(d+1,m(1:d+1),k(1:d+1))
                                    Q%gamma_new(new_idx) = Q%gamma_new(new_idx) &
                                                         + Q%gamma_old(old_idx) &
                                                         * Q%k1d(d+1)%f3(kp1)%f(x(d+1),Q%k1d(d+1)%par)
                                    
                                    !second, calculate the f2k part (uses gamma3k)
                                    m(d+1) = 1
                                    old_idx = Q%get_sig_idx(d+1,m(1:d+1),k(1:d+1))
                                    Q%gamma_new(new_idx) = Q%gamma_new(new_idx) &
                                                         + Q%gamma_old(old_idx) &
                                                         * Q%k1d(d+1)%f2(kp1)%f(x(d+1),Q%k1d(d+1)%par)   
                                end if                                                      
                            end do
                                                  
                            !increment m
                            m(d) = m(d) + 1
                            do i = d,2,-1
                                if(m(i) > 1) then
                                    m(i)   = 0
                                    m(i-1) = m(i-1) + 1
                                end if
                            end do
                        end do
                        !increment k
                        k(d) = k(d) + 1
                        do i = d,2,-1
                            if(k(i) > Q%k1d(i)%M2) then
                                k(i)   = 1
                                k(i-1) = k(i-1) + 1
                            end if
                        end do
                    end do  
                    
                    !update old gamma for next iteration, but only the values that are actually used
                    kp1 = 2**d*product(Q%k1d(1:d)%M2) 
                    Q%gamma_old(1:kp1) = Q%gamma_new(1:kp1)
                end do
                    
                !Now, finally, the gamma value is updated and the energy can be computed
                d2f(deriv,deriv2) = 0d0 !initialize to 0
                !loop over possible k values in dimension 1
                do kp1 = 1,Q%k1d(1)%M2
                    if(deriv == 1 .and. deriv2 == 1) then !second derivative is needed
                        !first, calculate the f3k part (uses gamma2k)        
                        new_idx = Q%get_sig_idx(1,(/0/),(/kp1/))
                        d2f(deriv,deriv2) = d2f(deriv,deriv2) + Q%gamma_new(new_idx) &
                                            * Q%k1d(1)%d2f3(kp1)%f(x(1),Q%k1d(1)%par)
                        
                        !second, calculate the f2k part (uses gamma3k)
                        new_idx = Q%get_sig_idx(1,(/1/),(/kp1/))
                        d2f(deriv,deriv2) = d2f(deriv,deriv2) + Q%gamma_new(new_idx) &
                                            * Q%k1d(1)%d2f2(kp1)%f(x(1),Q%k1d(1)%par)                     
                    else if(deriv == 1 .or. deriv2 == 1) then !first derivative is needed
                        !first, calculate the f3k part (uses gamma2k)        
                        new_idx = Q%get_sig_idx(1,(/0/),(/kp1/))
                        d2f(deriv,deriv2) = d2f(deriv,deriv2) + Q%gamma_new(new_idx) &
                                            * Q%k1d(1)%df3(kp1)%f(x(1),Q%k1d(1)%par)
                        
                        !second, calculate the f2k part (uses gamma3k)
                        new_idx = Q%get_sig_idx(1,(/1/),(/kp1/))
                        d2f(deriv,deriv2) = d2f(deriv,deriv2) + Q%gamma_new(new_idx) &
                                            * Q%k1d(1)%df2(kp1)%f(x(1),Q%k1d(1)%par)                  
                    else
                        !first, calculate the f3k part (uses gamma2k)        
                        new_idx = Q%get_sig_idx(1,(/0/),(/kp1/))
                        d2f(deriv,deriv2) = d2f(deriv,deriv2) + Q%gamma_new(new_idx) &
                                            * Q%k1d(1)%f3(kp1)%f(x(1),Q%k1d(1)%par)
                        
                        !second, calculate the f2k part (uses gamma3k)
                        new_idx = Q%get_sig_idx(1,(/1/),(/kp1/))
                        d2f(deriv,deriv2) = d2f(deriv,deriv2) + Q%gamma_new(new_idx) &
                                            * Q%k1d(1)%f2(kp1)%f(x(1),Q%k1d(1)%par)  
                    end if                              
                end do
            end do
        end do
    end if 
    return 
end subroutine evaluate_kernel_fast

!-----------------------------------------------------------------------
!> @author
!> Oliver T. Unke, University of Basel
!
!> @brief
!> Evaluates the kernel using the naive (slow) algorithm and returns the function value (and derivatives).
!
!> @details
!> This method should never be called and instead, evaluate_kernel_fast should be used. However,
!! if not enough RAM is available to store the lookup table of precomputed sums, this may be the only
!! way to evaluate the kernel at all.
!-----------------------------------------------------------------------
subroutine evaluate_kernel_slow(Q,x,f,df,d2f,df_mask,d2f_mask)
    implicit none
    class(kernel)                                                      :: Q
    real(kind(0d0)), dimension(:), intent(in)                          :: x !point at which to evaluate
    real(kind(0d0)), intent(out)                                       :: f !function value at the evaluation point
    real(kind(0d0)), dimension(size(x)), intent(out), optional         :: df !partial derivatives at that point
    real(kind(0d0)), dimension(size(x),size(x)), intent(out), optional :: d2f !Hessian at that point
    logical,         dimension(size(x)), intent(in),  optional         :: df_mask !for calculating only some derivatives
    logical,         dimension(size(x),size(x)), intent(in),  optional :: d2f_mask !for calculating only some entries
    integer, dimension(size(x))                                        :: dimindx
    real(kind(0d0))                                                    :: prod
    integer                                                            :: i, d, d2, counter
    
    f = 0d0
    counter = 0
    dimindx = 1 !set dimindx = 1
    do while(dimindx(1) <= size(Q%grid(1)%x))
        counter = counter + 1
        prod = 1d0
        do i = 1,Q%Ndim
            prod = prod*Q%k1d(i)%k(x(i),Q%grid(i)%x(dimindx(i)),Q%k1d(i)%par)
        end do

        f = f + Q%alpha(counter)*prod
                
        !increase dimindx
        dimindx(Q%Ndim) = dimindx(Q%Ndim) + 1
        do i = Q%Ndim,2,-1
            if(dimindx(i) > size(Q%grid(i)%x)) then
                dimindx(i)   = 1
                dimindx(i-1) = dimindx(i-1) + 1
            end if
        end do
    end do
    
    if(present(df)) then
        do d = 1,Q%Ndim
            !skip values if a mask is provided
            if(present(df_mask)) then
                if(.not.df_mask(d)) cycle
            end if      
            df(d) = 0d0
            counter = 0
            dimindx = 1
            do while(dimindx(1) <= size(Q%grid(1)%x))
                counter = counter + 1
                prod = 1d0
                do i = 1,Q%Ndim
                    if(i == d) then
                        prod = prod*Q%k1d(i)%dk(x(i),Q%grid(i)%x(dimindx(i)),Q%k1d(i)%par)
                    else
                        prod = prod*Q%k1d(i)%k(x(i),Q%grid(i)%x(dimindx(i)),Q%k1d(i)%par)
                    end if
                end do
                df(d) = df(d) + Q%alpha(counter)*prod
                        
                !increase dimindx
                dimindx(Q%Ndim) = dimindx(Q%Ndim) + 1
                do i = Q%Ndim,2,-1
                    if(dimindx(i) > size(Q%grid(i)%x)) then
                        dimindx(i)   = 1
                        dimindx(i-1) = dimindx(i-1) + 1
                    end if
                end do
            end do    
        end do
    end if
    
    if(present(d2f)) then
        do d = 1,Q%Ndim
            do d2 = 1,Q%Ndim
                !skip values if a mask is provided
                if(present(d2f_mask)) then
                    if(.not.d2f_mask(d,d2)) cycle
                end if      
                !because the Hessian is symmetric, we don't need to recalculate values
                !that were already calculated
                if(d2 > d) then
                    if(present(d2f_mask)) then
                        if(d2f_mask(d2,d)) then !we can only do this if the value was calculated before
                            d2f(d,d2) = d2f(d2,d)
                            cycle
                        end if
                    else
                        d2f(d,d2) = d2f(d2,d)
                        cycle
                    end if
                end if
                                
                d2f(d,d2) = 0d0
                counter = 0
                dimindx = 1
                do while(dimindx(1) <= size(Q%grid(1)%x))
                    counter = counter + 1
                    prod = 1d0
                    do i = 1,Q%Ndim
                        if (i == d .and. i == d2) then !second derivative needed
                            prod = prod*Q%k1d(i)%d2k(x(i),Q%grid(i)%x(dimindx(i)),Q%k1d(i)%par)
                        else if(i == d .or. i == d2) then !first derivative needed
                            prod = prod*Q%k1d(i)%dk(x(i),Q%grid(i)%x(dimindx(i)),Q%k1d(i)%par)
                        else
                            prod = prod*Q%k1d(i)%k(x(i),Q%grid(i)%x(dimindx(i)),Q%k1d(i)%par)
                        end if
                    end do
                    d2f(d,d2) = d2f(d,d2) + Q%alpha(counter)*prod
                            
                    !increase dimindx
                    dimindx(Q%Ndim) = dimindx(Q%Ndim) + 1
                    do i = Q%Ndim,2,-1
                        if(dimindx(i) > size(Q%grid(i)%x)) then
                            dimindx(i)   = 1
                            dimindx(i-1) = dimindx(i-1) + 1
                        end if
                    end do
                end do 
            end do   
        end do
    end if
    return 
end subroutine evaluate_kernel_slow

!-----------------------------------------------------------------------
!> @author
!> Oliver T. Unke, University of Basel
!
!> @brief
!> Calculates the error bound of the kernel at position x
!
!> @details
!> The calculation of the error bound is rather expensive, but it can be
!! useful in order to check the quality of the interpolation (or where
!! it would be useful to introduce more points). This is only possible
!! if coefficients were calculated using the slow method, as the full
!! inverse of the kernel matrix is needed for this (the reason the fast
!! method to calculate coefficients is fast is precisely because the
!! inverse matrix is never computed!)
!-----------------------------------------------------------------------
subroutine calculate_error_bound(Q,x,errorbound)
    implicit none
    class(kernel)                             :: Q
    real(kind(0d0)), dimension(:), intent(in) :: x !point at which to evaluate
    real(kind(0d0)), intent(out)              :: errorbound
    integer, dimension(size(x))               :: dimindx
    real(kind(0d0))                           :: prod
    real(kind(0d0))                           :: Qxx, cardinal, f !kernel function at point xx
    real(kind(0d0)), dimension(Q%Npoint)      :: Qvec !vector that contains all xix
    integer                                   :: i, j, counter
    
    
    if(allocated(Q%invQ)) then 
        !loop over all possible grid points to construct Qvec
        counter = 0
        dimindx = 1 !set dimindx = 1
        do while(dimindx(1) <= size(Q%grid(1)%x))
            counter = counter + 1
            if(Q%valueIsPresent(counter)) then
                Qvec(counter) = 1d0
                do i = 1,Q%Ndim
                    Qvec(counter) = Qvec(counter)*Q%k1d(i)%k(x(i),Q%grid(i)%x(dimindx(i)),Q%k1d(i)%par)
                end do   
            else
                Qvec(counter) = 0d0 
            end if
                    
            !increase dimindx
            dimindx(Q%Ndim) = dimindx(Q%Ndim) + 1
            do i = Q%Ndim,2,-1
                if(dimindx(i) > size(Q%grid(i)%x)) then
                    dimindx(i)   = 1
                    dimindx(i-1) = dimindx(i-1) + 1
                end if
            end do
        end do
        
        !calculate Qxx
        Qxx = 1d0
        do i = 1,Q%Ndim
            Qxx = Qxx * Q%k1d(i)%k(x(i),x(i),Q%k1d(i)%par)
        end do
        
        !calculate f
        f = sqrt(sum(pack(Q%alpha,Q%valueIsPresent)*pack(Q%values,Q%valueIsPresent)))
        
        errorbound = 0d0
        do i = 1,Q%Ndim
            if(.not.Q%valueIsPresent(i)) cycle
            !calculate the regularized cardinal function
            !cardinal = sum(Q%invQ(i,:)*Qvec)
            cardinal = sum(pack(Q%invQ(i,:),Q%valueIsPresent)*pack(Qvec,Q%valueIsPresent))
            !calculate the sum of Qvec times cardinal function
            errorbound = errorbound + sum(cardinal*Qvec)
        end do
        errorbound = abs(Qxx - errorbound)*f    
        return  
    else
        Qxx = -huge(0d0)
        errorbound = sqrt(Qxx) !return NaN on purpose
        return
    end if
    return 
end subroutine calculate_error_bound

end module RKHS

!!this is the naive implementation of the summation, which scales O(Npoint^2) => redundant
!subroutine calculate_single_sum_slow(Q,sigma,m,k)
!    implicit none
!    type(kernel) :: Q
!    real(kind(0d0)), dimension(Q%Npoint), intent(out) :: sigma
!    integer, dimension(Q%Ndim), intent(in) :: m,k
!    type(f_ptr), dimension(Q%Ndim)     :: fmk !array of function pointers
!    real(kind(0d0)), dimension(Q%Ndim) :: p2k !array of p values
!    real(kind(0d0)) :: prod  !intermediate product
!    integer, dimension(Q%Ndim) :: z !index vector
!    integer, dimension(Q%Ndim) :: loop_start, loop_end, loop_idx
!    integer :: d,idx
!       
!    !build the array of function pointers and p values for this combination of m and k
!    !and also initialize loop_start, loop_end and loop_step for each dimension (depends on m)
!    !NOTE: m = 0 encodes f2 functions, m=1 encodes f3 functions
!    do d = 1,Q%Ndim
!        p2k(d) = Q%k1d(d)%p2(k(d)) !p2k is independent of m
!        if(m(d) == 0) then !this is like m=2 in the paper (f2k functions)
!            fmk(d)%f => Q%k1d(d)%f2(k(d))%f !assign to f2k function
!            loop_start(d) = 0 !loop start is always 1 for f2k functions
!            !loop_end(d) needs to be reassigned everytime for f2k functions        
!        else               !this is like m=3 in the paper (f3k functions)
!            fmk(d)%f => Q%k1d(d)%f3(k(d))%f !assign to f3k function
!            !loop_start(d) needs to be reassigned every time for f3k functions
!            loop_end(d) = size(Q%grid(d)%x) !loop end is always N for f3k functions                         
!        end if    
!    end do  
!    
!    z   = 0   !initialize all indices to 0 (below first point in each grid)
!    idx = 0   !global index (for the right entry of sigma)
!    
!    !==>outer loop over all possible grid points
!    do while(z(1) <= size(Q%grid(1)%x)) 
!        idx        = idx + 1 !increment index
!        sigma(idx) = 0d0 !initialize sum to 0
!        
!        !assign appropriate loop_start/loop_end, depending on z
!        do d = 1,Q%Ndim
!            if(m(d) == 0) then !this is like m=2 in the paper (f2k functions)
!                loop_end(d) = z(d)
!            else !this is like m=3 in the paper (f3k functions)
!                loop_start(d) = z(d)+1
!            end if  
!            loop_idx(d) = loop_start(d) !initialize loop_idx   
!        end do
!        
!        !==> inner loop over all necessary grid points
!        do while(loop_idx(1) <= loop_end(1))
!            !only do this if the loop_idx is meaningful
!            if(.not.any(loop_idx > loop_end) .and. &
!               .not.any(loop_idx < 1)) then
!                prod = 1d0 !initialize product
!                do d = 1,Q%Ndim
!                    prod = prod * p2k(d)*fmk(d)%f(Q%grid(d)%x(loop_idx(d)),Q%k1d(d)%par)
!                end do
!                sigma(idx) = sigma(idx) + Q%alpha(Q%get_alpha_idx(loop_idx)) * prod !add to sum
!            end if
!                            
!            !increase loop_idx successively (loops over all necessary grid points this way)
!            loop_idx(Q%Ndim) = loop_idx(Q%Ndim) + 1
!            do d = Q%Ndim,2,-1
!                if(loop_idx(d) > loop_end(d)) then
!                    loop_idx(d)   = loop_start(d)
!                    loop_idx(d-1) = loop_idx(d-1) + 1
!                end if
!            end do
!        end do
!        
!        !increase z successively (loops over all possible grid points this way)
!        z(Q%Ndim) = z(Q%Ndim) + 1
!        do d = Q%Ndim,2,-1
!            if(z(d) > size(Q%grid(d)%x)) then
!                z(d)   = 0
!                z(d-1) = z(d-1) + 1
!            end if
!        end do
!    end do
!    return
!end subroutine calculate_single_sum_slow