program test_derivative
implicit none
integer :: np
integer, parameter :: fp=kind(1.0d0)
real(kind=fp), allocatable :: x(:), y(:)
real(kind=fp) :: der_left2, der_right2
real(kind=fp) :: der_left3, der_right3
real(kind=fp) :: der_left4, der_right4
real(kind=fp) :: err2_left, err3_left, err4_left
real(kind=fp) :: err2_right, err3_right, err4_right
integer :: i

do np = 4, 40
allocate(x(np), y(np))

do i =1, np
 x(i) = 1._fp/3._fp + real(i-1, fp)/real(np-1, fp)
 y(i) = my_function( x(i) )
enddo

select case(np)
    case(:2)
        write(*,*) 'Error: np < 2 , np= ' , np

    case default

    der_left2  = (y(1)-y(2)) / (x(1)-x(2))
    der_right2 = (y(np-1)-y(np)) / (x(np-1)-x(np))

    der_left3=der_left2+ (x(3)*(y(1)-y(2))+x(1)*(y(2)-y(3))+x(2)*(y(3)-y(1)))/ ( (x(1)-x(3))*(x(3)-x(2) ))

    der_right3=der_right2+ (x(np)*(y(np-2)-y(np-1))+x(np-2)*(y(np-1)-y(np))+x(np-1)*(y(np)-y(np-2))) / &
                                                            ( (x(np-2)-x(np-1))*(x(np-2)-x(np) ))

    der_left4 = (y(1) - y(2))/(x(1) - x(2)) + (x(1) - x(2))*(((-y(1) + y(2))/(x(1) - x(2)) + &
                (y(2) - y(3))/(x(2) - x(3)))/(-x(1) + x(3)) + ((x(1) - x(3))*(-(((-y(1) + y(2))/(x(1) - x(2)) &
                + (y(2) - y(3))/(x(2) - x(3)))/(-x(1) + x(3))) + ((-y(2) + y(3))/(x(2) - x(3)) + &
                (y(3) - y(4))/(x(3) - x(4)))/(-x(2) + x(4))))/(-x(1) + x(4)))

    der_right4 = (y(np) - y(np-1))/(x(np) - x(np-1)) + (x(np) - x(np-1))*(((-y(np) + y(np-1))/(x(np) - x(np-1)) + &
                (y(np-1) - y(np-2))/(x(np-1) - x(np-2)))/(-x(np) + x(np-2)) + &
                ((x(np) - x(np-2))*(-(((-y(np) + y(np-1))/(x(np) - x(np-1)) &
                + (y(np-1) - y(np-2))/(x(np-1) - x(np-2)))/(-x(np) + x(np-2))) + ((-y(np-1) + y(np-2))/(x(np-1) - x(np-2)) + &
                (y(np-2) - y(np-3))/(x(np-2) - x(np-3)))/(-x(np-1) + x(np-3))))/(-x(np) + x(np-3)))

end select

    err2_left = log( abs(der_left2-my_function_der1(x(1))) ) /log(10._fp)
    err2_right = log( abs(der_right2-my_function_der1(x(np))) ) /log(10._fp)

    err3_left = log( abs(der_left3-my_function_der1(x(1))) ) /log(10._fp)
    err3_right = log( abs(der_right3-my_function_der1(x(np))) ) /log(10._fp)

    err4_left = log( abs(der_left4-my_function_der1(x(1))) ) /log(10._fp)
    err4_right = log( abs(der_right4-my_function_der1(x(np))) ) /log(10._fp)

    write(*,'(i10, 10F20.14)') np, err2_left, err2_right, err3_left, err3_right, err4_left, err4_right
deallocate(x,y)

enddo

contains

! for 6 points of less, in [1/3, 1] , the 4th order derivative is worse than the others
! for  1._fp / x**2 and the right-end derivative
pure function my_function(x)
    real(kind=fp), intent(in) :: x
    real(kind=fp) :: my_function
    my_function = 1._fp / x**2
    return
end function my_function

pure function my_function_der1(x)
    real(kind=fp), intent(in) :: x
    real(kind=fp) :: my_function_der1
    my_function_der1 = -2._fp / x**3
    return
end function my_function_der1


end program test_derivative




! given an array x(i) and y(i), gives the derivative at the leftmost and rightmost point
! computed by interpolations
