program trapoprince
  implicit none

  real :: a, b, integral
  integer :: m

  ! los limites, pos que lea de a-b
  print*, "Dame los limites de la integracion (a, b):"
  read*, a, b

  print*, "Cuantos trapecios quieres usar?"
  read*, m

  ! necesito llamar la subrutina porque pos lol no voy a rehacer otro programa de trapezios, voy a reusar el que ya esta
  call trapecio(m, a, b, integral)
  print*, "la integral deseada es", integral

end program trapoprince

! Subrutina 
subroutine trapecio(m, a0, b0, resultado)
  implicit none

  real, intent(in) :: a0, b0
  real, intent(inout) :: resultado
  integer, intent(in) :: m
  integer :: i
  real :: h, suma
  real :: f

  h = (b0 - a0) / m
  suma = 0.0

  do i = 1, m - 1
     suma = suma + f(a0 + i * h)
  end do

  resultado = (h / 2.0) * (f(a0) + f(b0) + 2 * suma)
end subroutine trapecio

! Función que tengo que integrar, ya estoy cansado jefe
real function f(x0)
  implicit none
  real, intent(in) :: x0

  f = 4*(3.1415)*171601074.4*(21.6**4)/((21.6+x0)**2)
end function f
