program DMH

  implicit none
  integer,parameter::k=selected_real_kind(30,3000)
  
  real(k),allocatable::r1(:),M(:),r2(:),Ma(:)
  real(k)::rh
  integer::i,N,lk

  N=100000
  lk = 10000
  allocate(r1(lk))
  allocate(m(lk))
  allocate(r2(N))
  allocate(Ma(N))
  
  rh=256
  do i=1,lk
     r1(i)=i*rh/real(lk)            !construccion de los radios

  enddo
  call rng (Ma,N)
  call trapecio(r1,M,lk)
  call Biseccion(Ma,r2,rh,N)
  call conv (ma,r2,N)
  call conv_2(m,r1,lk)

  open(unit=22,file="grafica_esfera.gnu",status="unknown")
  write(22,*) "set view equal xy"
  write(22,*) "plot 'lista.dat' every :::0::999 using 5:4 with dots title 'MvsR'"
  call system('gnuplot -p grafica_esfera.gnu')
  close(22)

  open(unit=22,file="grafica_esfera.gnu",status="unknown")
  write(22,*) "set view equal xy"
  write(22,*) "plot 'lista2.dat' every :::0::999 using 5:4 with lines title 'MvsR_analitico'"
  call system('gnuplot -p grafica_esfera.gnu')
  close(22)
  
  stop
contains  
  
!!!!!Trapecio  
  subroutine trapecio (r1,M,n)
    implicit none  
    real(k),intent(in)::r1(:)
    real(k),intent(out)::M(:)
    real(K)::a=0,p,sum
    integer::i,j,n

    
    do j=1,n
       sum = 0.
       p=(r1(j)-a)/real(n)
              
       do i=1,n-1
          sum=sum+ f(a+i*p)
       enddo

       M(j)=(p/2.)*(f(a)+f(r1(j))+2*sum)
       
    end do
  endsubroutine trapecio
  
  function f(r)
    implicit none
    real(k)::f,r          !valores de entrada/salida
    real(k)::pi=3.141592653589, ro=8562653.6767, a=21.6       !ctes.para no tener que escribirlas alm
    
    f=(4*pi*ro*(a**4))/((a+r)**2)
    
  endfunction f

!!!!!!!!!Biseccion!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 subroutine Biseccion(ma,r2,rh,N)
    
   implicit none
   real(k),intent(in)::ma(:)
   real(k),intent(out)::r2(:)
   real(k),intent(in)::rh
   integer,intent(in)::N
   real(k)::a,b
   integer::i,j,m

   do j=1,N
      a=0
      b=rh
      do i=1,200
         if((g(a,ma(j)))*(g(b,ma(j)))<0)then          !1          
            r2(j)=(a+b)/2.       
            if(g(a,ma(j))*g(r2(j),ma(j))<0)then       !2
               b=r2(j)
            else
               a=r2(j)  
            endif
         else
            exit      
         endif
      enddo
   end do
   
 endsubroutine Biseccion

 function g(t,h)
   implicit none
   real(k)::g,t,h
   real(k)::pi=3.141592653589, ro=8562653.6767, a=21.6   

   g=-((4*pi*ro*(a**4))/(a+t)) + 4*pi*ro*(a**3) - h
   
 end function g

!!!!RNG!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
 subroutine rng(ma,n)
   implicit none
   real(k),intent(out)::ma(:)
   integer,intent(in)::n
   integer :: i, j
   real, allocatable :: numeros(:)
   real(k) :: temp
   real(k) :: a=0, b=1e12
   

   ! Reservar memoria para el arreglo
   allocate(numeros(n))
   
   ! Inicializar el generador de números aleatorios
   call random_seed()
  
   ! Generar números aleatorios en el rango [a, b]
   do i = 1, n
      call random_number(temp)
      numeros(i) = b*temp
   end do
   
      
   ! Pasar resultados
   do i = 1, n
      ma(i)=numeros(i)
   end do
   ! Liberar memoria
   deallocate(numeros)
   
 end subroutine rng

!!!convertir las coordenadas!!!!!!!!!!!!!!!!!!!!
   subroutine conv(ma,r2,N)
     implicit none
     real(k),intent(in)::ma(:),r2(:)
     integer,intent(in)::N
     real(k)::x,y,z
     real(k)::theta,phi
     integer::i
     
     open(unit=10,file="lista.dat",status="replace")
     
     do i=1,N
        call random_number(theta)
        call random_number(phi)
        theta=theta*2.*acos(-1.)
        phi=acos(2.*phi-1)
        
        x=r2(i)*sin(phi)*cos(theta)
        y=r2(i)*sin(phi)*sin(theta)
        z=r2(i)*cos(phi)
        write(10,*) x,y,z, ma(i), r2(i)
     enddo
     close(10)
     open(unit=22,file="grafica_esfera.gnu",status="unknown")
     write(22,*) "set view equal xyz"
     write(22,*) "splot 'lista.dat' every :::0::999 using 1:2:3 with d&
          &ots title 'DMH'"
     call system('gnuplot -p grafica_esfera.gnu')
     close(22)
   endsubroutine conv
   
   subroutine conv_2(m,r1,N)
     implicit none
     real(k),intent(in)::m(:),r1(:)
     Integer,intent(in)::N
     real(k)::x,y,z
     real(k)::theta,phi
     integer::i
     
     open(unit=23,file="lista2.dat",status="replace")
     
     do i=1,N
        call random_number(theta)
        call random_number(phi)
        theta=theta*2.*acos(-1.)
        phi=acos(2.*phi-1.)
        
        x=r1(i)*sin(phi)*cos(theta)
        y=r1(i)*sin(phi)*sin(theta)
        z=r1(i)*cos(phi)
        write(23,*) x,y,z, m(i),r1(i)
     enddo
     close(23)
     open(unit=24,file="grafica_esfera_analitica.gnu",status="unknown")
     write(24,*) "set view equal xyz"
     write(24,*) "splot 'lista2.dat' every :::0::999 using 1:2:3 with d&
          &ots title 'DMH_analitico'"
     call system('gnuplot -p grafica_esfera_analitica.gnu')
     close(24)
   endsubroutine conv_2

     
   endprogram DMH
   
