
program bis1
  
  implicit none
  Integer,parameter::k=selected_real_kind(30,3000)
  real(k)::tol,x,err,a,b
  integer::bandera
  write(*,*) "introduce tu intervalo (de menor a mayor, primero uno y &
       &despues el otro)"
  read*, a
  read*, b
  write(*,*) "introduce tu tolerancia"
  read*, tol
  bandera=0
  call Biseccion(tol,a,b,x,err,bandera)
  if(bandera==0)then
     write(*,*) "no hay raíces en el intervalo o hay una cantidad par de &
          &raices"

  else
     write(*,*) "la aproximacion es:", x
     write(*,*) "el error aproximado porcentual es:",err
  
endif

contains
  subroutine Biseccion(tol,a,b,x,err,bandera)
    
    implicit none
    Integer,parameter::k=selected_real_kind(30,3000)
    real(k),intent(in)::tol
    real(k),intent(inout)::a
    real(k),intent(inout)::b
    real(k),intent(out)::x
    real(k),intent(out)::err
    integer,intent(out)::bandera
    real(k)::y
    integer::n,i
    y=0
    n=ceiling((log((b-a)/tol))/log(2.))    
    do i=1,n
       if(f(a)*f(b)<0)then          !1          
          x=(a+b)/2.
          err=abs((x-y)/x)*100
          y=x
          bandera=1       
          if(f(a)*f(x)<0)then       !2
             b=x
          else
             a=x          
          endif
       else
          exit      
       endif
    enddo
  endsubroutine Biseccion

  function f(u)
  
    implicit none
  
    integer,parameter::k=selected_real_kind(30,3000)
    real(k)::f,u
   
    f=u**2-13.247
   
  endfunction f
endprogram bis1
