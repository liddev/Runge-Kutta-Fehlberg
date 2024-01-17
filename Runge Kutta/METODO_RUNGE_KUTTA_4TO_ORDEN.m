
#######  UNIVERSIDAD NACIONAL MAYOR DE SAN MARCOS   #######
#######     METODO DE RUNGE KUTTA DE 4TO ORDEN      #######
#######        PROBLEMA DEL VALOR INICIAL           #######
## (x*(e^(3*x)))-(2*y) / y(0)=0 / h=0.1 / 0<x<1


function [A]=RungeKutta_Orden4;
  fprintf('\n');
  fprintf('METODO DE RUNGE KUTTA DE 4TO ORDEN \n');
  funcion=input('Ingrese la funcion despejada, F(x,y): ','s');
  x=input('Ingrese el valor inicial x: ');
  y=input('Ingrese el valor inicial F(x) = y: ');

  a=input('Ingrese el inicio del intervalo: ');
  b=input('Ingrese el final del intervalo: ');
  h=input('Ingrese el tamaño de paso: ');
  f=inline(funcion);

 
#Creando matriz para almacenar datos:
  B(1,1)=x;
  B(2,1)=y;
  
#INICIO DEL BUCLE:
for i=1:((b-a)/h);
  k1=h*(f(x,y));
  k2=h*(f(x+(h/2),y+(k1/2)));
  k3=h*(f(x+(h/2),y+(k2/2)));
  k4=h*(f(x+h,y+k3));
 
  y=y+((k1+(2*k2)+(2*k3)+k4)/6);
  x=a+(i*h);
  
  B(1,i+1)=x;
  B(2,i+1)=y;

endfor
A(:,1)=B(1,:);
A(:,2)=B(2,:);
plotmatrix (A)
endfunction