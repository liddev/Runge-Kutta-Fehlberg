#######  UNIVERSIDAD NACIONAL MAYOR DE SAN MARCOS   #######
#######     METODO DE RUNGE KUTTA DE 2DO ORDEN      #######
#######        PROBLEMA DEL VALOR INICIAL           #######




function [A]=RungeKutta_Orden2;
  fprintf('\n');
  fprintf('METODO DE RUNGE KUTTA DE 2DO ORDEN - METODO DE HEUN \n');
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
  k2=h*(f(x+((2/3)*h),y+((2/3)*k1)));

  y=y+((1/4)*(k1+(3*k2)));
  x=a+(i*h);
  
  B(1,i+1)=x;
  B(2,i+1)=y;

endfor
A(:,1)=B(1,:);
A(:,2)=B(2,:);
plotmatrix (A)
endfunction