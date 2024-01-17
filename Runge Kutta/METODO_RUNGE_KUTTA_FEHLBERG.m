function [B]=RungeKuttaFehlberg;
  fprintf('\n');
  fprintf('METODO DE CONTROL DE ERROR - RUNGE KUTTA FEHLBERG \n');
  funcion=input('Ingrese la funcion despejada, F(x,y): ','s');
  x=input('Ingrese el valor inicial x: ');
  y=input('Ingrese el valor inicial F(x) = y: ');
  f=inline(funcion);
  funcionsolucion=input('Ingrese la funcion solucion: ', 's');
  g=inline(funcionsolucion);
  
  a=input('Ingrese el inicio del intervalo: ');
  b=input('Ingrese el final del intervalo: ');
  hmax=input('Ingrese el tamaño de paso maximo: ');
  hmin=input('Ingrese el tamaño de paso minimo: ');
  tol=input('Ingrese la tolerancia de error: ');
  
  w=y;
  h=hmax;
  band=1;
  
  fprintf('\n');
  fprintf('     xi');
  fprintf('         yi=y(xi)');
  fprintf('        wi');
  fprintf('        hi');
  fprintf('         |yi-wi|');
  fprintf('        wii');
  fprintf('       |yi-wii|');
  fprintf('\n');
  fprintf('%12.3f',x);
  fprintf('%12.3e',g(x));
  fprintf('%12.3f',w);
  fprintf('         ');
  fprintf('%12.3e',g(x)-w);
  fprintf('\n');
  
  while band==1
    k1=h*f(x,w);
    k2=h*f(x+(1/4)*h,w+(1/4)*k1);
    k3=h*f(x+(3/8)*h,w+(3/32)*k1+(9/32)*k2);
    k4=h*f(x+(12/13)*h,w+(1932/2197)*k1-(7200/2197)*k2+(7296/2197)*k3);
    k5=h*f(x+h,w+(439/216)*k1-8*k2+(3680/513)*k3-(845/4104)*k4);
    k6=h*f(x+(1/2)*h,w-(8/27)*k1+2*k2-(3544/2565)*k3+(1859/4104)*k4-(11/40)*k5);
    R=(1/h)*abs((1/360)*k1-(128/4275)*k3-(2197/75240)*k4+(1/50)*k5+(2/55)*k6);
  if R<=tol
    x=x+h;
    wi=w+(16/135)*k1+(6656/12825)*k3+(28561/56430)*k4-(9/50)*k5+(2/55)*k6;
    w=w+(25/216)*k1+(1408/2565)*k3+(2197/4104)*k4-(1/5)*k5;
    
    fprintf('%12.7f',x);
    fprintf('%12.7f',g(x));
    fprintf('%12.7f',w);
    fprintf('%12.7f',h)
    fprintf(' %12.7e',abs(g(x)-w));
    fprintf('%12.7f',wi);
    fprintf(' %12.7e',abs(g(x)-wi));
    fprintf('\n');
  end
  delta=0.84*(tol/R)^(1/4);
  if delta<=0.1
    h=0.1*h;
  elseif delta>=4
    h=4*h;
  else h=delta*h;
  end
  if h>hmax
    h=hmax;
  end
  if x>=b
    band=0;
  elseif (x+h)>b
    h=b-x;
  elseif h<hmin
    band=0;
  end
  endwhile
endfunction