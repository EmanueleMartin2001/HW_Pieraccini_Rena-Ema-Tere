clc
clear


% This section prints the Rosenbrock function and the approximation
% parboloid w.r.t. a generic starting point, (which is the function that has to be
% minimized at the first step)


x0=500;
y0=100;

f_Rosenbrock = @(x) 100*(x(2)-x(1)^2)^2 + (1-x(1))^2;

gradf_Rosenbrock = @(x) [400*x(1)^3 - 400*x(1)*x(2) + 2*x(1) - 2 ; 200*( x(2) - x(1)^2 )];

Hessianf_Rosenbrock = @(x) [ 1200*x(1)^2 - 400*x(2) + 2 , -400*x(1); -400*x(1) , 200];


f= @(x,y)100*(y-x.^2).^2 + (1-x).^2;
gradf= @(x,y) [-400*x.*(y-x.^2)-2*(1-x); 200*(y-x.^2)];
hessf= @(x,y) [-400*(y-x.^2)+800*x.^2 , -400*x ; -400*x , 200];


[Xfunz,Yfunz]= meshgrid(linspace(-1000,1000,1000),linspace(-1000,1000,1000));
[Xdis,Ydis]= meshgrid(linspace(-200,1000,500),linspace(-200,200,500));
Zdis= f(Xfunz,Yfunz);
 


surf(Xfunz,Yfunz,Zdis)
shading interp;         % Applica un'ombreggiatura continua
colormap(hsv);          % Applica il colormap "parula"
xlabel('x_1');
ylabel('x_2');
zlabel('f(x)');
hold on;

plot3(x0, y0, f(x0,y0), 'ro', 'MarkerSize', 15, 'MarkerFaceColor', 'b');

hold on;

fxk=f(x0,y0);
grad=gradf(x0,y0);
a=grad(1);
b=grad(2);
px=Xdis-x0;
py=Ydis-y0;
apx=a*px;
bpy=b*py;
gradp=apx+bpy;
pianotang=fxk+gradp;

% surf(Xdis,Ydis,pianotang)
% shading interp;         % Applica un'ombreggiatura continua
% colormap(jet); 

hold on;



hess=hessf(x0,y0);
e=hess(1,1);
g=hess(1,2);
c=hess(2,1);
d=hess(2,2);

hessp= px.*px*e+px.*py*c+px.*px*g+py.*py*d;

parabolatang=pianotang+hessp;

surf(Xdis,Ydis,parabolatang)
% shading interp;         % Applica un'ombreggiatura continua
% colormap(gray);

