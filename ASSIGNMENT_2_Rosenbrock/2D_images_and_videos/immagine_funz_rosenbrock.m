clc
clear

%This section prints the graph of Rosenbrock function for 2 dimensions

figure("Color", "w")

%f = @(x) 100*(x(2)-x(1)^2)^2+(1-x(1))^2; %2d rosenbrock
f_16 = @(x) (1-cos(x(1)) - sin(x(2))) + 2* ((1-cos(x(2))) + sin(x(1))) ; %2d problem 16
f_16_ = @(x,y) 1-cos(x) - sin(y) + 2*(1-cos(y) + sin(x)); 

%X = [1.2,1.2; -1.2,1; 0,-8];
%[x_opt, fx_opt, k, xseq, singular] = Nelder_Mead(X, f);
x_opt_16 = [-atan(2), atan(1/2)];  %should be a minimum
fx_opt_16 = f_16_(x_opt_16(1), x_opt_16(2));


%per il disegno
g = @(x,y) 100*(y-x.^2).^2+(1-x).^2;
[Xdis,Ydis]= meshgrid(linspace(-5,5,1000),linspace(-2,2,1000));
Zdis= f_16_(Xdis,Ydis);


surf(Xdis,Ydis,Zdis)
shading interp;         % Applica un'ombreggiatura continua
colormap(parula);          % Applica il colormap "jet"
xlabel('x_1');
ylabel('x_2');
zlabel('f(x)');

hold on; % Permette di aggiungere altri elementi sullo stesso grafico

plot3(x_opt_16(1), x_opt_16(2), fx_opt_16, 'ro', 'MarkerSize', 15, 'MarkerFaceColor', 'r');

% Creazione dell'etichetta con x_opt
label = ['x ','opt = ','(' num2str(x_opt_16(1), '%.2f') ', ' num2str(x_opt_16(2), '%.2f') ', ' num2str(fx_opt_16, '%.2f') ')'];

% Aggiunta della label al punto
text(x_opt_16(1), x_opt_16(2), fx_opt_16+4.4, label,'Color', 'black','HorizontalAlignment',...
    'center','VerticalAlignment', 'bottom','BackgroundColor', 'white','EdgeColor', 'black', 'clipping', 'off' ); 