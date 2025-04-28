clc
clear

%This section prints the graph of Rosenbrock function for 2 dimensions


f = @(x) 100*(x(2)-x(1)^2)^2+(1-x(1))^2;
X = [1.2,1.2; -1.2,1; 0,-8];


[x_opt, fx_opt, k, xseq, singular] = Nelder_Mead(X, f);


%per il disegno
g = @(x,y) 100*(y-x.^2).^2+(1-x).^2;
[Xdis,Ydis]= meshgrid(linspace(-100,100,1000),linspace(-100,100,1000));
Zdis= g(Xdis,Ydis);


surf(Xdis,Ydis,Zdis)
shading interp;         % Applica un'ombreggiatura continua
colormap(parula);          % Applica il colormap "jet"
xlabel('x_1');
ylabel('x_2');
zlabel('f(x)');

hold on; % Permette di aggiungere altri elementi sullo stesso grafico

plot3(x_opt(1), x_opt(2), fx_opt, 'ro', 'MarkerSize', 15, 'MarkerFaceColor', 'r');

% Creazione dell'etichetta con x_opt
label = ['x ','opt = ','(' num2str(x_opt(1), '%.2f') ', ' num2str(x_opt(2), '%.2f') ', ' num2str(fx_opt, '%.2f') ')'];

% Aggiunta della label al punto
text(x_opt(1), x_opt(2), fx_opt+200000, label,'Color', 'black','HorizontalAlignment',...
    'center','VerticalAlignment', 'bottom','BackgroundColor', 'white','EdgeColor', 'black' ); 
     