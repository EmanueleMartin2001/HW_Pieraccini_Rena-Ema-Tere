%function that generates a video of the evolution of the symplex in the 2D 
%space w.r.t the iterations. This is helpful to show how this method
%searches the optimum value


%NB!
% The function used in this section, are slight modifications of the
% original Nelder Mead function, in order to show graphs and plots, to
% represent visually how the solution is searched by the algorithm.

% These changes are not already implemented in the original function
% because, in order to make the plots, it's necessary to save all the
% values found during the process, wich makes the resolution much slower
% and is opposite of the aim of our code.


function [x_opt, fx_opt, k, x_seq, singular] = Funz_video_nelder(X, f, varargin) 


%% input parsing
p = inputParser;

% adding optional parameters with default values

addOptional(p, 'kmax', 200);  
addOptional(p, 'tol', 10^(-6));
addOptional(p, 'rho', 1);
addOptional(p, 'chi', 2);
addOptional(p, 'gamma', 0.5);
addOptional(p, 'sigma', 0.5);

parse(p, varargin{:}); % add those variables as input
kmax= p.Results.kmax;
tol = p.Results.tol;
rho = p.Results.rho;
chi= p.Results.chi;
gamma = p.Results.gamma;
sigma = p.Results.sigma;


% inizialization for default output values
x_opt = [];
fx_opt = [];
k = 0;
x_seq = [];
singular = false;

%% 0 singularity check
%at first we have to check that the simplex is not singular
n = size(X,2);
simplex_vectors = ones(n);
for i = 2:n+1
    simplex_vectors(i-1,:) = X(i,:)-X(1,:);
end
if (det(simplex_vectors) == 0) %checks the non singularity
    singular = true;
    return
end


%% 1 initialization
%first we have to reorder: we can save the values of the function evalutations
f_evaluation = ones(n+1,1);
for i = 1:n+1
    f_evaluation(i) = f(X(i,:)); %evaluates the function in the point x_i
end

[f_evaluation,new_indices] = sort(f_evaluation); %in ascending order
%new_indices(k) is the index of the k-th best point of the X matrix
%so new_indices(n+1) is the worst point (x_{n+1} in the slides) 
%f_evaluation(k) is f(x_{k})
values_indices(:,1) = f_evaluation;
values_indices(:,2) = new_indices; %this matrix contains 2 columns;
%a row of this matrix contains f(x) and his corresponding index in the
%matrix X

%-------------------------------------------------------------------------------
%-------------------------------------------------------------------------------
%-------------------------------------------------------------------------------
%inizializzo video
video = VideoWriter('symplex_evolution.avi');  % Nome del file video
open(video);
figure;
hold on;
xlabel('x');
ylabel('y');
grid on;
[xmin, xmax, ymin, ymax] = deal(-15, 15, -15, 15);
axis([xmin, xmax, ymin, ymax]);
active_symplex=[]; %voglio avere a video solo 4 simplessi alla volta, quando questo vettore ha 
                  %lunghezza >=4, elimino il più vecchio
colore_attivo=[];

%lista colori
colors = {'red', 'blue', 'green', 'cyan', 'magenta'};
n_tot_colori = length(colors);
%-------------------------------------------------------------------------------
%-------------------------------------------------------------------------------
%-------------------------------------------------------------------------------

%% the algorithm
while k < kmax &&... % check ix x_{1} and x_{n+1} are close or f(x_{1}) is close to f(x_{n+1})
        (abs((values_indices(1,1)-values_indices(n+1,1))/values_indices(1,1)) > tol...
        || norm(X(values_indices(n+1,2),:)-X(values_indices(1,2),:),1)/norm(X(values_indices(1,2),:),1) > tol)
    x_baricenter = mean(X([1:values_indices(n+1,2)-1, values_indices(n+1,2)+1:end],:),1);
    %centroid of the simplex (not taking in account x_{k+1}, the worst point)

    xr = x_baricenter + rho*(x_baricenter - X(values_indices(n+1,2),:));
    f_xr = f(xr); %point reflected evaluated
    k = k + 1;

    %----------------------------------------------------------------------------------------
    %----------------------------------------------------------------------------------------
    %----------------------------------------------------------------------------------------
    %plot del simplesso per le prime iterazioni, che saranno usate come
    %frame
    if k<=100
        X_finale1=X;
        newrow1=X_finale1(1,:);
        X_finale1=[X_finale1;newrow1];%aggiungo il primo vertice anche in fondo alla matrice, così 
                                      %da poter chiudere il poligono con tutti i
                                      %lati
        x_point1 = X_finale1(:,1); % Coordinata x del punto
        y_point1 = X_finale1(:,2); % Coordinata y del punto
        
        %plot
        %scelgo il colore della figura corrente
        indice=mod(k,n_tot_colori)+1;
        colore_attuale=colors{indice};
        
        %sto creando l'oggetto grafico symplex(con plot) che viene aggiunto
        %a figure
        symplex=plot(x_point1, y_point1,'-o', 'LineWidth', 1.5,'MarkerSize', 5, 'MarkerFaceColor',colore_attuale );
        colore=fill(x_point1, y_point1, colore_attuale, 'FaceAlpha', 0.3); % Riempie l'area con un colore
        title(sprintf('Iterazione %d', k));
        
        if mod(k,10)==0
            xmin = min(x_point1); % Calcola il valore minimo di x
            xmax = max(x_point1); % Calcola il valore massimo di x
            ymin = min(y_point1); % Calcola il valore minimo di y
            ymax = max(y_point1); % Calcola il valore massimo di y
            axis([xmin-20/k, xmax+20/k, ymin-20/k, ymax+20/k]);
        end
        




        %aggiorno i simplessi attivi
        active_symplex=[active_symplex,symplex]; %aggiungo il simplesso agli attivi
        colore_attivo=[colore_attivo,colore];
        if length(active_symplex)>4
            delete(active_symplex(1)); %sto eliminando da figure il plot dei 4 più vecchio
            active_symplex(1)=[]; %tolgo l'elemento anche dal vettore dei simplessi attivi
            delete(colore_attivo(1)); 
            colore_attivo(1)=[];
        end
        
        % Acquisizione del frame
        frame = getframe(gcf);
        writeVideo(video, frame);
        if k<=10
            pause(0.6);  % Pausa per visualizzare l'evoluzione
        elseif k>10 & k<25
            pause(0.4)
        elseif k>=25 & k<=50
            pause(0.2)
        else
            pause(0.08)
        end

    end
    %---------------------------------------------------------------------------------------------
    %----------------------------------------------------------------------------------------
    %----------------------------------------------------------------------------------------



    if values_indices(1,1) <= f_xr && f_xr < values_indices(n,1) %%% REFLECTION %%%
        X(values_indices(n+1,2),:) = xr; %x{n+1} --> xr
        values_indices(n+1,1) = f_xr;    %f_x{n+1} --> f_{xr}  
        values_indices = sortrows(values_indices,1); %points reordered
        x_seq(k,:) = X(values_indices(1,2),:); % best point found so far
        % we know that x_seq(k,:) = x_seq(k-1,:) due to the if above, but
        % we still want to save this point in order to have an
        % understanding about the convergence speed of the method
        continue
    elseif f_xr < values_indices(1,1) %%% EXPANTION %%%
        xe = x_baricenter + chi*(xr-x_baricenter);
        f_xe = f(xe); %point expanded evaluated
        if f_xe < f_xr
            X(values_indices(n+1,2),:) = xe; %x_{n+1} --> xe
            values_indices(n+1,1) = f_xe;    %f_x{n+1} --> f_{xe}  
            values_indices = sortrows(values_indices,1); %points reordered
            x_seq(k,:) = X(values_indices(1,2),:); % best point found so far
            continue
        else
            X(values_indices(n+1,2),:) = xr; %x_{n+1} --> xr
            values_indices(n+1,1) = f_xr;    %f_x{n+1} --> f_{xr}  
            values_indices = sortrows(values_indices,1); %points reordered
            x_seq(k,:) = X(values_indices(1,2),:); % best point found so far
            continue
        end
    else %%% CONTRACTION %%%  (f_xr >= values_indices(n,1))
        if values_indices(n+1,1) < f_xr %x_{n+1} is better then xr
            xc = x_baricenter - gamma*(x_baricenter - X(values_indices(n+1,2),:));
        else %xr is better then x_{n+1}
            xc = x_baricenter - gamma*(x_baricenter-xr);
        end
        f_xc = f(xc);
        if f_xc < values_indices(n+1,1)
            X(values_indices(n+1,2),:) = xc; %x_{n+1} --> xc
            values_indices(n+1,1) = f_xc;    %f_x{n+1} --> f_{xc}  
            values_indices = sortrows(values_indices,1); %points reordered
            x_seq(k,:) = X(values_indices(1,2),:); % best point found so far
            continue
        else %%% SHRINKAGE %%%
            for i = 2:n+1
                X(values_indices(i,2),:) = X(values_indices(1,2),:) ...
                + sigma*(X(values_indices(i,2),:) - X(values_indices(1,2),:));
                values_indices(i,1) = f(X(values_indices(i,2),:));
            end
            values_indices = sortrows(values_indices,1);
            x_seq(k,:) = X(values_indices(1,2),:); % best point found so far
    
    
        end
    
   
    end

    

%chiusura while loop
end
x_opt = X(values_indices(1,2),:); % best point
fx_opt = f(x_opt);

%---------------------------------------------------------------------------------------------
%----------------------------------------------------------------------------------------
%----------------------------------------------------------------------------------------

% Chiusura del file video
close(video);
disp('Video generato con successo.');

%---------------------------------------------------------------------------------------------
%----------------------------------------------------------------------------------------
%----------------------------------------------------------------------------------------




%chiusura funzione
end