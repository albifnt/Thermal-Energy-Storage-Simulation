%n_cells initial
n_cells=10;
grid_refinement=8;
r=2; % Coefficient that multiplies the number of cells

% Load
B=load("Manufactured_Solid.txt");
B=B';
Norms=load('Solid_Norm.txt');
Refinement=load('Refinement.txt');

Norm_1=Norms(:,1);
Norm_inf=Norms(:,2);

% Plot
figure(grid_refinement+1)
plot(Refinement, Norm_1,'b-o',Refinement, Norm_inf,'r-o');
axis equal;
title('Solid GRS');
xlabel('log(h)');
ylabel('log(E)');
legend('E_1', 'E_i_n_f');


% Ovs Plot
for i=1:grid_refinement-1
    p_1(i)=(Norm_1(i+1)-Norm_1(i))/(Refinement(i+1)-Refinement(i))
end
for i=1:grid_refinement-1
    p_inf(i)=(Norm_inf(i+1)-Norm_inf(i))/(Refinement(i+1)-Refinement(i))
end

figure(grid_refinement+2)
plot(Refinement(2:grid_refinement),p_1,'b-o', Refinement(2:grid_refinement),p_inf, 'r-o');
title('Solid GRS');
xlabel('log(h)');
ylabel('log(\DeltaE)/log(\Deltah)');
legend('E_1', 'E_i_n_f');
ylim([0.5 2.5]);


% Initialization
start=1;
endl=0;
gap=0;

% For
for i=0:(grid_refinement-1)
    
    % Number of cells
    gap=n_cells*(r^i);
    endl=endl+gap;
    % Re-initialization
    start=start + gap;   
end

start= start-gap;

% Plot Exact solution
Vector_Ex=B(2,start:endl);
Vector_Ex_Ax=B(3,start:endl);
%plot(Vector_Ex_Ax,Vector_Ex);

% Initialization
start=1;
endl=0;
gap=0;
 
% For
for i=0:(grid_refinement-1)
    
    % Number of cells
    gap=n_cells*2^i;
    endl=endl+gap;
    % Vector to plot
    Vector_Tf=B(1,start:endl);
    Vector_Ax=B(3,start:endl);
    
    % Plot
    figure(i+1)
    plot(Vector_Ax, Vector_Tf);
    hold on;
    plot(Vector_Ex_Ax, Vector_Ex);
    title(sprintf('N_c_e_l_l_s = %d', n_cells*(r^i)));
    %%title('N=20480');
    xlabel('x axis');
    legend('Manufactured','Exact solution');
    ylim([-1 1]);
    hold off;
    
    % Re-initialization
    start=start + gap;  
end

