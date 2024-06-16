% Remove
clear all

% Loading Diamater
D=load('Diameter.txt');
E=load('Efficiency.txt');
C=load('CapacityFactor.txt');
K=load('Temperature_difference.txt');

% Iteration
iteration=E(:,2);
start=1;
B=0;
gap=0;

% Plot Energy efficiency
figure(1)
for i=1:(length(iteration)-1)
    if iteration(i+1) - iteration(i) < 1
       gap=gap+iteration(i);
       plot(iteration(start:gap),E(start:gap,1),'linewidth',2);
       hold on;
       start=1+gap;
    end   
end

plot(iteration(start:end),E(start:end,1),'linewidth',2);
Q=title('Storage Design Study');
T=xlabel('Cycles')
Y=ylabel('Exergy Efficiency');
U=legend('Diameter = 4 m','Diameter = 5 m','Diameter = 6 m','Diameter = 7 m','Diameter = 8 m')
Q.FontSize=14;
T.FontSize=14;
Y.FontSize=14;
U.FontSize=14;
hold off;

% Iteration
start=1;
B=0;
gap=0;
t=1;
% Plot Capacity Factor
figure(2)
for i=1:(length(iteration)-1)
    if iteration(i+1) - iteration(i) < 1
       gap=gap+iteration(i);
       plot(iteration(start:gap),C(start:gap,1),'linewidth',2);
       hold on;
       start=1+gap;
    end   
end

plot(iteration(start:end),C(start:end,1),'linewidth',2);
Q=title('Storage Design Study');
T=xlabel('Cycles')
Y=ylabel('Capacity Factor');
U=legend('Diameter = 4 m','Diameter = 5 m','Diameter = 6 m','Diameter = 7 m','Diameter = 8 m')
Q.FontSize=14;
T.FontSize=14;
Y.FontSize=14;
U.FontSize=14;
hold off;

% Iteration
iteration=E(:,2);
start=1;
B=0;
gap=0;

% Plot Temperature Increase
figure(3)
for i=1:(length(iteration)-1)
    if iteration(i+1) - iteration(i) < 1
       gap=gap+iteration(i);
       plot(iteration(start:gap),K(start:gap,1),'linewidth',2);
       hold on;
       start=1+gap;
    end   
end

plot(iteration(start:end),K(start:end,1),'linewidth',2);
Q=title('Storage Design Study');
T=xlabel('Cycles')
Y=ylabel('Temperature Increase [K]');
U=legend('Diameter = 4 m','Diameter = 5 m','Diameter = 6 m','Diameter = 7 m','Diameter = 8 m')
Q.FontSize=14;
T.FontSize=14;
Y.FontSize=14;
U.FontSize=14;
hold off;
