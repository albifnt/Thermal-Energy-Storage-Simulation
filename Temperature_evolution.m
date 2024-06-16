% Remove
clear all

% Loading
F=load('Temperature_Fluid.txt');
S=load('Temperature_Solid.txt');
x=load('Space_Vector.txt');


% Simulation
for t=1:length(F(:,1))
 
     % update figure
     % If the visualisation proceeds to slowly,in  " mod(t,*)"  you have to
     % increase the star_variable
     if (mod(t,1)==0)
         
         plot(fliplr(F(t,:)),x,'b-',fliplr(S(t,:)),x,'r-','LineWidth',2)
         hold on;
         title(sprintf('Visualisation = %d', t));
         ylabel('height')
         xlabel('Temperature (K)')
         legend('Fluid','Solid');
         xlim([250 873]);
         hold off;
         
         drawnow
     end
end