%Remove
clear all

% Loading
STATUS=load('System_Status.txt');
TIME=load('Time_vector.txt');

% Plot
figure(1)
plot(TIME,STATUS,'Linewidth',2)
xlabel('Time (s)')
ylabel('Status')
yticks([1 2 3 4])
%xticks([20000 40000 60000 80000])
