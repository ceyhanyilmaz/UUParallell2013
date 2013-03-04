% Timing array:
% Rows: num fish
% Columns: num processors.

timing = [
0.000789    0.002772    0.033044    1.636971 
0.000979    0.001987    0.021890    1.140199    
0.006679    0.004934    0.020894    0.675296
0.017210    0.016930    0.026923    0.591835
0.046270    0.044872    0.045339    0.513553
0.069881    0.069573    0.069902    0.387441
];

nprocs = [1 3 9 23 49 73]; 
nelements = [1000 10000 100000 10000000];

% Possible plots:
figure(1);
plot(nprocs,timing(:,1), 'g-');

hold on;
plot(nprocs,timing(:,2), 'r--');
plot(nprocs,timing(:,3), 'b:');
legend('n = 1 000', 'n = 10 000', 'n = 100 000', 'Location','SouthEast');

hold off;
title('Comparison');
xlabel('Num threads'); 
ylabel('Time in seconds');

figure(2);
plot(nprocs,timing(:,4), 'g-');

hold on;
plot(nprocs,timing(:,4), 'r--');
legend('n = 1 000 000', 'n = 10 000 000', 'Location','NorthEast');

hold off;
title('Comparison');
xlabel('Num threads'); 
ylabel('Time in seconds');

