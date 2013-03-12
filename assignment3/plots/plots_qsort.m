

timing = [
0.018303    0.996843    100.75
0.009575    0.534549    52.25
0.005738    0.285060    25.5
0.005011    0.192109    15.340202
0.007466    0.216006    10.325731
0.033227    0.307130    9.155618
0.006699    0.153766    9.016023
];

nprocs = [1 2 4 8 16 32 64]; 
nelements = [100000 1000000 10000000];

speedup = zeros(7,3);
speedup(:,1) = timing(1,1)./timing(:,1);
speedup(:,2) = timing(1,2)./timing(:,2);
speedup(:,3) = timing(1,3)./timing(:,3);


% Possible plots:
figure(1);
plot(nprocs,timing(:,1), 'r-');
title('Quick sort with OpenMP');
xlabel('Number of threads'); 
ylabel('Time in seconds');
legend('n = 100 000', 'Location','SouthEast');

figure(2);
plot(nprocs,timing(:,2), 'r-');
title('Quick sort with OpenMP');
xlabel('Number of threads'); 
ylabel('Time in seconds');
legend('n = 1 000 000', 'Location','SouthEast');

figure(3);
plot(nprocs,timing(:,3), 'r-');
title('Quick sort with OpenMP');
xlabel('Number of threads'); 
ylabel('Time in seconds');
legend('n = 10 000 000', 'Location','SouthEast');

figure(4);

plot(nprocs,speedup(:,1), 'g--');
hold on;
plot(nprocs,speedup(:,2), 'r:');
plot(nprocs,speedup(:,3), 'b-');

title('Quick sort speedup with OpenMP');
xlabel('Number of threads'); 
ylabel('Speedup');
legend('n = 100 000', 'n = 1 000 000', 'n = 10 000 000', 'Location','SouthEast');


