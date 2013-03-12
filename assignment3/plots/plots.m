

timing = [
1.994249    19.870712    65.174862    143.5       308.75 
1.021339    12.072147    40.806183    94.852069   184.935649
0.493782    6.633971     28.322395    66.75       146.5
0.259408    4            19.450961    64.25       129.5 
0.187886    2.808593     19.182807    63.5        123.830462 
22.174041   40.024980    66.444143    116.844250  198.75
1.096868    3.058758     18.25        66.5        136.25 
];

nprocs = [1 2 4 8 16 32 64]; 
nelements = [1000 2000 3000 4000 5000];

speedup = zeros(7,5);
speedup(:,1) = timing(1,1)./timing(:,1);
speedup(:,2) = timing(1,2)./timing(:,2);
speedup(:,3) = timing(1,3)./timing(:,3);
speedup(:,4) = timing(1,4)./timing(:,4);
speedup(:,5) = timing(1,5)./timing(:,5);

% Possible plots:
figure(1);
plot(nprocs,timing(:,1), 'r-');
title('Gram-Schmidt');
xlabel('Number of threads'); 
ylabel('Time in seconds');
legend('n = 1 000', 'Location','SouthEast');

figure(2);
plot(nprocs,timing(:,2), 'r-');
title('Gram-Schmidt');
xlabel('Number of threads'); 
ylabel('Time in seconds');
legend('n = 2 000', 'Location','SouthEast');

figure(3);
plot(nprocs,timing(:,3), 'r-');
title('Gram-Schmidt');
xlabel('Number of threads'); 
ylabel('Time in seconds');
legend('n = 3 000', 'Location','SouthEast');

figure(4);
plot(nprocs,timing(:,4), 'r-');
title('Gram-Schmidt');
xlabel('Number of threads'); 
ylabel('Time in seconds');
legend('n = 4 000', 'Location','SouthEast');

figure(5);
plot(nprocs,timing(:,5), 'r-');
title('Gram-Schmidt');
xlabel('Number of threads'); 
ylabel('Time in seconds');
legend('n = 5 000', 'Location','SouthEast');

figure(6);
plot(nprocs,speedup(:,1), 'g--');
hold on;
plot(nprocs,speedup(:,2), 'r:');
plot(nprocs,speedup(:,3), 'b-');
plot(nprocs,speedup(:,4), 'b-.');
plot(nprocs,speedup(:,5), 'r-');
title('Gram-Schmidt speedup');
xlabel('Number of threads'); 
ylabel('Speedup');
legend('n = 1000', 'n = 2000', 'n = 3000', 'n = 4000', 'n = 5000', 'Location','SouthEast');

