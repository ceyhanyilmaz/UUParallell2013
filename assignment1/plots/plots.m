% Timing array:
% Rows: num fish
% Columns: num processors.


timing = [
0.0123180002  0.0123180002
0.5687389988  0.5687389988
2.4870360009  2.1265500002
5.8614530005  5.0209519994
13.2836399991 6.0105939992
23.6162660010 11.3004159983
40.4612520002 19.2244499996
];

nprocs = [1 1 4 4 4 4 4]; 
nelements = [100 400 600 800 1000 1200 1400];

for i = 1:length(nelements)
   speedup(:,i) = timing(i,1) ./ timing(i,2);
   efficiency(:,i) = speedup(:,i)/nprocs(i);
end

% Possible plots:
figure(1)
plot(nelements,speedup)
title('Speed plot')
xlabel('Num elements') 
ylabel('Speedup')

figure(2)
plot(nelements,efficiency)
title('Efficiency plot')
xlabel('Num elements') 
ylabel('Efficiency')

