% 4 processors
timing=[0.0046380000 0.0044649999 0.0317310000 0.0537430000
0.3055450000 0.0940820000 0.0937680000 0.1174500000
2.6802310000 0.6996980000 0.3857760000 0.2979840000
7.7199790000 1.4716900000 0.6704180000 0.4979620000
13.3761460000 2.8211280000 1.3280860001 1.1920820000
25.1697130000 4.2472240000 2.3558620000 1.5529640000
];


nelements = [100 400 800 1000 1200 1440];

for i = 1:length(nelements)
  speedup(i,1) = timing(i,1) ./ timing(i,2);
  speedup(i,2) = timing(i,1) ./ timing(i,3);
  speedup(i,3) = timing(i,1) ./ timing(i,4);
end

% Possible plots:
figure(1)
plot(nelements,speedup)
title('Speed plot')
xlabel('Num elements')
ylabel('Speedup')
legend('proc = 4','proc = 9','proc = 16')
