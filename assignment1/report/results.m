% 4 processors
tab1=[
100 0.0135580003
400 0.1406709999
600 0.3969019987
800 1.0913990010
1000 2.3293339983
1200 4.2203289997
1400 10.0848089997
];

clf;
xlabel('Time in seconds');
ylabel('Number of elements');

plot(tab(:,1),tab(:,2),'b');
hold on;


