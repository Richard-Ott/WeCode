k = 2.3028e-6;

dt = 1e2;
t = 1:dt:1e6;

P = 30;

N0 = 2e4;

Nt = nan(1,length(t)); N_nodecay = Nt;
Nt(1) = N0;  N_nodecay(1) = N0;
for i = 2:length(t)
    Nt(i) = Nt(i-1) + P*dt;
    Nt(i) = Nt(i)*exp(-k*dt);
    
    N_nodecay(i) =  N_nodecay(i-1) + P*dt;
end

% analytical solution
c = N0 - P/k;
Ntana = c .*exp(-k.*t) + P/k;


%% plot
figure()
plot(t,Nt);
hold on
plot(t,Ntana);
plot(t,N_nodecay);
legend("numerical","analytical","no_decay")