%% README
%integrate & fire neuron:
%part a,b,c generate spikes of a neuron using poisson process.
%part d,e,f investigates the effect of delleting all but kth spikes on the
%statistical parameters (renewal process)
%part g add refractory period effect
%the leaky integrate and fire neuron:
%at firstgenerting the membrane voltage by LIF model with constant current
%input then change the imput to time varing one using EPSC kernel and
%investigate the statistical properties 
%%  
clear all
clc
%% Integration & Fire Neuron
dt = 1 ; % ms
FR = 100 ; 
T_end = 1000 ; % ms
T = 0:dt:T_end-dt ; % ms
N_trail = 500 ;
threshold = FR*dt/1000 ;
spike = rand(N_trail , T_end) < threshold ;
%spike = double(spike) ;
%% part a
plotRaster(spike , T) ;
xlim([0 , T_end]) ;
xlabel('Time(ms)') ;
ylabel('Trial Number') ;
title('Raster Plot of Spike Train(Poisson Process)') ;
%% part a limited
plotRaster(spike(1:100 , 1:500) , T(1:500)) ;
%xlim([0 , T_end]) ;
xlabel('Time(ms)') ;
ylabel('Trial Number') ;
title('Raster Plot of Spike Train(Poisson Process)') ;
%% part b
spike_count = sum(spike , 2) ;
a = min(spike_count) ;
b = max(spike_count) ;
n = b-a+1 ;
histogram(spike_count , n , 'Normalization' , 'pdf');
hist = histcounts(spike_count , n , 'Normalization' , 'pdf');
hold on ;
x = linspace(a , b , n) ;
pdf = poisspdf(x , FR) ;
plot(x , pdf) ;
xlim([a , b]) ;
title('Spike Count Histogram vs Poisson Distribution PDF') ;
legend('Spike Count Histogram' , 'Poisson Distribution PDF') ;
%% part c
ISI = [] ;
for i = 1:N_trail
    dis = find(spike(i , :)) ;
    dis = diff(dis) ;
    ISI = [ISI , dis] ;
end
ISI = ISI./1000 ; 
%% part c ploting
a = min(ISI) ;
b = max(ISI) ;
n = (b-a+1)*1000 ;
n = 150 ;
histogram(ISI , n , 'Normalization' , 'pdf');
hist = histcounts(ISI , n , 'Normalization' , 'pdf');
hold on ;
x = linspace(a , b , n) ;
pdf = exppdf(x , 1/(FR)) ;
%plot(x , pdf/(11*pdf(1))) ;
%plot(x , pdf./(sum(pdf))) ;
plot(x , 1.8*pdf) ;
title('Inter-Spike Interval Histogram vs Exponential Distribution') ;
legend('Inter-Spike Interval Histogram' , 'Exponential Distribution') ;
xlabel('time(s)') ;
%% repeat part a
k = [4 , 8 , 20] ;
s = spike ;
T_end = 10000 ;
T = 0:dt:T_end-dt ;
spike = spike_generator(FR , dt , N_trail , T_end) ;
spike_new1 = integrate(spike , k(1)) ;
spike_new2 = integrate(spike , k(2)) ;
spike_new3 = integrate(spike , k(3)) ;
%%
figure ;
plotRaster(spike_new1(1:100 , 1:500) , T(1:500)) ;
xlabel('Time(ms)') ;
ylabel('Trial Number') ;
title('Raster Plot For Poisson Process(Deleting Spike Procedure (k=4))') ;
figure ;
plotRaster(spike_new2(1:100 , 1:500) , T(1:500)) ;
xlabel('Time(ms)') ;
ylabel('Trial Number') ;
title('Raster Plot For Poisson Process(Deleting Spike Procedure (k=8))') ;
figure ;
plotRaster(spike_new3(1:100 , 1:500) , T(1:500)) ;
xlabel('Time(ms)') ;
ylabel('Trial Number') ;
title('Raster Plot For Poisson Process(Deleting Spike Procedure (k=20))') ;
%% repeat part b
spike_count_new1 = sum(spike_new1 , 2) ;
a = min(spike_count_new1) ;
b = max(spike_count_new1) ;
n = b-a+1 ;
n1 = 30 ;
figure ;
histogram(spike_count_new1 , n1 , 'Normalization' , 'pdf');
hist_new = histcounts(spike_count_new1 , n1 , 'Normalization' , 'probability');
hold on ;
x_new = linspace(a , b , n) ;
y = linspace(k(1)*a , k(1)*b , n) ;
pdf = k(1)*poisspdf(y , T_end*FR/1000) ;
plot(x_new , pdf) ;
title('Spike Count Histogram(Deleting Spike Procedure k=4) vs Poisson Distribution') ;
legend('Spike Count Histogram' , 'Poisson Distribution') ;
xlabel('Number of Spikes') ;
%%
spike_count_new2 = sum(spike_new2 , 2) ;
a = min(spike_count_new2) ;
b = max(spike_count_new2) ;
n = b-a+1 ;
n1 = 30 ;
figure ;
histogram(spike_count_new2 , n1 , 'Normalization' , 'probability');
hist_new = histcounts(spike_count_new2 , n1 , 'Normalization' , 'probability');
hold on ;
x_new = linspace(a , b , n) ;
y = linspace(k(2)*a , k(2)*b , n) ;
pdf = k(2)*poisspdf(y , T_end*FR/1000) ;
plot(x_new , pdf) ;
title('Spike Count Histogram(Deleting Spike Procedure k=8) vs Poisson Distribution') ;
legend('Spike Count Histogram' , 'Poisson Distribution') ;
xlabel('Number of Spikes') ;
%%
spike_count_new3 = sum(spike_new3 , 2) ;
a = min(spike_count_new3) ;
b = max(spike_count_new3) ;
n = b-a+1 ;
n1 = 30 ;
figure ;
histogram(spike_count_new3 , n1 , 'Normalization' , 'probability');
hist_new = histcounts(spike_count_new3 , n1 , 'Normalization' , 'probability');
hold on ;
x_new = linspace(a , b , n) ;
y = linspace(k(3)*a , k(3)*b , n) ;
pdf = k(3)*poisspdf(y , T_end*FR/1000) ;
plot(x_new , pdf) ;
title('Spike Count Histogram(Deleting Spike Procedure k=20) vs Poisson Distribution') ;
legend('Spike Count Histogram' , 'Poisson Distribution') ;
xlabel('Number of Spikes') ;
%% repeat part c 
ISI_new1 = ISI_compute(spike_new1) ;
ISI_new2 = ISI_compute(spike_new2) ;
ISI_new3 = ISI_compute(spike_new3) ; 
%% repeat part c ploting
a = min(ISI_new1) ;
b = max(ISI_new1) ;
n = (b-a+1)*1000 ;
n = 160 ;
histogram(ISI_new1 , n , 'Normalization' , 'pdf');
hist = histcounts(ISI_new1 , n , 'Normalization' , 'pdf');
hold on ;
x = linspace(a , b , n) ;
pdf = exppdf(x , 1/(FR)) ;
pdf = ((FR.*x).^(k(1)-1)).*(pdf)./(factorial(k(1)-1)) ;
plot(x , pdf) ;
xlabel('\Deltat(s)') ;
title('Inter-Spike Interval Histogram(Deleting Spike Procedure k=4) vs Erlang Distribution') ;
legend('Inter-Spike Interval Histogram'  , 'Erlang Distribution') ;
%%
a = min(ISI_new2) ;
b = max(ISI_new2) ;
n = (b-a+1)*1000 ;
n = 100 ;
figure
histogram(ISI_new2 , n , 'Normalization' , 'pdf');
hist = histcounts(ISI_new2 , n , 'Normalization' , 'pdf');
hold on ;
x = linspace(a , b , n) ;
pdf = exppdf(x , 1/(FR)) ;
pdf = ((FR.*x).^(k(2)-1)).*(pdf)./(factorial(k(2)-1)) ;
plot(x , pdf) ;
xlabel('\Deltat(s)') ;
title('Inter-Spike Interval Histogram(Deleting Spike Procedure k=8) vs Erlang Distribution') ;
legend('Inter-Spike Interval Histogram'  , 'Erlang Distribution') ;
%%
a = min(ISI_new3) ;
b = max(ISI_new3) ;
n = (b-a+1)*1000 ;
n = 120 ;
histogram(ISI_new3 , n , 'Normalization' , 'pdf');
hist = histcounts(ISI_new3 , n , 'Normalization' , 'pdf');
hold on ;
x = linspace(a , b , n) ;
pdf = exppdf(x , 1/(FR)) ;
pdf = ((FR.*x).^(k(3)-1)).*(pdf)./(factorial(k(3)-1)) ;
plot(x , pdf) ;
xlabel('\Deltat(s)') ;
title('Inter-Spike Interval Histogram(Deleting Spike Procedure k=20) vs Erlang Distribution') ;
legend('Inter-Spike Interval Histogram'  , 'Erlang Distribution') ;
%% part d
K = 1:1:50 ;
CV = [] ;
CV_ref = [] ;
for i = K 
    S = integrate(spike , i) ;
    S = ISI_compute(S) ;
    CV = [CV , std(S)/mean(S)] ;
    CV_ref = [CV_ref , 1/sqrt(i)] ;
end
stem(K , CV) ;
hold on ;
stem(K , CV_ref) ;
%%
final_CV = zeros(1 , 50) ;
for j = 1:100
    spike = spike_generator(FR , dt , N_trail , T_end) ;
    K = 1:1:50 ;
    CV = [] ;
    CV_ref = [] ;
    for i = K
        S = integrate(spike , i) ;
        S = ISI_compute(S) ;
        CV = [CV , std(S)/mean(S)] ;
    end
    final_CV = final_CV + CV ;
end
final_CV = final_CV/100 ;
%% Part d again
final_CV1 = zeros(1 , 50) ;
FR = 100 ;
CV2 = [] ;
dt = 0.04 ;
for j = 1:100
    spike5 = spike_generator(FR , dt , N_trail , T_end) ;
    S = ISI_compute(spike5) ;
    CV2 = [CV2 , std(S)/mean(S)] ;
end
figure ;
yline(mean(CV2)) ;
hold on ;
scatter(1:100 , CV2 , 'filled') ;
xlabel('Experiment') ;
ylabel('CV') ;
title('Computing CV for generated spike trian') ;
legend('average CV = ' , 'CV for each experiment') ;
%%
dt = 1 ;
%%
plot(K , final_CV) ;
hold on ;
plot(K , 1./(sqrt(K))) ;
title('Computed CV vs 1/\surdk') ;
legend('Computed CV' , '1/\surdk') ;
xlabel('k') ;
%plot(K , -final_CV+1./(sqrt(K)))
%% part g
CV1 = CV ;
%FR = 50:1:300 ;
t1 = 0.001*(1:0.5:45) ;
N = [1 , 4 , 51] ;
ref = [1 , 2 , 4] ;
dt = 1 ;
T = 5000 ;
trail = 50 ;
CV = zeros(length(N) , length(t1) , length(ref)) ;
for i = t1
    spike = spike_generator(1/i , dt , trail , T) ;
    for j = ref
        spike1 = refractory_compute(j , spike) ;
        for r = N
            spike2 = integrate(spike1 , r) ;
            CV(find(N == r) , find(t1 == i) , find(ref == j)) = CV_compute(spike2) ;
        end
    end
end
a1 = zeros(1 , length(t1)) ;
b1 = a1 ;
b2 = b1 ;
b3 = b1 ;
c1 = a1 ;
c2 = a1 ; 
c3 = a1 ;
a1 = CV(1 , : , 1) ;
b1 = CV(1 , : , 2) ;
c1 = CV(1 , : , 3) ;
a2 = zeros(1 , length(t1)) ;
a2 = CV(2 , : , 1) ;
b2 = CV(2 , : , 2) ;
c2 = CV(2 , : , 3) ;
a3 = zeros(1 , length(t1)) ;
a3 = CV(3 , : , 1) ;
b3 = CV(3 , : , 2) ;
c3 = CV(3 , : , 3) ;
%%
scatter(1000*t1 , a1 , 15 , [1 1 0] , 'filled') ;
hold on ;
scatter(1000*t1 , b1 , 15 , [1 1 0] , 'filled') ;
scatter(1000*t1 , c1 , 15 , [1 1 0] , 'filled') ;
yline(1/sqrt(N(1)) , 'color' , [1 1 0]) ;
scatter(1000*t1 , a2 , 15 , [0 1 0] , 'filled') ;
scatter(1000*t1 , b2 , 15 , [0 1 0] , 'filled') ;
scatter(1000*t1 , c2 , 15 , [0 1 0] , 'filled') ;
yline(1/sqrt(N(2)) , 'color' , [0 1 0]) ;
scatter(1000*t1 , a3 , 15 , [0 0 1] , 'filled') ;
scatter(1000*t1 , b3 , 15 , [0 0 1] , 'filled') ;
scatter(1000*t1 , c3 , 15 , [0 0 1] , 'filled') ;
yline(1/sqrt(N(3)) , 'color' , [0 0 1]) ;
title('Comparison of CV for integrator models') ;
xlabel('\Deltat(ms)') ;
ylabel('CV') ;
%% Leaky Integrate & Fire Neuron
clear all ;
clc ;
%% part a
T_end = 100 ; % ms
dt = 0.01 ; % ms
T = 0:dt:T_end ; 
Vth = 15 ; % mV
Vmax = 60 ;
Vrest = 0 ; % mV
IR = 20 ; % mV
Is = IR*ones(1 , T_end/dt + 1) ; 
Vm1 = zeros(1 , T_end/dt + 1) ;
r = 40 ;
tau = 10 ;
p = 0 ;
a = 2:1:T_end/dt ;
refractory = r*dt ;
for i = 2:1:T_end/dt+1
    if ((Vm1(i-1) > Vth)||(Vm1(i-1) == Vth))&&(p == 0)
        a = linspace(Vmax , Vrest , refractory/dt+1) ;
        Vm1(i:min(i+refractory/dt , T_end/dt + 1)) = a(1:min(refractory/dt+1 , T_end/dt + 2-i)) ;
        p = 1 ;
    elseif (p == 0)||(Vm1(i-1) == Vrest)
        Vm1(i) = Vm1(i-1) -dt*Vm1(i-1)/tau + dt*IR/tau ;
        p = 0 ;
    end
end
plot(T , Vm1) ;
%% part a
tau_all = [2 , 5 , 10] ;
refractory_all = dt.*[5 , 10 , 15] ;
Vm = Vrest*ones(length(tau_all) , length(refractory_all) , T_end/dt + 1) ;
for tau = tau_all
    for refractory = refractory_all
        Vm1 = Vrest*ones(1 , T_end/dt + 1) ;
        for i = 2:1:T_end/dt+1
            if ((Vm1(i-1) > Vth)||(Vm1(i-1) == Vth))&&(p == 0)
                a = linspace(Vmax , Vrest , refractory/dt+1) ;
                Vm1(i:min(i+refractory/dt , T_end/dt + 1)) = a(1:min(refractory/dt+1 , T_end/dt + 2-i)) ;
                p = 1 ;
            elseif (p == 0)||(Vm1(i-1) == Vrest)
                Vm1(i) = Vm1(i-1) -dt*Vm1(i-1)/tau + dt*IR/tau ;
                p = 0 ;
            end
        end
%         Vm(find(tau_all == tau) , find(refractory_all == refractory) , :) = Vm1 ;
%         plot(T , Vm1) ;
%         hold on ;
    end
end
%%
figure ;
for i = 1:1:length(refractory_all)
     a = zeros(1 , T_end/dt + 1) ;
     a(1 , :) = Vm(tau_all(1) , i , :) ;
     plot(T , a) ;
     hold on ;
end
title('Leaky Integrate and Fire(comparing refractory period)') ;
xlabel('time(ms)') ;
ylabel('potential(mV)') ;
legend('t_0 = 5' , 't_0 = 10' ,'t_0 = 15') ;
%%
figure ;
for i = 1:1:length(tau_all)
     a = zeros(1 , T_end/dt + 1) ;
     a(1 , :) = Vm(i , 2 , :) ;
     plot(T , a) ;
     hold on ;
end
title('Leaky Integrate and Fire(comparing time constant(\tau_m))') ;
xlabel('time(ms)') ;
ylabel('potential(mV)') ;
legend('\tau_m = 2' , '\tau_m = 5' ,'\tau_m = 10') ;
%% part c
FR = 100 ;
trail = 50 ;
T_end = 1000 ; % ms
dt = 0.01 ; % ms
T = 0:dt:T_end-dt ;
tau_peak = 1 ;
Vth = 10 ;
Vmax = 40 ;
spike = spike_generator(FR , dt , trail , T_end/dt) ;
t1 = 0:dt:5 ;
w = 1.5 ;
kernel = (w).*t1.*exp(-t1./tau_peak) ;
% plot(t1 , kernel) ;
% xlabel('Time(ms)') ;
% title('The kernel') ;
input = [] ;
for i = 1:1:trail
    a = conv(spike(i , :) , kernel) ;
    a = a(1 , 1:T_end/dt) ;
    input = [input ; a] ;
end
%%
input_prime = sum(input) ;
for i = 1:1:T_end/dt
    if(input_prime(i) > Vth)||(input_prime(i) == Vth)
       input_prime(i) = Vmax ;
    end
end
plot(T , input_prime) ;
xlabel('Time(ms)') ;
ylabel('Potential(mV)') ;
title('stimulating neuron by a time-varying input current') ;
%% part d
n = 20 ;
inhibit = randi([1 50] , 1 , (n/100)*trail) ;
input_new = input ;
for i = inhibit
    input_new(i , :) = -spike(i , :) ;
end
input_prime = sum(input_new) ;
%
for i = 1:1:T_end/dt
    if(input_prime(i) > Vth)||(input_prime(i) == Vth)
       input_prime(i) = Vmax ;
    end
end
plot(T , input_prime) ;
xlabel('Time(ms)') ;
ylabel('Potential(mV)') ;
title('stimulating neuron by a time-varying input current(EPSP & IPSP)') ;
%% part d details
inhibitory = 0:5:40 ;
CV_IPSP = [] ;
for n = inhibitory
    inhibit = randi([1 50] , 1 , round((n/100)*trail)) ;
    input_new = input ;
    for i = inhibit
        input_new(i , :) = -spike(i , :) ;
    end
    input_prime = sum(input_new) ;
    a = zeros(1 , T_end/dt) ;
    for i = 1:1:T_end/dt
        if(input_prime(i) > Vth)||(input_prime(i) == Vth)
            input_prime(i) = Vmax ;
            a(i) = 1 ;
        else
            a(i) = 0 ;
        end
    end
    CV_IPSP = [CV_IPSP , std(a)/mean(a)] ;
    plot(T , input_prime) ;
    hold on ;
end
xlabel('Time(ms)') ;
ylabel('Potential(mV)') ;
title('stimulating neuron by a time-varying input current(EPSP & IPSP)') ;
hold off ;
figure ;
plot(inhibitory , CV_IPSP) ;
title('CV') ;
xlabel('percentage of inhibitory neurons') ;
%% part e 
FR = 100 ;
T_end = 1000 ; % ms
dt = 0.01 ; % ms
T = 0:dt:T_end-dt ;
N = 60 ;
n = 10 ;
M = round(N*n/100) ;
spike_set = spike_generator(FR , dt , N , T_end/dt) ;
refractory = 0.04 ; %ms
r = refractory/dt ;
new_spike = zeros(1 , T_end/dt) ;
D = 0.2 ; %ms
D_len = D/dt ;
for i = 1:1:T_end/dt
    a = sum(sum(spike_set(: , i:min(i+D_len , end)))) ;
    new_spike(1 , i) = double(a > M-1) ;
end
figure ;
plot(T , new_spike) ;
hold on ; 
plot(T , spike_set(1 , :)) ;
title('Generated spike') ;
xlabel('time(ms)') ;
legend('post' , 'pre')
% figure ;
% plotRaster(spike_set, (T+dt)/dt) ;
%% dependence to N
N_all = 60:10:150 ;
CV_N = [] ;
    for j = N_all
        a = coincidence1(FR , dt , j , T_end , D , 20 , 0) ;
        CV_N = [CV_N , std(a)/mean(a)] ;
    end
plot(CV_N)
%% dependence to M

%% dependence to D
D_all = 0.01:0.05:1 ;
CV_coincidence = [] ;
percent_inhibit = 0 ;
for i = D_all
    a = coincidence(FR , dt , N , T_end , i , n , percent_inhibit) ;
    CV_coincidence = [CV_coincidence , std(a)/mean(a)] ;
end
%%
plot(D_all , CV_coincidence) ;
title('dependece to D') ;
xlabel('length of the window(ms)') ;
ylabel('CV') ;
%% part f
percent_inhibit_all = 0:2:50 ;
CV_inhibit = [] ;
vector = [] ;
N = 100 ;
D = 0.5 ;
for i = percent_inhibit_all
    a = coincidence(FR , dt , N , T_end , D , n , i) ;
    CV_inhibit = [CV_inhibit , std(a)/mean(a)] ;
    if mean(a) ~= 0
        vector = [vector , i] ;
    end
end
%%
plot(vector , CV_inhibit) ;
title('dependece to the percentage of inhibitory neurons') ;
xlabel('the percentage') ;
ylabel('CV') ;
%% functions 
function [] = plotRaster(spikeMat, tVec)
hold all;
for trialCount = 1:size(spikeMat,1)
    spikePos = tVec(spikeMat(trialCount, :));
    for spikeCount = 1:length(spikePos)
        plot([spikePos(spikeCount) spikePos(spikeCount)], ...
        [trialCount-0.4 trialCount+0.4], 'k');
    end
end
ylim([0 size(spikeMat, 1)+1]);
end
function [spike_new] = integrate(spike , k) 
[a , b] = size(spike) ;
spike_new = zeros(a , b) ;
for i = 1:a
    t = 0 ;
    for j = 1:b
        if spike(i , j) == 1
            t = t+1 ;
        end
        if t == k
            spike_new(i , j) = 1 ;
            t = 0 ;
        end
    end
end
spike_new = logical(spike_new) ;
end
function [ISI_new] = ISI_compute(spike)
ISI_new = [] ;
[a b] = size(spike) ;
for i = 1:a
    dis = find(spike(i , :)) ;
    dis = diff(dis) ;
    ISI_new = [ISI_new , dis] ;
end
ISI_new = ISI_new -1 ; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ISI_new = ISI_new./1000 ; 
end
function [spike] = spike_generator(FR , dt , trail , T) 
threshold = FR*dt/1000 ;
spike = rand(trail , T) < threshold ;
spike = double(spike) ;
end
function [spike] = spike_generator1(FR , dt , trail , T) 
threshold = FR*dt/1000 ;
spike = rand(trail , T) < threshold ;
spike = double(spike) ;
end
function [spike_new] = refractory_compute(delta_t , spike)
[a , b] = size(spike) ;
spike_new = spike ;
for i = 1:a
    loc = find(spike(i , :)) ;
    for j = loc
        if spike_new(i , j) == 1
            spike_new(i , j+1:min(j+delta_t , b)) = 0 ;
        end
    end
end
end
function [CV] = CV_compute(spike) 
ISI = ISI_compute(spike) ;
CV = std(ISI)/mean(ISI) ;
end
function [new_spike] = coincidence(FR , dt , N , T_end , D , n , percent_inhibit) 
M = round(N*n/100) ;
spike_set = spike_generator(FR , dt , N , T_end/dt) ;
refractory = 0.04 ; %ms
r = refractory/dt ;
new_spike = zeros(1 , T_end/dt) ;
D_len = D/dt ;
%
spike_new = spike_set ;
if percent_inhibit ~= 0
    inhibit = randi([1 50] , 1 , (percent_inhibit/100)*N) ;
    for i = inhibit
        spike_new(i , :) = -spike_set(i , :) ;
    end
end
%
for i = 1:1:T_end/dt
    a = sum(sum(spike_new(: , i:min(i+D_len , end)))) ;
    new_spike(1 , i) = double(a > M-1) ;
end
%plot(T , new_spike) ;
end

function [new_spike] = coincidence1(FR , dt , N , T_end , D , number , number_inhibit) 
M = number ;
spike_set = spike_generator(FR , dt , N , T_end/dt) ;
refractory = 0.04 ; %ms
r = refractory/dt ;
new_spike = zeros(1 , T_end/dt) ;
D_len = D/dt ;
%
spike_new = spike_set ;
if number_inhibit ~= 0
    inhibit = randi([1 50] , 1 , (number_inhibit)) ;
    for i = inhibit
        spike_new(i , :) = -spike_set(i , :) ;
    end
end
%
for i = 1:1:T_end/dt
    a = sum(sum(spike_new(: , i:min(i+D_len , end)))) ;
    new_spike(1 , i) = double(a > M-1) ;
end
%plot(T , new_spike) ;
end