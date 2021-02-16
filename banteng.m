% program for analysis of banteng population projection

% ----------------------------------------
% base model - rates from Choquenot (1993)
% giving stability
% ----------------------------------------

format long g

% define max age
age_max=17;

% define survival - fitted function (Choquenot 1993)
%s0=0.7402; % fitted from Choquenot (1993)
s0=0.55; % to give stability
s1=0.9278;
s2=0.9346;
s3=0.9370;
s4=0.9382;
s5=0.9389;
s6=0.9394;
s7=0.9397;
s8=0.9400;
s9=0.9402;
s10=0.9403;
s11=0.9405;
s12=0.9406;
s13=0.9407;
s14=0.9408;
s15=0.9408;
s16=0.9409;
s17=0.9409;

% define fecundity parameters (Choquenot 1993)
%m1=0;
%m2=0;
%m3=0;
%m4=0.220;
%m5=0.351;
%m6=0.500;
%m7=0.450;
%m8=0.200;
%m9=0.325;
%m10=0.125;
%m11=0.225;
%m12=m11;
%m13=m11;
%m14=m11;
%m15=m11;
%m16=m11;

% fitted fecundity parameters
m0=0.0000;
m1=0.0001;
m2=0.0133;
m3=0.0947;
m4=0.2303;
m5=0.3440;
m6=0.3964;
m7=0.3921;
m8=0.3525;
m9=0.2979;
m10=0.2414;
m11=0.1903;
m12=0.1471;
m13=0.1122;
m14=0.0849;
m15=0.0638;
m16=0.0478;
m17=0.0358;

primi=4; % age at first reproduction (females)

% initial population sizes
% Nmin=4000; % J. Christopherson (2004)
% Nmax=5000; % J. Christopherson (2004)
Nmin=7000; % K. Saalfeld (2004)
Nmax=8000; % K. Saalfeld (2004)
Navg = mean([Nmin Nmax]);

% Use 2 scenarios - pop_init at 4000 and 8000
N1 = 4000;
N2 = 8000;

% Density estimates
Area_Cob = 220000; % Area of Cobourg Peninsula in hectares
D_avg1 = N1/(Area_Cob*0.01);
D_avg2 = N2/(Area_Cob*0.01);

% sex ratios
sr = 0.5; % overall sex ratio
x = 0.5; % foetal sex ratio (Choquenot 1993)

%total females in population
f1=N1*sr;
f2=N2*sr;

%total males in population
mal1=N1*(1-sr);
mal2=N2*(1-sr);

%Initial population size vector
N=N1;

% the normal matrix
a=[0 s0*m1 s0*m2 s0*m3 s0*m4 s0*m5 s0*m6 s0*m7 s0*m8 s0*m9 s0*m10 s0*m11 s0*m12 s0*m13 s0*m14 s0*m15 s0*m16 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
    s1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
    0 s2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
    0 0 s3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
    0 0 0 s4 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
    0 0 0 0 s5 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
    0 0 0 0 0 s6 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
    0 0 0 0 0 0 s7 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 s8 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 0 s9 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 0 0 s10 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 0 0 0 s11 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 0 0 0 0 s12 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 0 0 0 0 0 s13 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 0 0 0 0 0 0 s14 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 s15 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 s16 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
    0 s0*m1 s0*m2 s0*m3 s0*m4 s0*m5 s0*m6 s0*m7 s0*m8 s0*m9 s0*m10 s0*m11 s0*m12 s0*m13 s0*m14 s0*m15 s0*m16 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 s1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 s2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 s3 0 0 0 0 0 0 0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 s4 0 0 0 0 0 0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 s5 0 0 0 0 0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 s6 0 0 0 0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 s7 0 0 0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 s8 0 0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 s9 0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 s10 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 s11 0 0 0 0 0 0
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 s12 0 0 0 0 0
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 s13 0 0 0 0
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 s14 0 0 0
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 s15 0 0
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 s16 0];

% eigenvalues & eigenvectors
[w,d] = eig(a);
v=conj(inv(w));

%lambda
lambda=diag(d);

% max lambdas
[maxlambda,imax]=max(diag(d));

%logmaxlambda=exponential rate of increase (r)
r=log(maxlambda);

% stable age distribution
w=w(:,imax);

% reproductive values
v=real(v(imax,:))';

% sensitivities and elasticities
senmat=(v*w')/(v'*w);
emat=(a.*senmat)/maxlambda;

%damping ratios
rho=lambda(1)./abs(lambda);

% periods of oscillation
period=2*pi./angle(lambda);

% stable size distribution (scaled to sum to 1)
ssd=w(:,1)/sum(w(:,1));

% size of matrix
k=size(a,1);

% ssd classes
    % female
    ssd_juvf = sum(ssd(1:primi-1,1));
    ssd_adf = sum(ssd(primi:age_max,1));

    % male
    ssd_juvm = sum(ssd(age_max+1:age_max+primi-1,1));
    ssd_adm = sum(ssd(age_max+primi:k,1));
    
% reproductive value (scaled so stage 1 = 1)
reprovalue=real(v(1,:)')/real(v(1,1));

% pick initial vectors
n=ones(k,1);
n=ssd*N;

% age vector
age_vec=zeros(age_max,1);
for j=1:age_max-1
    age_vec(j+1,1)=j;
end % j loop
age_vec=age_vec+1;

%Calculate Quasi-extinction times
thresh = 50; % Define quasi-exinction threshold

Q=(log(thresh/sum(n)))/log(maxlambda);
if Q < 0;
    Q='infinity';
else Q=Q;
end

% do a simulation
% first specify the initial condition and length of simulation
tlimit=10;

%set population size year step vector
pop_vec=ones(tlimit+1,1);
pop_vec(1,1)=sum(n);

%set year step vector
yr_vec=ones(tlimit+1,1);
for c=1:tlimit
    yr_vec(c+1,1)=c;
    yr_vec(1,1)=0;
end % c loop

% assume stable-age distribution
N_stable = ssd*N;

%then iterate
for i=1:tlimit;
    n=a*n;
    s1_vec(i+1,1)=(s1*(1-s1))/n(1,1);
    s2_vec(i+1,1)=(s2*(1-s2))/n(2,1);
    pop_vec(i+1,1)=(sum(n));    
end

log_pop_vec=log10(pop_vec);

%total population size after 'tlimit' years
pop_st=N
pop_end=sum(n)

%total population size after 'tlimit' years
N
tlimit
pop_end=sum(n)
maxlambda
r
Q
format short

% construct survival vector for display
surv=ones(age_max,1);
surv(1,1)=s0;
surv(2,1)=s1;
surv(3,1)=s2;
surv(4,1)=s3;
surv(5,1)=s4;
surv(6,1)=s5;
surv(7,1)=s6;
surv(8,1)=s7;
surv(9,1)=s8;
surv(10,1)=s9;
surv(11,1)=s10;
surv(12,1)=s11;
surv(13,1)=s12;
surv(14,1)=s13;
surv(15,1)=s14;
surv(16,1)=s15;
surv(17,1)=s16;

% construct fertility vector for display
fert=ones(age_max,1);
fert(1,1)=m1;
fert(2,1)=m2;
fert(3,1)=m3;
fert(4,1)=m4;
fert(5,1)=m5;
fert(6,1)=m6;
fert(7,1)=m7;
fert(8,1)=m8;
fert(9,1)=m9;
fert(10,1)=m10;
fert(11,1)=m11;
fert(12,1)=m12;
fert(13,1)=m13;
fert(14,1)=m14;
fert(15,1)=m15;
fert(16,1)=m16;
fert(17,1)=m16;

% continue displays
%a
ssd_juv = ssd_juvf + ssd_juvm
ssd_ad = ssd_adf + ssd_adm
%emat

%Make density independent plots
subplot(2,2,1), plot(yr_vec,pop_vec);
axis square;
axis([-0.5 tlimit+1 0 (max(pop_vec)+(0.25*(max(pop_vec))))]);
xlabel('year');
ylabel('N females');
subplot(2,2,2), plot(yr_vec,log_pop_vec);
axis square;
axis([-0.5 tlimit+1 0 (max(log_pop_vec)+(0.25*(max(log_pop_vec))))]);
xlabel('year');
ylabel('log N females');

%Make survival and fecundity plots
subplot(2,2,3), plot(age_vec,surv);
axis square;
axis([0.5 age_max+0.5 0 1]);
xlabel('age in years');
ylabel('survival probability');
subplot(2,2,4), plot(age_vec,fert);
axis square;
axis([0.5 age_max+0.5 0 1]);
xlabel('age in years');
ylabel('m (fertility)');


% *******************************************
% Density-dependence in survival & fertility
% Set male kill rates
% *******************************************

format long g

% ***********************************************************************************************
%Initial population size vector
Nd=4000;

% Set time limit for projection
tlimit = 100;

% Set kill range
min_kill = 0; % Set minimum male (>5 year) kill rate
max_kill = 50; % Set maximum male (>5 year) kill rate

% Minimum operative reproductive sex ratio
orsr_min = 0.1; % Proportion of population that must be mature males to allow females to breed
% ************************************************************************************************

karray = linspace(min_kill,max_kill,100); % male kill rate array

min_popd_vec = zeros(length(karray),1); % Minimum population size vector
min_mal6p_vec = zeros(length(karray),1); % Minimum male (>5 years) vector
lambdad_vec = zeros(length(karray),1); % mean population rate of change (r)

for q=1:length(karray);

% define start survival
        % survival parameters with Nmax = 4000
        % s0
        a_s0=0.913;
        b_s0=6.284615518;
        x0_s0=4224.122758;

        % s1+
        a_s1p=0.995;
        b_s1p=3.9151;
        x0_s1p=7286.656;

        % Redefine survival probabilities
        s0_d = (a_s0/(1+((sum(Nd))/x0_s0)^b_s0));
            if s0_d < 0.55;
                s0_d = 0.55;
            else
                s0_d = s0_d;
            end
            
        s1p_d = (a_s1p/(1+((sum(Nd))/x0_s1p)^b_s1p));
            if s1p_d < 0.92;
               s1p_d = 0.92;
            else
                s1p_d = s1p_d;
            end
        
        % fecundity parameters with Nmax = 4000
        a_ma = 0.7013;
        b_ma = 8.0804;
        x0_ma = 3959.0825;
        a_md = (a_ma/(1+((sum(Nd))/x0_ma)^b_ma));
            if a_md < 0.4004;
                a_md = 0.4004;
            else
                a_md = a_md;
            end
            
        % Re-define fertility vector
        b_fert = 0.4453;
        x0_fert = 6.3898;
        for j=1:age_max;
            md=a_md*exp(-0.5*(log(j/x0_fert)/b_fert)^2);
            m_vec(j) = md;
        end

% sex ratios
sr = 0.5; % overall sex ratio
x = 0.5; % foetal sex ratio (Choquenot 1993)

%total females in population
f=Nd*sr;

%total males in population
mal=Nd*(1-sr);

% the normal matrix
        ad=[0 s0_d*m_vec(1) s0_d*m_vec(2) s0_d*m_vec(3) s0_d*m_vec(4) s0_d*m_vec(5) s0_d*m_vec(6) s0_d*m_vec(7) s0_d*m_vec(8) s0_d*m_vec(9) s0_d*m_vec(10) s0_d*m_vec(11) s0_d*m_vec(12) s0_d*m_vec(13) s0_d*m_vec(14) s0_d*m_vec(15) s0_d*m_vec(16) 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
            s1p_d 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
            0 s1p_d 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
            0 0 s1p_d 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
            0 0 0 s1p_d 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
            0 0 0 0 s1p_d 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 s1p_d 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 s1p_d 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 s1p_d 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 s1p_d 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 s1p_d 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 0 s1p_d 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 0 0 s1p_d 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 0 0 0 s1p_d 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 0 0 0 0 s1p_d 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 s1p_d 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 s1p_d 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
            0 s0_d*m_vec(1) s0_d*m_vec(2) s0_d*m_vec(3) s0_d*m_vec(4) s0_d*m_vec(5) s0_d*m_vec(6) s0_d*m_vec(7) s0_d*m_vec(8) s0_d*m_vec(9) s0_d*m_vec(10) s0_d*m_vec(11) s0_d*m_vec(12) s0_d*m_vec(13) s0_d*m_vec(14) s0_d*m_vec(15) s0_d*m_vec(16) 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 s1p_d 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 s1p_d 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 s1p_d 0 0 0 0 0 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 s1p_d 0 0 0 0 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 s1p_d 0 0 0 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 s1p_d 0 0 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 s1p_d 0 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 s1p_d 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 s1p_d 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 s1p_d 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 s1p_d 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 s1p_d 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 s1p_d 0 0 0 0
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 s1p_d 0 0 0
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 s1p_d 0 0
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 s1p_d 0];

% eigenvalues & eigenvectors
[wd,dd] = eig(ad);
vd=conj(inv(wd));

%lambda
lambdad=diag(dd);

% max lambdas
[maxlambdad,imax]=max(diag(dd));

%logmaxlambda=exponential rate of increase (r)
rd=log(maxlambdad);

% stable age distribution
wd=wd(:,imax);

% stable size distribution (scaled to sum to 1)
ssdd=wd(:,1)/sum(wd(:,1));

% size of matrix
k=size(a,1);

% ssd classes
    % female
    ssdd_juvf = sum(ssdd(1:primi-1,1));
    ssdd_adf = sum(ssdd(primi:age_max,1));

    % male
    ssdd_juvm = sum(ssdd(age_max+1:age_max+primi-1,1));
    ssdd_adm = sum(ssdd(age_max+primi:k,1));
    
    % males 6+
    ssdd_mal6p = sum(ssdd(23:34));
    
% pick initial vectors
nd=ones(k,1);
nd=ssdd*Nd;

%set population size year step vector
popd_vec=ones(tlimit+1,1);
popd_vec(1,1)=sum(nd);

%set female population size year step vector
femd_vec=ones(tlimit+1,1);
femd_vec(1,1)=sum(nd(1:17));

%set male population size year step vector
mald_vec=ones(tlimit+1,1);
mald_vec(1,1)=sum(nd(18:34));

% set large male year step vector
male6p_vec = ones(tlimit+1,1);
male6p_vec(1,1)=sum(nd(23:34));

%set year step vector
yr_vec=ones(tlimit+1,1);
for c=1:tlimit
    yr_vec(c+1,1)=c;
    yr_vec(1,1)=0;
end % c loop

%then iterate
for i=1:tlimit;

    nd=ad*nd;
    
        % Set negative density feedback function for the matrix

        % survival parameters with Nmax = 4000
        % s0
        a_s0=0.913;
        b_s0=6.284615518;
        x0_s0=4224.122758;

        % s1+
        a_s1p=0.995;
        b_s1p=3.9151;
        x0_s1p=7286.656;

        % Redefine survival probabilities
        s0_d = (a_s0/(1+((sum(nd))/x0_s0)^b_s0));
            if s0_d < 0.55;
                s0_d = 0.55;
            else
                s0_d = s0_d;
            end
            
        s1p_d = (a_s1p/(1+((sum(nd))/x0_s1p)^b_s1p));
            if s1p_d < 0.92;
               s1p_d = 0.92;
            else
                s1p_d = s1p_d;
            end
        
        % fecundity parameters with Nmax = 4000
        a_ma = 0.7013;
        b_ma = 8.0804;
        x0_ma = 3959.0825;
        a_md = (a_ma/(1+((sum(nd))/x0_ma)^b_ma));
            if a_md < 0.4004;
                a_md = 0.4004;
            else
                a_md = a_md;
            end
            
        % Re-define fertility vector
        for j=1:age_max;
            md=a_md*exp(-0.5*(log(j/x0_fert)/b_fert)^2);
            m_vec(j) = md;
        end

        % Minimum operative reproductive sex ratio
        if sum(nd(18:34)) / sum(nd(1:17)) < orsr_min;
            %m_vec = m_vec * (sum(nd(18:34)) / sum(nd(1:17)));
            m_vec = 0;
        else
            m_vec = m_vec;
        end
                
        % *****************************************************************
        % Kill large males
        % *****************************************************************
        
        male_kill = karray(q); % set annual number of large males to kill
        %male_kill = 20;
        nm_kill = male_kill/12; % uniform distribution of kill in ages > 5
        kill_vec = zeros(34,1);
        kill_vec(23:34) = nm_kill;
            if sum(nd) < male_kill;
                nd = nd;
            else
                nd = nd - nm_kill; % Adjust nd vector for kill
            end
        
        % *****************************************************************
        
        % the new matrix
        ad=[0 s0_d*m_vec(1) s0_d*m_vec(2) s0_d*m_vec(3) s0_d*m_vec(4) s0_d*m_vec(5) s0_d*m_vec(6) s0_d*m_vec(7) s0_d*m_vec(8) s0_d*m_vec(9) s0_d*m_vec(10) s0_d*m_vec(11) s0_d*m_vec(12) s0_d*m_vec(13) s0_d*m_vec(14) s0_d*m_vec(15) s0_d*m_vec(16) 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
            s1p_d 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
            0 s1p_d 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
            0 0 s1p_d 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
            0 0 0 s1p_d 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
            0 0 0 0 s1p_d 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 s1p_d 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 s1p_d 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 s1p_d 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 s1p_d 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 s1p_d 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 0 s1p_d 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 0 0 s1p_d 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 0 0 0 s1p_d 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 0 0 0 0 s1p_d 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 s1p_d 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 s1p_d 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
            0 s0_d*m_vec(1) s0_d*m_vec(2) s0_d*m_vec(3) s0_d*m_vec(4) s0_d*m_vec(5) s0_d*m_vec(6) s0_d*m_vec(7) s0_d*m_vec(8) s0_d*m_vec(9) s0_d*m_vec(10) s0_d*m_vec(11) s0_d*m_vec(12) s0_d*m_vec(13) s0_d*m_vec(14) s0_d*m_vec(15) s0_d*m_vec(16) 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 s1p_d 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 s1p_d 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 s1p_d 0 0 0 0 0 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 s1p_d 0 0 0 0 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 s1p_d 0 0 0 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 s1p_d 0 0 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 s1p_d 0 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 s1p_d 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 s1p_d 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 s1p_d 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 s1p_d 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 s1p_d 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 s1p_d 0 0 0 0
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 s1p_d 0 0 0
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 s1p_d 0 0
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 s1p_d 0];

     
     popd_vec(i+1,1)=(sum(nd));    
     femd_vec(i+1,1)=sum(nd(1:17));
     mald_vec(i+1,1)=sum(nd(18:34));
     male6p_vec(i+1,1)=sum(nd(23:34));
     
end % i loop

log_popd_vec=log10(popd_vec);

%total population size after 'tlimit' years
%popd_st=Nd
%popd_end=sum(nd)
%tlimit
%maxlambdad
%rd
%format short
%ssdd_juv = ssdd_juvf + ssdd_juvm
%ssdd_ad = ssdd_adf + ssdd_adm

%Make density dependent plots
subplot(2,2,1), plot(yr_vec,popd_vec);
axis square;
axis([-0.5 tlimit+1 0 (max(popd_vec)+(0.25*(max(popd_vec))))]);
xlabel('year');
ylabel('N');
subplot(2,2,2), plot(yr_vec,log_popd_vec);
axis square;
axis([-0.5 tlimit+1 0 (max(log_popd_vec)+(0.25*(max(log_popd_vec))))]);
xlabel('year');
ylabel('log N');
subplot(2,2,3), plot(yr_vec,male6p_vec);
axis square;
axis([-0.5 tlimit+1 0 (max(male6p_vec)+(0.25*(max(male6p_vec))))]);
xlabel('year');
ylabel('N bulls >5 years');

% Place data in q vectors
min_popd_vec(q) = min(popd_vec);
minmal6p_vec(q) = min(male6p_vec);

min_male6p_sub = find(male6p_vec == min(male6p_vec));
min_mald = mald_vec(min_male6p_sub);

minpmal6p_vec(q) = min(male6p_vec) / min_mald;

% Calculate stochastic r
popd1_vec=zeros(length(popd_vec)+1,1);
popd1_vec(2:length(popd_vec)+1,1)=popd_vec;
lambda_m = mean(popd_vec(2:length(popd_vec),1) ./ popd1_vec(2:length(popd_vec),1));
lambdad_vec(q) = lambda_m;

end % q loop

% Fix min_popd_vec
reals1 = find(min_popd_vec>0);
reals1_sub = find(diff(reals1)>1);
    if sum(reals1_sub) == 0;
        reals1_sub = length(min_popd_vec);
    else
        reals1_sub = reals1_sub;
    end

    if reals1_sub(1) ~= length(min_popd_vec);
            min_popd_vec((reals1_sub(1)+1):length(min_popd_vec)) = 0;
        else
            min_popd_vec = min_popd_vec;
    end

% Fix minmal6p
reals2 = find(minmal6p_vec>0);
reals2_sub = find(diff(reals2)>1);
    if sum(reals2_sub) == 0;
        reals2_sub = length(minmal6p_vec);
    else
        reals2_sub = reals2_sub;
    end

    if reals2_sub(1) ~= length(minmal6p_vec);
            minmal6p_vec((reals2_sub(1)+1):length(minmal6p_vec)) = 0;
        else
            minmal6p_vec = minmal6p_vec;
    end

% Fix minpmal6p
reals3 = find(minpmal6p_vec>0);
reals3_sub = find(diff(reals3)>1);
    if sum(reals3_sub) == 0;
        reals3_sub = length(minpmal6p_vec);
    else
        reals3_sub = reals3_sub;
    end

    if reals3_sub(1) ~= length(minpmal6p_vec);
            minpmal6p_vec((reals3_sub(1)+1):length(minpmal6p_vec)) = 0;
        else
            minpmal6p_vec = minpmal6p_vec;
    end

% Fix stochastic r for real numbers only
reals4 = find(lambdad_vec>0.7);
reals4_sub = find(diff(reals4)>1);
    if sum(reals4_sub) == 0;
        reals4_sub = length(lambdad_vec);
    else
        reals4_sub = reals4_sub;
    end

    if reals4_sub(1) ~= length(lambdad_vec);
            lambdad_vec((reals4_sub(1)+1):length(lambdad_vec)) = 0;
        else
            lambdad_vec = lambdad_vec;
    end
        
% Make plots
subplot(2,2,1), plot(karray,minmal6p_vec);
axis square;
axis([-0.5 max(karray) 0 (max(minmal6p_vec)+(0.25*(max(minmal6p_vec))))]);
xlabel('Mature males killed per year');
ylabel('Min mature males');
subplot(2,2,2), plot(karray,minpmal6p_vec);
axis square;
axis([-0.5 max(karray) 0 (max(minpmal6p_vec)+(0.25*(max(minpmal6p_vec))))]);
xlabel('Mature males killed per year');
ylabel('Min prop mature males');
subplot(2,2,3), plot(karray,min_popd_vec);
axis square;
axis([-0.5 max(karray) 0 (max(min_popd_vec)+(0.25*(max(min_popd_vec))))]);
xlabel('Mature males killed per year');
ylabel('Min N');
subplot(2,2,4), plot(karray,lambdad_vec);
axis square;
axis([-0.5 max(karray) 0.95 1.05]);
xlabel('Mature males killed per year');
ylabel('Mean lambda');
hold
plot(karray,1,'r-');
hold

% *****************************************
% Allometrically derived vital rates
% *****************************************

% Bos taurus - Barry's database:
% based on 725 kg mass
% primariparity = 18 months
% fecundity = 1 calf/year
% longevity = 12 years
% average density = 5 km-2

% Allometry equations
% Average female weight: 200-300 kg
% Average male weight: 500-600 kg

mass_avg_f_cob = 250;
mass_avg_m_cob = 550;

mass_avg_f_lit = 350;
mass_avg_m_lit = 700; % Wootton (1987); Nowak & Paradiso (1983)

% Maximum rate of population increase (Hennemann 1983; Caughley & Krebs 1983)
rmax1 = 18 * (mass_avg_f_cob*1000)^(-0.36); % Caughley & Krebs (1983)
rmax2 = 4.9 * (mass_avg_f_cob*1000)^(-0.2622); % Hennemann (1983); go with this one

lamb_max1 = exp(rmax1);
lamb_max2 = exp(rmax2);

% Average population density (Damuths 1981; Freeman 1990; Damuths 1987: Biol. J. Linn. Soc. 31:193), per km2
log10_pop_d1 = 4.23 - (0.75 * log10(mean([mass_avg_f_cob mass_avg_m_cob])*1000)); % Damuths (1981)
pop_d1 = 10^log10_pop_d1

log10_pop_d2 = 4.196 - (0.74 * log10(mean([mass_avg_f_cob mass_avg_m_cob])*1000)); % Freeland (1990)
pop_d2 = 10^log10_pop_d2

% For tropical herbivores
pop_d3 = 12000 * (mass_avg_f_cob*1000)^-0.73 % Damuths (1987) ****Corrected, but check

% Home range size (Harestad & Bennell 1979: Ecol. 60:389)
hr1 = 0.002 * (mass_avg_f_cob*1000)^1.02; % area in ha

% Maximum lifespan (Weston 1979)
lspd1 = 2140 * (mass_avg_f_cob)^0.22; % in days
lspy1 = lspd1/365

% Gestation time (Weston 1979)
gestd1 = 117.5 * (mass_avg_f_cob)^0.16; % in days
gestm1 = gestd1/30

% Birth rate (Millar 1981)
log_br1 = 0.888 * (log(mass_avg_f_cob)) - 3.097; 
br1 = exp(log_br1)

% Litter size (Millar 1981: Evolution 35:1149; Sacher & Staffeldt 1974: Am. Nat. 108:593)
ls = 3.034 - (0.298 * log(mass_avg_f_cob)) % Millar (1981)

% Max m
m_max = ls/2; % if all animals survived; per sex

% Age to maturity (Weston 1979: Afr. J. Ecol. 17:185; Bluewiess 1978: Oecologia 37:257)
matd1 = 257 * (mass_avg_f_cob)^0.24; % McDonald (1984)
maty1 = matd1/365

% Individual growth rate (Case 1978: Quart. Rev. Biol. 53:243)
log10_grf_d1 = (0.52 * (log10(mass_avg_f_cob*1000))) - 0.25; % Case (1978)
gr_f_d1 = 10^log10_grf_d1 % g/day

log10_grm_d1 = (0.52 * (log10(mass_avg_m_cob*1000))) - 0.25; % Case (1978)
gr_m_d1 = 10^log10_grm_d1 % g/day

