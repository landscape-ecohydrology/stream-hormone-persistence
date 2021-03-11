% % ----------------------------------------------------------------------- 
% % ----------------------------------------------------------------------- 
% % Hormone Compounds in a Riverine Environment
% % Frederick Cheng
% % 
% % This code simulates the contaminant transport of hormones with a
% % stream and its hyporheic zone. The transient behaviour is developed
% % using the finite volume and Euler's methods
% % ----------------------------------------------------------------------- 


function HormoneADRE() 
clc; clear;  close all; 
clc, format compact
% Additions:    E2 and E1 species
%               E2 to E1 in MC
%               Exchange between MC and TS
%               Retardation factor
%               Permanent decay

 
% -------------------------------------------------------------------------
% Model parameters (User Defined)------------------------------------------
N        = 500;              % number of control volumes
NCyc     = 1;                % number of cycles
duration = 365;              % duration (days)
dt       = 0.01;             % time step (days)
L        = 100;              % domain length [km]
beta     = 0.5;              % spatial weighting (0.5 for central weighting)
theta    = 0.5;              % temporal weighting (crank-nicolson for 0.5)
dx       = L/N;              % length of control volume

% Physical parameters (User Defined)---------------------------------------
C_in_2 = xlsread('Network Data.xlsx','Velocities','G3:H682'); % ug/m3 
C_in_1 = xlsread('Network Data.xlsx','Velocities','E3:F682'); % ug/m3

temp         = interpC(dt, duration/NCyc, C_in_2);  % initial concentration of E2 at x=0 (type I) (mg/m3)
C_1_bc_orig  = interpC(dt, duration/NCyc, C_in_1);  % initial concentration of E1 at x=0 (type I) (mg/m3)

% QA/QC
C_2_bc = temp;
C_1_bc = C_1_bc_orig;
C_2_bc(isnan(C_2_bc)) = 0;
C_1_bc(isnan(C_1_bc)) = 0;
C_1_bc(abs(C_1_bc - 0.1148)<=0.001) = 0;

v       = [35 20];         % stream velocity [km/d] 
D       = 8.6;             % diffusion coefficient [m2/s]
ex      = [1.5 1.5];       % exchange rate between MC and TS [1/d]

R       = [1.4 1.4];       % ratio of MC/TS areas
pb      = 1.3;             % bulk density of TS sediment [g/L]
n       = 0.4;             % porosity [-]
kd_2    = 75;              % partition coefficient for E2 [L/g] 
kd_1    = 65;              % partition coefficient for E1 [L/g]
Css     = 2.6e-4;          % suspended solids concentration [g/L]
k_sett  = 180;             % settling rate constant [1/d]


k_2P_MC = 2.25;            % permanent decay in MC [1/d] 
k_21_MC = 0.25;            % decay of E2 to E1 in MC [1/d]
k_2P_TS = 0;               % permanent decay in MC [1/d]
k_21_TS = [0.9 0.9];       % decay of E2 to E1 in TS [1/d]
k_12_TS = [0.06 0.06]*0;   % decay of E2 to E1 in TS [1/d]

k_1P_MC = 12;              % permanent decay in MC [1/d]
k_1P_TS = 0.4;             % permanent decay in MC [1/d]


HiLo = zeros(duration/dt+2,1);
HiLo(1:180*(1/dt),1) = 1;      % 1 = high; 2 = low
HiLo(180*(1/dt)+1:288*(1/dt),1) = 2;
HiLo(288*(1/dt)+1:365*(1/dt)+2,1) = 1;

C_2_bc = C_2_bc * (1+kd_2*Css);
C_1_bc = C_1_bc * (1+kd_1*Css);


% Calculate statistics for numerical stability
fprintf('\nGrid Peclet Number:     %8.2f', v*dx/D)   
fprintf('\nPeclet Number:          %8.2f', v*L/D)  
fprintf('\nCourant Number:         %8.2f', v*dt/dx)  
fprintf('\nLet''s go\n\n')

% Initialize matrices -----------------------------------------------------
C       = zeros(4*N,1);     % unknown concntrations (column vector; first quarter
                            % for E2 in MC; second quarter for E2 in TS; third for
                            % E1 in MC; fourth for E1 in TS)
                            
A_E2_MC = zeros(N,N);       % coefficient matrix for E2 in MC
A_E2_TS = zeros(N,N);       % coefficient matrix for E2 in TS
A_E1_MC = zeros(N,N);       % coefficient matrix for E1 in MC
A_E1_TS = zeros(N,N);       % coefficient matrix for E1 in TS

B_E2_MC = zeros(N,1);       % right hand side (columns vector) for E2 in MC
B_E2_TS = zeros(N,1);       % right hand side (columns vector) for E2 in TS
B_E1_MC = zeros(N,1);       % right hand side (columns vector) for E1 in MC
B_E1_TS = zeros(N,1);       % right hand side (columns vector) for E1 in TS

timeindex = 0;

for t=0:dt:duration
    
    if mod(t,10) == 0
        disp(t)
    end

    
AQ_MC_2 = 1/(1+kd_2*Css);   
AQ_MC_1 = 1/(1+kd_1*Css);   

SD_MC_2 = 1-1/(1+kd_2*Css);
SD_MC_1 = 1-1/(1+kd_1*Css);

AQ_TS_2 = 1/(n+kd_2*pb);   
AQ_TS_1 = 1/(n+kd_1*pb);   

% SD_TS_2 = 1-1/(1+kd_2*Css);
% SD_TS_1 = 1-1/(1+kd_1*Css);

% -------------------------------------------------------------------------
%Build the system of equations LHS matrix ---------------------------------

% FOR E2 in MAIN CHANNEL---------------------------------------------------
a_E2_MC =     D/dx/dx + v(HiLo(timeindex+1))*beta      /dx;
b_E2_MC = - 2*D/dx/dx + v(HiLo(timeindex+1))*(1-2*beta)/dx - ex(HiLo(timeindex+1))*AQ_MC_2 - k_2P_MC*AQ_MC_2 - k_21_MC*AQ_MC_2 - k_sett*SD_MC_2;
c_E2_MC =     D/dx/dx - v(HiLo(timeindex+1))*(1-beta)  /dx; 

a_E2_MC_prime = - a_E2_MC * theta;
b_E2_MC_prime = - b_E2_MC * theta +1/dt;
c_E2_MC_prime = - c_E2_MC * theta;

S_E2_MC_ex = ex(HiLo(timeindex+1))*AQ_TS_2;

MC_E2_bi_ss  = - D/dx/dx - v(HiLo(timeindex+1))/2/dx - ex(HiLo(timeindex+1))*AQ_MC_2 - k_2P_MC*AQ_MC_2 - k_21_MC*AQ_MC_2;

% for U/S BC
MC_E2_b1_s   = -3*D/dx/dx - v(HiLo(timeindex+1))*(beta)  /dx - ex(HiLo(timeindex+1))*AQ_MC_2 - k_2P_MC*AQ_MC_2 - k_21_MC*AQ_MC_2;
MC_E2_c1_s   =    D/dx/dx - v(HiLo(timeindex+1))*(1-beta)/dx;

MC_E2_a1_s   = -2*D/dx/dx - v(HiLo(timeindex+1))/dx;

for i=1:N
    if (i>1) && (i<N) %internal cells
        A_E2_MC(i,i-1) = a_E2_MC_prime;
        A_E2_MC(i,i  ) = b_E2_MC_prime;
        A_E2_MC(i,i+1) = c_E2_MC_prime;
    elseif (i==1) %left boundary - Dirichlet condition - C=C_hat
        A_E2_MC(1,1) = 1;
        
    elseif (i==N) %right boundary - outflow condition -D dC/dx=0
        A_E2_MC(N,N)   = -0.5*MC_E2_bi_ss + 1/dt;
        A_E2_MC(N,N-1) = a_E2_MC_prime;
    end
end  

% FOR E2 in TRANSIENT STORAGE----------------------------------------------
a_E2_TS = 0;
b_E2_TS = (-ex(HiLo(timeindex+1))*R(HiLo(timeindex+1))*AQ_TS_2 - k_21_TS(HiLo(timeindex+1))*AQ_TS_2 - k_2P_TS*AQ_TS_2); 
c_E2_TS = 0; 

a_E2_primeTS = -a_E2_TS * theta;
b_E2_primeTS = -b_E2_TS * theta +1/dt;
c_E2_primeTS = -c_E2_TS * theta;

S_E2_TS_ex = ex(HiLo(timeindex+1))*R(HiLo(timeindex+1))*AQ_MC_2 + k_sett*R(HiLo(timeindex+1))*SD_MC_2; 
S_E2_TS_k12 = k_12_TS(HiLo(timeindex+1))*AQ_TS_1;

TS_E2_bi_ss  = (-ex(HiLo(timeindex+1))*R(HiLo(timeindex+1))*AQ_TS_2 - k_21_TS(HiLo(timeindex+1))*AQ_TS_2 - k_2P_TS*AQ_TS_2);
TS_E2_bi_s   = (-ex(HiLo(timeindex+1))*R(HiLo(timeindex+1))*AQ_TS_2 - k_21_TS(HiLo(timeindex+1))*AQ_TS_2 - k_2P_TS*AQ_TS_2);

for i=1:N
    if (i>1) && (i<N) %internal cells
        A_E2_TS(i,i-1) = a_E2_primeTS;
        A_E2_TS(i,i  ) = b_E2_primeTS;
        A_E2_TS(i,i+1) = c_E2_primeTS;
    elseif (i==1) %left boundary - Dirichlet condition - C=C_hat
        A_E2_TS(i,i  ) = b_E2_primeTS;
        A_E2_TS(i,i+1) = c_E2_primeTS;

    elseif (i==N) %right boundary - outflow condition -D dC/dx=0
        A_E2_TS(N,N)   = -0.5*TS_E2_bi_ss + 1/dt;
        A_E2_TS(N,N-1) = a_E2_primeTS;
    end
end



% FOR E1 in MAIN CHANNEL---------------------------------------------------
a_E1_MC =     D/dx/dx + v(HiLo(timeindex+1))*beta      /dx;
b_E1_MC = - 2*D/dx/dx + v(HiLo(timeindex+1))*(1-2*beta)/dx - ex(HiLo(timeindex+1))*AQ_MC_1 - k_1P_MC*AQ_MC_1 - k_sett*SD_MC_1;
c_E1_MC =     D/dx/dx - v(HiLo(timeindex+1))*(1-beta)  /dx; 

a_E1_MC_prime = - a_E1_MC * theta;
b_E1_MC_prime = - b_E1_MC * theta +1/dt;
c_E1_MC_prime = - c_E1_MC * theta;

S_E1_MC_ex  = ex(HiLo(timeindex+1))*AQ_TS_1;
S_E1_MC_k21 = k_21_MC*AQ_MC_2;

MC_E1_bi_ss  = - D/dx/dx - v(HiLo(timeindex+1))/2/dx - ex(HiLo(timeindex+1))*AQ_MC_1 - k_1P_MC*AQ_MC_1 - k_sett*SD_MC_1;

% for U/S BC
MC_E1_b1_s   = -3*D/dx/dx - v(HiLo(timeindex+1))*(beta)  /dx - ex(HiLo(timeindex+1))*AQ_MC_1 - k_1P_MC*AQ_MC_1 - k_sett*SD_MC_1;
MC_E1_c1_s   =    D/dx/dx - v(HiLo(timeindex+1))*(1-beta)/dx;

MC_E1_a1_s   = -2*D/dx/dx - v(HiLo(timeindex+1))/dx;

for i=1:N
    if (i>1) && (i<N) %internal cells
        A_E1_MC(i,i-1) = a_E1_MC_prime;
        A_E1_MC(i,i  ) = b_E1_MC_prime;
        A_E1_MC(i,i+1) = c_E1_MC_prime;
    elseif (i==1) %left boundary - Dirichlet condition - C=C_hat
        A_E1_MC(1,1) = 1;
              
    elseif (i==N) %right boundary - outflow condition -D dC/dx=0
        A_E1_MC(N,N)   = -0.5*MC_E1_bi_ss + 1/dt;
        A_E1_MC(N,N-1) = a_E1_MC_prime;
    end
end  

% FOR E1 in TRANSIENT STORAGE----------------------------------------------
a_E1_TS = 0;
b_E1_TS = (-ex(HiLo(timeindex+1))*R(HiLo(timeindex+1))*AQ_TS_1 - k_12_TS(HiLo(timeindex+1))*AQ_TS_1 - k_1P_TS*AQ_TS_1);
c_E1_TS = 0; 

a_E1_primeTS = -a_E1_TS * theta;
b_E1_primeTS = -b_E1_TS * theta +1/dt;
c_E1_primeTS = -c_E1_TS * theta;

S_E1_TS_ex  = ex(HiLo(timeindex+1))*R(HiLo(timeindex+1))*AQ_MC_1 + k_sett*R(HiLo(timeindex+1))*SD_MC_1;
S_E1_TS_k21 = k_21_TS(HiLo(timeindex+1))*AQ_TS_2;

TS_E1_bi_ss  = (-ex(HiLo(timeindex+1))*R(HiLo(timeindex+1)) - k_12_TS(HiLo(timeindex+1))- k_1P_TS)*AQ_TS_1;  % downstream boundary condition
TS_E1_bi_s   = (-ex(HiLo(timeindex+1))*R(HiLo(timeindex+1)) - k_12_TS(HiLo(timeindex+1))- k_1P_TS)*AQ_TS_1;

for i=1:N
    if (i>1) && (i<N) %internal cells
        A_E1_TS(i,i-1) = a_E1_primeTS;
        A_E1_TS(i,i  ) = b_E1_primeTS;
        A_E1_TS(i,i+1) = c_E1_primeTS;
    elseif (i==1) %left boundary - Dirichlet condition - C=C_hat
        A_E1_TS(i,i  ) = b_E1_primeTS;
        A_E1_TS(i,i+1) = c_E1_primeTS;

    elseif (i==N) %right boundary - outflow condition -D dC/dx=0
        A_E1_TS(N,N)   = -0.5*TS_E1_bi_ss + 1/dt;
        A_E1_TS(N,N-1) = a_E1_primeTS;
    end
end


% -------------------------------------------------------------------------
% Build the system of equations RHS vector B ------------------------------

% FOR E2 in MAIN CHANNEL-----------------------------------------------
for i=1:N
    if (i==1)           % upstream B.C. (type I)
         B_E2_MC(i)   = C_2_bc(timeindex+1);

    elseif (i==N)       % downstream B.C. (type II)
         B_E2_MC(i)   = (0.5 * a_E2_MC    * C(i-1) ...
                             + MC_E2_bi_ss* C(i) ...
                             +  C(i)/dt)...
                             + -ex(HiLo(timeindex+1))*AQ_MC_2        * C(i+N);
    else                %-All internal CVs
         B_E2_MC(i)   = (1-theta) * (a_E2_MC * C(i-1) ...
                                  +  b_E2_MC * C(i) ... 
                                  +  c_E2_MC * C(i+1) ...
                                  +  S_E2_MC_ex * C(i+N)) ... 
                                  +  1/dt    * C(i);
    end
end

% FOR E2 in TRANSIENT STORAGE----------------------------------------------
    for i=1:N
        if     (i==1)       %-left hand B.C. (type I)
            B_E2_TS(i)    =   0.5 * b_E2_TS           * C(i+N) ...
                                + 1/dt                * C(i+N) ...
                                + 0.5 * (-ex(HiLo(timeindex+1))*R(HiLo(timeindex+1)))*AQ_MC_2  * C(i);
        elseif (i==N)       %-right hand B.C. (type II)
            B_E2_TS(i)    = 0.5*(-ex(HiLo(timeindex+1))*R(HiLo(timeindex+1))*AQ_TS_2 * C(i+N) + C(i+N)/dt ...  
                                - ex(HiLo(timeindex+1))*R(HiLo(timeindex+1))*AQ_MC_2 * C(i));   
        else               %-All internal CVs
            B_E2_TS(i)    = (1-theta) * b_E2_TS       * C(i+N) ...
                                      + 1/dt          * C(i+N) ...
                                      + (1-theta) * S_E2_TS_ex * C(i)...
                                      + (1-theta) * S_E2_TS_k12 * C(i+3*N); 
        end
    end
    
% FOR E1 in MAIN CHANNEL-----------------------------------------------
for i=1:N
    if (i==1)           % upstream B.C. (type I)
         B_E1_MC(i)   = C_1_bc(timeindex+1); 

    elseif (i==N)       % downstream B.C. (type II)
         B_E1_MC(i)   = (0.5 * a_E1_MC    * C(i-1+2*N) ...
                             + MC_E1_bi_ss* C(i+2*N) ...
                             +  C(i+2*N)/dt)...
                             + -ex(HiLo(timeindex+1)) *AQ_TS_1       * C(i+3*N)...
                             + -k_21_MC *AQ_MC_2  * C(i);
    else                %-All internal CVs
         B_E1_MC(i)   = (1-theta) * (a_E1_MC * C(i-1+2*N) ...
                                  +  b_E1_MC * C(i+2*N) ... 
                                  +  c_E1_MC * C(i+1+2*N) ...
                                  +  S_E1_MC_ex * C(i+3*N)...
                                  +  S_E1_MC_k21* C(i)) ... 
                                  +  1/dt    * C(i+2*N);
    end
end

% FOR E1 in TRANSIENT STORAGE----------------------------------------------
    for i=1:N
        if     (i==1)       %-left hand B.C. (type I)
            B_E1_TS(i)    =   0.5 * b_E1_TS           * C(i+3*N) ...
                                + 1/dt                * C(i+3*N) ...
                                + 0.5 * (-ex(HiLo(timeindex+1))*R(HiLo(timeindex+1)))*AQ_MC_1 * C(i+2*N); 
        elseif (i==N)       %-right hand B.C. (type II)
            B_E1_TS(i)    = 0.5*(-ex(HiLo(timeindex+1))*R(HiLo(timeindex+1))*AQ_TS_1 * C(i+3*N) + C(i+3*N)/dt ...    
                                - ex(HiLo(timeindex+1))*R(HiLo(timeindex+1))*AQ_MC_1 * C(i+2*N));   
        else               %-All internal CVs
            B_E1_TS(i)    = (1-theta) * b_E1_TS       * C(i+3*N) ...
                                      + 1/dt          * C(i+3*N) ...
                                      + (1-theta) * S_E1_TS_ex * C(i+2*N)...
                                      + (1-theta) * S_E1_TS_k21* C(i+N); 
        end
    end
    
    
    
% APPEND MATRICES TOGETHER-------------------------------------------------
    % Creating diagonal exchange matrix
    BLANK = zeros(N,N);  
    
    % E2 in Main Channel
    I_E2_MC = -S_E2_MC_ex*theta * eye(N);

    % E2 in Transient Storage
    I_E2_TS_ex = -S_E2_TS_ex*theta * eye(N);
    I_E2_TS_k12 = -S_E2_TS_k12 *theta * eye(N);

    % E1 in Main Channel
    I_E1_MC     = -S_E1_MC_ex  *theta * eye(N);
    I_E2_MC_k21 = -S_E1_MC_k21 *theta * eye(N);

    % E1 in Transient Storage
    I_E1_TS_ex  = -S_E1_TS_ex  *theta * eye(N);
    I_E2_TS_k21 = -S_E1_TS_k21 *theta * eye(N);


    %Append matrices together
    A_new = [A_E2_MC     I_E2_MC     BLANK       BLANK;
             I_E2_TS_ex  A_E2_TS     BLANK       I_E2_TS_k12;
             I_E2_MC_k21 BLANK       A_E1_MC     I_E1_MC;
             BLANK       I_E2_TS_k21 I_E1_TS_ex  A_E1_TS];
         
    B_new = [B_E2_MC; B_E2_TS; B_E1_MC; B_E1_TS];

    C = A_new \ B_new; 
    timeindex = timeindex +1;
    
end % end time loop


end



function C_2_bc= interpC(dt, duration, C)

time = linspace(dt, duration, duration/dt);
time = [0 time]';
time(:,2) = -9999;


[numrow_C, ~] = size(C);
[numrow_t, ~] = size(time);

for k = 1:numrow_C
    index = find(time >= C(k,1),1);
    time(index,2) = C(k,2);
    
end

for k = 1:numrow_t
   if time(k,2) == -9999      
  time(k,2) =  interp1(C(:,1),C(:,2),time(k,1));
   end   
end

C_2_bc = time(:,2);
end