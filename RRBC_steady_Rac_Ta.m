%clc;
clear, clear global
global sol options

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rotating Rayleigh Bénard Convection 
% Compute critical Rayleigh number and wavenumber at the onset of steady 
% convection looping through a range of Taylor numbers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set boundary conditions
% 1: FT, Stress-free
% 2: FF, Stress-free
% 3: FT, No-slip
% 4: FF, No-slip
% 5: FT, mixed velocity 
% 6: FF, mixed velocity
% 7: mixed thermal and mixed velocity (FT,SF on same boundary)
BCs = 1; 

% Set fixed parameters
L = 2; % domain length
Pr = 1; % Prandtl number
S=0; % growth rate; set to zero for marginal solutions
fprintf('L = %9.2f, Pr = %9.2f, s=%9.5f\n', L, Pr, S);

% define domain mesh
zmesh =  linspace(0,L,50);

% Define a logarithmically spaced interval for Ta
Ta = logspace(1, 14, 200);   % adjust start/end as needed


% Set initial 'guess' of Rac and kc:
R=43.0268; k= 1.14879;   % solutions for RRBC with Ta=1 and L=2


% Initialise the problem
sol = bvpinit(zmesh,@(z) rbc_init(z,BCs,L), R);

% Loop through Ta 
for i = 1:length(Ta)
    
    sqTa = sqrt(Ta(i));

    % Extrapolate the results already calculated to determine a better
    % 'guess' for this value of Ta
    if i==3
        % Extrapolate linearly 
        k = extrapol(struc.kc(i-2:i-1),log10(Ta(i-2:i-1)),log10(Ta(i)),false);
        R = extrapol(struc.Rac(i-2:i-1),log10(Ta(i-2:i-1)),log10(Ta(i)),false);
        sol.parameters = R;
        fprintf('#L %3i-- Ta=%9.3e: Rae=%9.3e at ke=%8.3f \n',i,Ta(i),R,k);
    end        
    if i>=4
        % Extrapolate quadratically
        k = extrapol(struc.kc(i-3:i-1),log10(Ta(i-3:i-1)),log10(Ta(i)),true);
        R = extrapol(struc.Rac(i-3:i),log10(Ta(i-3:i-1)),log10(Ta(i)),true);
        sol.parameters = R;
        fprintf('#Q %3i-- Ta=%9.3e: Rae=%9.3e at ke=%8.3f \n',i,Ta(i),R,k);
    end

    % Determine Rac and kc 3 times, each time with a smaller tolerance for
    % higher accuracy of the solution
    options = bvpset('RelTol',1e-2, 'Nmax', 10000);
    [kc,Rac] = comp_Rac(k,sqTa,Pr,S,BCs,0.1);
    k = kc; sol.parameters = Rac;
    options = bvpset('RelTol',1e-6, 'Nmax', 10000);
    [kc,Rac] = comp_Rac(k,sqTa,Pr,S,BCs,0.001);
    k = kc; sol.parameters = Rac;
    options = bvpset('RelTol',1e-9, 'Nmax', 10000);
    [kc,Rac] = comp_Rac(k,sqTa,Pr,S,BCs,0.0001);
    fprintf('#1 %3i-- Ta=%9.3e: Rac=%9.3e at kc=%8.3f \n',i,Ta(i),Rac,kc);

    % save data in a structure
    struc.Ta(i) = Ta(i);
    struc.Rac(i) = Rac;
    struc.kc(i) = kc;
    struc.solution(i) = sol;
end


% Save the structure 
savename = "RRB_SF_FT.mat";
save(savename,'struc')

% Plot of Rac vs Ta
figure(1)
loglog(Ta, Rac, 'x-', 'Linewidth', 2); hold off
xlabel('$Ta$', 'Interpreter','latex','Fontsize', 16); 
ylabel('$Ra^c$', 'Interpreter','latex','Fontsize', 16)

% Plot of kc vs Ta
figure(2)
loglog(Ta, kc, 'x-', 'Linewidth', 2); hold off
xlabel('$Ta$', 'Interpreter','latex','Fontsize', 16); 
ylabel('$k^c$', 'Interpreter','latex','Fontsize', 16)

%% ==================================================================
% Compute Rac and kc by finding the minimum of the MSC
function [kc,Rac]=comp_Rac(k,sqTa,Pr,S,BCs,dk)
    TolX=1.e-8; % set the tolerance
    tst_bnd=true;
    lbound=k*(1-dk);ubound=k*(1+dk); % set upper and lower boundaries of k
    while tst_bnd == true
        tst_bnd=false;

        % determine Rac and kc
        [kc, Rac]=fminbnd(@(k) comp_Rak(k,sqTa,Pr,S,BCs),lbound,ubound,optimset('TolX',TolX));
    
        % if Rac is found at either boundary of k then extend the range of
        % k and calculate Rac and kc again
        if (kc-lbound<2*TolX)
            fprintf('Extremum at lower bound: %10.5f, %10.5f, %10.5f\n',...
               lbound,kc,ubound)
            ubound = lbound; lbound =lbound-2*kc*dk;
            tst_bnd=true;
        elseif (ubound-kc<2*TolX)
            fprintf('Extremum at upper bound: %10.5f, %10.5f, %10.5f\n',...
                lbound,kc,ubound)
            lbound = ubound; ubound =ubound+2*kc*dk;
            tst_bnd=true;        
        end
    end
end
% 
%% ==================================================================
% Compute Ram for a given k
function Rak=comp_Rak(k,sqTa,Pr,S,BCs)
global sol options
k2=k*k;
sol = bvp5c(@(z,y,R)rbc_ode(z,y,R,sqTa,k2,Pr,S),@(ya,yb,p)rbc_bc(ya,yb,p,k,BCs),sol,options);
Rak = sol.parameters;
end

%% ==================================================================
% Extrapolate the curve to obtain a more accurate guess for the next Rac
% and kc
function F=extrapol(Fi,Xi,X,quadratic)
if quadratic
    F = Fi(1)*(X-Xi(2))/(Xi(1)-Xi(2))*(X-Xi(3))/(Xi(1)-Xi(3))+Fi(2)*(X-Xi(1))/(Xi(2)-Xi(1))*(X-Xi(3))/(Xi(2)-Xi(3))+Fi(3)*(X-Xi(1))/(Xi(3)-Xi(1))*(X-Xi(2))/(Xi(3)-Xi(2));
else
    F = Fi(1)*(X-Xi(2))/(Xi(1)-Xi(2)) +Fi(2)*(X-Xi(1))/(Xi(2)-Xi(1));
end
end
%% ==================================================================
% Set the odes
function dydx = rbc_ode(z,y,R,sqTa,k2,Pr,S)
% Zeta = y(1) D1Zeta=y(2) W=y(3) D1W=y(4) D2W=y(5) D3W=y(6) T=y(7) D1T=y(8) 
dydx = [ y(2)
         (k2+S/Pr)*y(1)-sqTa*y(4)
         y(4)
         y(5)
         y(6)
         sqTa*y(2)-k2*(k2+S/Pr)*y(3)+(2*k2+S/Pr)*y(5)+R*k2*y(7)    
         y(8)
         -y(3)+(k2+S)*y(7)];
end
%% ==================================================================
% Set the initial conditions
function yinit = rbc_init(zi,BCs,L)
% y[1]:zeta, y[3]:w, y[7]:T

if BCs == 1 || BCs == 2  % Stress-free
yinit = [      cos(pi*zi/L)
            pi*sin(pi*zi/L)
               sin(pi*zi/L)
            pi*cos(pi*zi/L)
         -pi^2*sin(pi*zi/L)
         -pi^3*cos(pi*zi/L)
          sin(pi*zi/L)
          pi*cos(pi*zi/L) ];

else % No-slip
yinit = [      sin(pi*zi/L)
            pi*cos(pi*zi/L)
               sin(pi*zi/L)
            pi*cos(pi*zi/L)
         -pi^2*sin(pi*zi/L)
         -pi^3*cos(pi*zi/L)
          sin(pi*zi/L)
          pi*cos(pi*zi/L) ];
end
end
%% ==================================================================
% Determine residuals with different boundary conditions
function res = rbc_bc(ya,yb,p,k,BCs)
% Zeta=y(1) D1Zeta=y(2) W=y(3) D1W=y(4) D2W=y(5) D3W=y(6) T=y(7) D1T=y(8)
% 1: FT, Stress-free: y(2)=y(3)=y(5)=y(7)=0
% 2: FF, Stress-free: y(2)=y(3)=y(5)=y(8)=0
% 3: FT, No-slip    : y(1)=y(3)=y(4)=y(7)=0
% 4: FF, No-slip    : y(1)=y(3)=y(4)=y(8)=0
% 5: FT, MV         : At z = 0: y(2)=y(3)=y(5)=y(7)=0    At z = 1: y(1)=y(3)=y(4)=y(7)=0
% 6: FF, MV         : At z = 0: y(2)=y(3)=y(5)=y(8)=0    At z = 1: y(1)=y(3)=y(4)=y(8)=0
% 7: MT, MV (FT,SF) : At z = 0: y(2)=y(3)=y(5)=y(7)=0    At z = 1: y(1)=y(3)=y(4)=y(8)=0

if BCs == 1 % FT, Stress-free
res = [  ya(2) 
         yb(2)
         ya(3)
         yb(3)
         ya(5)
         yb(5)
         ya(7)
         yb(7)
         ya(4)-pi];
elseif BCs == 2 % FF, Stress-free
res = [  ya(2) 
         yb(2)
         ya(3)
         yb(3)
         ya(5)
         yb(5)
         ya(8)
         yb(8)
         ya(4)-1 ];
elseif BCs == 3 % FT, No-slip 
res = [  ya(1) 
         yb(1)
         ya(3)
         yb(3)
         ya(4)
         yb(4)
         ya(7)
         yb(7)
         ya(5)-1 ];
elseif BCs == 4 % FF, No-slip
res = [  ya(1) 
         yb(1)
         ya(3)
         yb(3)
         ya(4)
         yb(4)
         ya(8)
         yb(8)
         ya(5)-1 ];
elseif BCs == 5 % FT, MV
res = [  ya(2) 
         yb(1)
         ya(3)
         yb(3)
         ya(5)
         yb(4)
         ya(7)
         yb(7)
         ya(4)-1 ];
elseif BCs == 6 % FF, MV
res = [  ya(2) 
         yb(1)
         ya(3)
         yb(3)
         ya(5)
         yb(4)
         ya(8)
         yb(8)
         ya(4)-1 ];
else % MT, MV (FT,SF)
res = [  ya(2) 
         yb(1)
         ya(3)
         yb(3)
         ya(5)
         yb(4)
         ya(7)
         yb(8)
         ya(4)-1 ];
end
end