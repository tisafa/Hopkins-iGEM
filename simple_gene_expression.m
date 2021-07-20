% Simple Gene Expression Model

tfin = 60*8; %simlulation final time.
step = 0.1; %simulation step time.
tspan = 0:step:tfin-step;

opti = odeset('AbsTol', 1e-8, 'RelTol', 1e-6); %ODE solver options:
Init = [0 0]; %initial conditions

[t0, x0] = ode23t(@(t,x) model_const_(t,x), tspan, Init, opti)

%Plot
loglog(t0, x0, 'LineWidth',2)
%axis([0 200 -10 300])
legend('mRNA', 'protein')

function [dxdt] = model_const_(t,x)
p_CN = 17; %plasmid nnumber, 17 copies/cell
p_d1 = log(2)/3; %mRNA degradation rate (1/min)
p_d2 = 0.02; %degradation rate (1/min)
p_k2 = 8.23; %translation rate (1/min)
p_k1 = 1.19; %transcription rate (1/min)

% x1 is mRNA
dxdt(1,1) = p_CN*p_k1 - p_d1*x(1);
% x2 is protein
dxdt(2,1) = p_k2*x(1) - p_d2*x(2);
end
