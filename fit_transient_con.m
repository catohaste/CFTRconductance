function [output] = fit_transient_con(Vm,G,G_nonCFTR_Cl,TAU_nonCFTR,dt)

%% INITIAL PARAMATERS %%
% Binding affinity of I- and Cl^{-} to YFP (mM)
    Ki              = 1.9;
    Kcl             = 85;
% Concentration inside and outside the cell (mM)
    K_in            = 4.7/(10^(Vm/62));
    K_out           = 4.7;
    I_out           = 100; 
    I_in            = 1.0000e-4;
    Cl_out          = 117.1;
    Cl_in           = 152/(10^(Vm/-62));        % [Cl^{-}]in at eq with 152 mM [Cl^{-}]out 
% Proportion of YFP bound (Ir, Clr) and unbound (r)
    P_Ir            = I_in/(Ki*(1+Cl_in/Kcl)+I_in);
    P_Clr           = Cl_in/(Kcl*(1+I_in/Ki)+Cl_in);
    P_r             = 1-P_Ir-P_Clr;
    P_Ir_Clr_ratio  = .83;                      % 0.83*GclMax Lindsdell01
    P_r_exp         =(1-P_Ir-P_Clr)/P_r;        %normalised to t0
% conductance (nS)
    G_CFTR_I        = P_Ir_Clr_ratio*G;
    G_K_leak        = 2.5;
%  G_nonCFTR_Cl   = 0;
% equilibrian potentials (mV)
    E_K             =  62*(log10(K_out/K_in));
    E_Cl            = -62*(log10(Cl_out/Cl_in));
    E_I             = -62*(log10(I_out/I_in));
% charachteristics HEK cells (L, cm^2, uF)
    HEK_VOL         = 4/3*pi*((7*10^-6)^3)*10^3; % CHECK DIAMETER AND RADIUS OF HEKS
    HEK_SURF        = 1.5*4*pi*((7*10^-6)^2)*10^4;
    C_HEK           = 1*HEK_SURF;
% other parameters
% TAU_nonCFTR   = 0;       % tau (s)
    time            = (0:dt:40)';
    exp_decay       = 1*exp(-time/TAU_nonCFTR);
    F               = 96500;    % Faraday constant (C/mol)
    T               = 310;      % Temperature (Kelvin)
    R               = 8.31;     % J*(deg*mol)^-1
    z_K             = 1; 
    z_I             = -1;
    z_Cl            = -1;
    Gclnonmax       = G_nonCFTR_Cl;
    Gclnon          = G_nonCFTR_Cl*exp_decay(1);
    Ginon           = G_nonCFTR_Cl*P_Ir_Clr_ratio*exp_decay(1);
% current carried by Cl^{-}, K+ and I-(pA)
    I_Cl            =  (G/140)*Vm*(Cl_in-Cl_out*exp(-z_Cl*F*0.001*Vm/(R*T)))...
                   /(1-exp(-z_Cl*F*0.001*Vm/(R*T)))+(Gclnon/140)*Vm*...
                   (Cl_in-Cl_out*exp(-z_Cl*F*0.001*Vm/(R*T)))...
                   /(1-exp(-z_Cl*F*0.001*Vm/(R*T)));
    I_I             =  (G_CFTR_I/140)*Vm*(I_in-I_out*exp(-z_I*F*0.001*Vm/(R*T)))...
                   /(1-exp(-z_I*F*0.001*Vm/(R*T)))+(Ginon/140)*Vm*...
                   (I_in-I_out*exp(-z_I*F*0.001*Vm/(R*T)))...
                   /(1-exp(-z_I*F*0.001*Vm/(R*T)));
    I_K             =  (G_K_leak/140)*Vm*(K_in-K_out*exp(-z_K*F*0.001*Vm/(R*T)))...
                   /(1-exp(-z_K*F*0.001*Vm/(R*T))); % Gk = 2.5 ns Rapedius06
               
    error           = Ki*Cl_in/Kcl;
    I_estimate      = I_in+error;

       
    
   
%% TIMECOURSE %%
for t = 2:length(time)
% New concentration inside the cell (mM)
    Cl_in(t,1)      = Cl_in(t-1)+((time(t)-time(t-1))...
                    *I_Cl(t-1)/F/HEK_VOL)*10^-9;
    I_in(t,1)       = I_in(t-1)+((time(t)-time(t-1))...
                    *I_I(t-1)/F/HEK_VOL)*10^-9;
    K_in(t,1)       = K_in(t-1)+((time(t)-time(t-1))...
                    *I_K(t-1)/F/HEK_VOL)*10^-9;
% New Vm (mV)
    Vm(t,1)         = Vm(t-1) +((time(t)-time(t-1))...
                    *-(I_Cl(t-1)+I_I(t-1)+I_K(t-1))/C_HEK)*10^-3; %changed -I_K to +I_K
% New equilibrian potentials (mV)
    E_K(t,1)        =  62*(log10(K_out/K_in(t)));
    E_Cl(t,1)       = -62*(log10(Cl_out/Cl_in(t)));
    E_I(t,1)        = -62*(log10(I_out/I_in(t)));
% New transient conductances (nS)    
    Gclnon          = G_nonCFTR_Cl*exp_decay(t);
    Ginon           = G_nonCFTR_Cl*P_Ir_Clr_ratio*exp_decay(t);
% New current carried by Cl^{-}, K+ and I-(pA)
    I_Cl(t,1)      = (G/140)*Vm(t)*(Cl_in(t)-Cl_out...
                    *exp(-z_Cl*F*0.001*Vm(t)/(R*T)))...
                   /(1-exp(-z_Cl*F*0.001*Vm(t)/(R*T)))+(Gclnon/140)...
                   *Vm(t)*(Cl_in(t)-Cl_out*exp(-z_Cl*F*0.001*Vm(t)/(R*T)))...
                   /(1-exp(-z_Cl*F*0.001*Vm(t)/(R*T)));
    I_I(t,1)       = (G_CFTR_I/140)*Vm(t)*(I_in(t)-I_out...
                    *exp(-z_I*F*0.001*Vm(t)/(R*T)))...
                   /(1-exp(-z_I*F*0.001*Vm(t)/(R*T)))+(Ginon/140)*Vm(t)*...
                   (I_in(t)-I_out*exp(-z_I*F*0.001*Vm(t)/(R*T)))...
                   /(1-exp(-z_I*F*0.001*Vm(t)/(R*T)));
    I_K(t,1)       =  (G_K_leak/140)*Vm(t)*(K_in(t)-K_out...
                    *exp(-z_K*F*0.001*Vm(t)/(R*T)))...
                   /(1-exp(-z_K*F*0.001*Vm(t)/(R*T))); % Gk = 2.5 ns Rapedius06
% New [I-] estimate and error           
    error(t,1)      = Ki*Cl_in(t)/Kcl;
    I_estimate(t,1) = I_in(t)+error(t);
% Proportion of YFP
    P_Ir(t,1)       = I_in(t)/(Ki*(1+Cl_in(t)/Kcl)+I_in(t));
    P_Clr(t,1)      = Cl_in(t)/(Kcl*(1+I_in(t)/Ki)+Cl_in(t));
    P_r_exp(t,1)    = (1-P_Ir(t)-P_Clr(t))/P_r(1);
    P_r(t,1)        = 1-P_Ir(t)-P_Clr(t);
% Rate of iodide entry (d[I-]/dt)
    rate_I_entry(t,1)           = (I_in(t)-I_in(t-1))/(time(t)-time(t-1));
    rate_I_entry_estimate(t,1)  = (I_estimate(t)-I_estimate(t-1))/(time(t)-time(t-1));
end

output = [time P_r_exp];
