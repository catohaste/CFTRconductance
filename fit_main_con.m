clear;
clc;
%% set global variables 
    global data_x data_y dt;
    
%% set time interval (s) %%
    dt        = .0002; 

%% load observations %%
%% data %%
%     data       = xlsread("C:\Users\Emily\Documents\Lab work\Imaging\2018\eMay\24-5-18 rare muts first half B fsk DMSO\data_180524.xlsx");
data = csvread('~/Desktop/sampleData.csv')';    
data_x     = [0,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,33,35,37,39]';
   data_tofit   = data(1:20,1);
%     data_G551D  = data(1:20,2);
%    data_300nM  = data(1:20,3);
%     data_1uM   = data(1:20,3);
%     data_3uM  = data(1:20,4);
%     data_10uM  = data(1:20,5);
%     data_30uM   = data(1:20,6);
%     data_WTFskVX770   = data(1:20,10);
%     data_3uM   = data(1:20,6);
%     data_10uM   = data(1:20,7);
%    data_30uM   = data(1:20,8);
%     data_VX770   = data(1:20,10);
    data_y     = data_tofit;             % define y values
    
%% set initial parameters %%

x_init    = [-50,10,15,3]';	% x_init = [Vm,G,G_nonCFTR_Cl,TAU_nonCFTR]' 
lb = [-inf,0,0,0]'; 
ub = [0,inf,25,inf]'; 
      
%% output %%
% disp('Running fit...')
%       params         = fminsearchbnd(@fit_error_con,x_init,lb,ub); %includes upper and lower bounds for all parameters
%       Vm            = params(1)
%       G             = params(2)
%       G_nonCFTR_Cl  = params(3) 
%       TAU_nonCFTR   = params(4)
%       error = fit_error_con(params)
     
%% output %%
disp('Running fit...')

% we want to impose the additional condition that A * x <= b

% where A is a matrix, x is our vector of parameters and b is a constant
% vector
tic;
A = [	-5 1 0 0];
b = 450;

Aeq = [];
beq = [];

[params,error] = fmincon(@fit_error_con,x_init,A,b,Aeq,beq,lb,ub);

Vm            = params(1)
G             = params(2)
G_nonCFTR_Cl  = params(3)
TAU_nonCFTR   = params(4)
toc
			
%%
model_out  = fit_transient_con(params(1),params(2),params(3),params(4),dt);
init_out = fit_transient_con(x_init(1),x_init(2),x_init(3),x_init(4),dt);

error_init = fit_error_con(x_init)

model_time = model_out(:,1);
model_P_r_exp = model_out(:,2);

init_P_r_exp = init_out(:,2);


plot(model_time,model_P_r_exp,'-b')
hold on
plot(model_time,init_P_r_exp,'-r')
plot(data_x,data_y,'ok')

%% heatmap
% paramsN = length(params);
% pointsN = 100;
% 
% infReplacement = 1000;
% lb_noInf = [-infReplacement,0,0,0]'; 
% ub_noInf = [0,infReplacement,25,infReplacement]';
% 
% x = zeros(pointsN,paramsN);
% for i=1:paramsN
% 	x(:,i) = linspace(lb_noInf(i),ub_noInf(i),pointsN);
% end

