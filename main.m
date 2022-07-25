clear;
clc;
%% set global variables
global data_x data_y dt;

%% set time interval (s) %%
dt        = .0002;

%% load data %%
data_filename = '~/Desktop/forskolin28/quench/meanTimeline_2018_08_25_2105.csv';

data = importdata(data_filename,',');
data.textdata(2:end,5:end) = num2cell(data.data);

data_num = data.data;

data_time = [0,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,33,35,37,39]';
data_YFP_unbound_norm = data_num(:,2:end);

data_x = [0,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,33,35,37,39]';

conditionN = size(data_num,1);
conditions = cell(conditionN,1);
for i=1:conditionN
	conditions(i) = strcat(data.textdata(1+i,2),{' '},data.textdata(1+i,3),{' '},data.textdata(1+i,4));
end

%% set initial parameters %%

% params = [ Vm , G , G_nonCFTR_Cl , TAU_nonCFTR ]

% initialize parameters
params	= zeros(conditionN,4);
error		= zeros(conditionN,1);

% set starting conditions
params_init    = [-50,10,15,3]';

% define upper and lower bounds for our parameter values
% lb = [-150,0,0,0]'; 
% ub = [0,200,25,15]'; 
lb =  [-150,0,13.2,4.9]';
ub =  [0,200,13.2,4.9]';

% we also want to impose the additional condition that A * params <= b
% where A is a matrix, x is our vector of parameters and b is a constant
% vector
A = [	-5 1 0 0];
b = 450;

% Aeq * params = beq is not needed
Aeq = [];
beq = [];

%% output
disp('Running fit...')

% run the fitting for each condition
for i=1:conditionN
	
	tic;
	
	data_y = data_YFP_unbound_norm(i,:)';

	[params(i,:),error(i,:)] = fmincon(@fit_error_con,params_init,A,b,Aeq,beq,lb,ub);

% 	Vm            = params(1);
% 	G             = params(2);
% 	G_nonCFTR_Cl  = params(3);
% 	TAU_nonCFTR   = params(4);

	msg = [num2str(i),' of ',num2str(conditionN),' conditions completed.'];
	disp(msg)
	toc
	
end

%% output results

data_header = data.textdata(1,1:5);
params_header = { 'Vm' , 'G' , 'G_nonCFTR_Cl' , 'TAU_nonCFTR' };

header = horzcat(data_header,params_header);
output = horzcat(data.textdata(2:end,1:5),num2cell(params));

tableOutput = cell2table(output,'VariableNames',header);
cellArrayOutput = vertcat(header,output);

if ispc == true
	outputFilename = '~/Desktop/resultsQuenchModel/all_exp_forskolin_modelOutput_2fixed.xlsx';
	xlswrite(outputFilename,cellArrayOutput)
elseif isunix == true
	outputFilename = '~/Desktop/resultsQuenchModel/all_exp_forskolin_modelOutput_2fixed.csv';
	writetable(tableOutput,outputFilename);
end

disp('Finished writing output to file')


%%
% init_out = fit_transient_con(x_init(1),x_init(2),x_init(3),x_init(4),dt);
% init_P_r_exp = init_out(:,2);

%% plot results

close all

model_time = 0:dt:40;
model_time_pointsN = length(model_time);
model_YFP_unbound_norm = zeros(conditionN,model_time_pointsN);

% subplotIdx = [1,8,2,9,3,10,4,11,5,12,6,13,7,14];
subplotIdx = 1:18;

for i=1:conditionN

	model_out  = fit_transient_con(params(i,1),params(i,2),params(i,3),params(i,4),dt);
	model_YFP_unbound_norm(i,:) = model_out(:,2);

	ax = subplot(2,7,subplotIdx(i));
		plot(ax,model_time,model_YFP_unbound_norm(i,:),'-b')
		hold on
		plot(ax,data_time,data_YFP_unbound_norm(i,:),'ok')
		ylim([0,1])
		title(ax,sprintf(conditions{i}))
% 		legend(ax,{'model','data'})

end

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

