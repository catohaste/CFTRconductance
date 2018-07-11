
function error = fit_error_con(F)

    global data_x data_y dt; % use global variables defined in fit_main
   
    Vm              = F(1);
    G               = F(2);
    G_nonCFTR_Cl    = F(3);
    TAU_nonCFTR     = F(4);
    pred            = fit_transient_con(Vm,G,G_nonCFTR_Cl,TAU_nonCFTR,dt); 


for i=1:length(data_x)
    x           = pred(:,1)==data_x(i);
    y           = sum(pred.*x);
    model_y(i,1)= y(2);
end
error           = sum((model_y-data_y).^2);
end
