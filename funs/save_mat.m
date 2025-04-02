%% save mat file

file_name = strcat(dataType,'_proj',num2str(proj),...
    '_ACC',num2str(result(1),4),'.mat');

result_file = fullfile(file_path, file_name);

%%
clear file_name;

save(result_file,'dataType','result','proj','y_pred','clusternum',...
    'alpha1','readme','dk');


