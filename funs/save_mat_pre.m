%% save mat file, preprocessing

[~,codeVersion] = system('git describe --all');
codeVersion = codeVersion(6:end);
if length(codeVersion) < 6
    codeVersion = strcat(codeVersion,'.0');
end

readme.codeVersion = codeVersion;
readme.matlabVersion = version;
readme.computerVersion = computer;
readme.codePath = mfilename('fullpath');
readme.currentTime = datetime;

% file_path = strcat('./results/',codeVersion,'/',dataType,...
%     '/denoisingK',num2str(dk),'/alpha',num2str(alpha1));
file_path = strcat('./results/',codeVersion,'/',dataType,...
    '/denoisingK',num2str(dk));

if ~exist(file_path, 'dir')
    mkdir(file_path)
end