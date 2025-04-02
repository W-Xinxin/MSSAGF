%% sum_1~V 1/av sum_i~n sum_j~m||Zv_ij-Mv_ij||_2^2 S_ij + alpha||S||_F^2 + lambda||GF'-S||_F^2
%% st, a1 = 1, av>0, G_ij \in{0,1}, ||G_i,:|| =1,  ( F'F = I )

clear all; clc;
addpath('Entropy Rate Superpixel Segmentation');
ds = {'HS-MS Houston2013','HS-SAR Berlin','HS-SAR-DSM Augsburg'}; %  KSC  Pavia  Botswana  Indian   SalinasA
dsPath = './Datasets/';
resPath = './results/';
bestResult = 'bestResults';


for dsi = 3
   fileName = ds{dsi}; disp(fileName);
   addpath( strcat(dsPath,fileName)); 
   data3D ={};
   if strcmp(fileName,'HS-MS Houston2013')
       load  data_HS_LR.mat ; load data_MS_HR.mat; load Houston2013_gt.mat
       data3D{1} = data_HS_LR./max(data_HS_LR(:)) ;  data3D{2}= data_MS_HR./max(data_MS_HR(:));
       gt2D = gnd_flag;
       clear data_HS_LR data_MS_HR  Houstan2013_gt
   elseif strcmp(fileName,'HS-SAR Berlin')
%        data1 = strcat(dsPath,fileName,'/data_HS_LR.mat');      
%        load(data1) ; 
       load data_HS_LR.mat; 
       load data_SAR_HR.mat; load TestImage.mat; load TrainImage.mat
       data3D{1} = data_HS_LR./max(data_HS_LR(:)) ;  data3D{2}= data_SAR_HR./max(data_SAR_HR(:));
       gt2D = TestImage + TrainImage;
       clear data_HS_LR  data_SAR_HR  TestImage  TrainImage
   elseif strcmp(fileName,'HS-SAR-DSM Augsburg')
       load data_DSM.mat; load  data_HS_LR.mat ; load data_SAR_HR.mat; load TestImage.mat; load TrainImage.mat
       data3D{1}=data_DSM./max(data_DSM(:)); data3D{2} = data_HS_LR./max(data_HS_LR(:)) ;  data3D{3}= data_SAR_HR./max(data_SAR_HR(:));
       gt2D = TestImage + TrainImage;
       clear data_DSM  data_HS_LR data_SAR_HR  TestImage TrainImage
   end
end 

txtpath = strcat(resPath,strcat(fileName,'.txt'));
bestpath = strcat(resPath,strcat(bestResult,'.txt'));
if (~exist(resPath,'file'))
    mkdir(resPath);
    addpath(genpath(resPath));
end
dlmwrite(txtpath, strcat('Dataset:',cellstr(fileName), 'Date:',datestr(now)),'-append','delimiter','','newline','pc');
dlmwrite(bestpath, strcat('Dataset:',cellstr(fileName), 'Date:',datestr(now),'-best results'),'-append','delimiter','','newline','pc');

%======================setup=======================
gt = double(gt2D(:));
ind = find(gt);                  % valid instances
numD = length(ind);              % number of samples
numC = length(unique(gt(ind)));  % number of class
numV = size(data3D,2);           % number of views

%% acquire the dim of each view
dim = cell(numV,1);              % band dimension
for iv= 1:numV
    if ndims(data3D{iv}) == 2
        dim{iv} = 1;
    else
        dim{iv} = size(data3D{iv},3);       
    end
end

%% parameter set
spLabel = cell(numV,1);
num_Pixel = cell(numV,1);
lambda_list = [1e4];    % denoising paremeter:        lambda    (for manifold regular term)
d_nk = 5;               % neighbor number for denoising
l_nklist = [40];        % neighbor number for local anchor graph construction
beta_set = [1e-5];      % model trade-off parameter:  beta
num_gM_set = [6*numC];  % anchor number for global anchor graph
g_nk =25;               % neighbor number for global anchor graph 

%% denoising processing 
for ib = 1: length(lambda_list)
    lambda = lambda_list(ib);
    for iv = 1: numV
       [dX{iv},spLabel{iv},num_Pixel{iv}] = preData(data3D{iv},lambda,d_nk);   % local smooth
    end

    %% construct anchor graph with different scale for each view 
    %newData2D = mat2cell(X,numD,cell2mat(dim)); % convert entire 2D matrix data into 2d matrix of each view
    graph =cell(numV,1);   % anchor graph for each view
    A = cell(numV,1);      % anchor for each view
    
    %% generated anchor from each ERS block by average operator and construct anchor graph
    for in = 1: length(l_nklist)
        l_nk = l_nklist(in);
        for iv = 1: numV
           newData2D = dX{iv};
           A{iv} = meanInd(newData2D,spLabel{iv});                % anchor in ERS block dim * m for each view
           [graph{iv}, ~] = Gen_Achor_Adj( newData2D,A{iv},l_nk,0); % issymmetric  local anchor graphs
        end
         
        %% main procedure     
        best_M2beta = zeros(1,6);  % record best beta
        for jm=1:length(num_gM_set);
           %% inilization for global anchors
            num_gM = num_gM_set(jm);  %% anchor number  for global anchor graph
            M = cell(numV,1);         %% global anchor matrix
            parfor iv = 1: numV
                rng(5489,'twister');
                [labels{iv}, M{iv}] = litekmeans(graph{iv}, num_gM, 'MaxIter', 100,'Replicates',1);%augsurb and berlin 100, houston 70
                M{iv} = M{iv}';     % dv * m         
            end
            
            for ia = 1: length(beta_set)
                numD = size(graph{1},1);    % n*dv        number of samples
                numA = size(M{1},2);        % dv *numM    number of global anchors      
                beta = beta_set(ia);        % parameter
                av = ones(numV,1)/numV;     % view weights                    
                k = g_nk;                     % set the neighbor number for global anchor graph 
                converged = 0;
                iter =0;                    
                Anchor = M;                 % acquire initial global anchor 
                                
                while converged ~=1 && iter < 30   % 10   
                    %% update S
                    sum_distX = zeros();
                    for iv = 1: numV
                         X = graph{iv}';  % dv * n
                         %anchor = M{iv}; % dv * m
                         anchor = Anchor{iv};
                         distX{iv} = L2_distance_1(X, anchor);
                        % [~, idx{iv}] = sort(distX{iv}, 2);
                         sum_distX = sum_distX + (1/av(iv))* distX{iv};
                    end
                    [~, idx] = sort(sum_distX, 2);

                    if iter < 1  % to initialize S
                       %% update S 
                       S = zeros(numD,numA);
                       for i =1 : numD
                           id = idx(i,1:k+1);
                           di = sum_distX(i,id);
                           S(i,id) =  (di(k+1)-di)/(k*di(k+1)-sum(di(1:k))+eps);
                       end
                    elseif iter > 0
                       %% to update S  
                        S = zeros(numD,numA);
                        U = G * F' ;    % ||GF' - S ||_F^2
                        for i = 1:numD
                            id = idx(i,1:k+1);
                            di = sum_distX(i, id);
                            numerator = di(k+1)-di+2*beta*U(i,id(:))-2*beta*U(i,id(k+1));
                            denominator1 = k*di(k+1)-sum(di(1:k));
                            denominator2 = 2*beta*sum(U(i,id(1:k)))-2*k*beta*U(i,id(k+1));
                            S(i,id) = max(numerator/(denominator1+denominator2+eps),0);
                        end;
                    end
                    
                    %% update F and G
                    if iter < 1   % to initialize G,F with initial S
                         rng(5489,'twister');
                         [Y, Center] = litekmeans(S, numC, 'MaxIter', 100,'Replicates',1);
                         F = Center'; % numA * numC
                         G = zeros(numD,numC);
                         for i = 1: numD
                             dex = Y(i);
                             G(i,dex)= 1;
                         end
                         G_old = G;
                         F_old = F;

                    elseif iter > 0
                         %% update F  center matrix F'F = I
                         F_old = F;
                         C = S'*G;
                         [U,~,V] = svd(C,'econ');
                         F = U*V';  % 
                         %% update G indicator matirx
                         G = zeros(numD,numC); %n * m
                         for i  = 1: numD
                             Dis=zeros(numC,1);
                             instance = S(i,:)';
                             for j=1:numC         
                                Dis(j)=(norm(instance-F(:,j)))^2;
                             end
                             [~,r]=min(Dis);
                             G(i,r(1))=1;           
                         end
                    end 
                    
                    %% update M{iv}  -> tr(ZPZ'-2ZSM'+ MHM')
                     Anchor_old = Anchor;  %% preserve the old value
                     H = diag(sum(S,1)); % m*m
                     for iv = 1: numV
%                          M{iv} = graph{iv}' * S * H^(-1) ;  % d*m  M =(Z'*S*H^-1)
                           Anchor{iv} = graph{iv}' * S * H^(-1) ;
                     end
                     
                    %% update av
                    av_old =av;  % original     
                    J = zeros(numV,1);
                    for iv = 1: numV
                         X = graph{iv}';  % dv * n
                         %anchor = M{iv};  % dv * m
                         anchor = Anchor{iv};
                         Dxa{iv} = L2_distance_1(X, anchor);
                         J(iv) = sqrt(trace(Dxa{iv}'* S));
                    end
                    av = J / sum(J);
                    res_av =  norm(av-av_old,'fro');
                    res_av_list(iter+1) = res_av;
                                        
                    value = norm(F-F_old,'fro');
                    obj(iter+1) = value;
                    %fprintf('loss: \t%.6f\n',value);
                    if  iter>3 && res_av < 1e-5;
                       converged = 1;
                    end 
                     iter = iter+ 1;                                
                end    %% for optimization iter range
                
                %% obtain indicators
                [~,Pre_idx]= max(G');
                results = evaluate_results_clustering(gt(ind),Pre_idx(ind));
                dlmwrite(txtpath ,[results,lambda,l_nk,beta,num_gM,g_nk],'-append','delimiter','\t','newline','pc');
                fprintf('result:\t beta:%f l_nk:%d lambda:%.6f num_gM:%d g_nk:%d ACC:%12.6f Kappa:%12.6f NMI:%12.6f Purity:%12.6f\n',[lambda,l_nk,beta,num_gM,g_nk, results(1),results(2),results(3),results(4)]);               
                if results(1) > best_M2beta(1)
                    best_M2beta = [results,lambda,l_nk,beta,num_gM,g_nk];
                end
            end  %% for range of lambda
        end   %% for range of numM
        dlmwrite(bestpath,[best_M2beta],'-append','delimiter','\t','newline','pc');
    end  %% for range of nk
end %% for range of beta
