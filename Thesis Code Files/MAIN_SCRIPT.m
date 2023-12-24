% MAIN_SCRIPT tests a variety of eigenvalue solvers (most of then inverse
% free) in their ability to solver for k << n eigenpairs of the symmetric
% positive definite generalized eigenvalue problem

% A*X = B*X*Lambda,

% where A and B are both symmetric matrices of order n, B is positive
% definite, and Lambda is a diagonal matrix of order k. The columns of the
% n x k matrix X are B-orthonormal, i.e.

% X'*B*X = I in R^{k x k}.

% Required files: RUN_COMPARISON, bleigifp.m, lobpcg.m, geigpen.m, and
% Matrix Pencils directory containing all the relevant Matrix Pencil
% structs.

% Copyright 2016, Frankie Camacho, Rice University
clear; close all; clc;

size_limit = 40000; %Limit the size of the problems that we tackle

% This directory stores MAIN_SCRIPT.m, RUN_COMPARISON.m, bleigifp.m,
% geigpen.m, lobpcg.m, and the folder Matrix_Pencils containing all the
% relevant matrix pencils mat files in folder titled Matrix_Pencil_j, where
% j is an integer between 1 and 32.

FunctionDir = '/Users/frankiecamacho/Desktop/Final_Test_Code_GEIGPEN';
% FunctionDir = '/home/fc16/Desktop/Final_Test_Code_GEIGPEN/';

% K is a vector containing the number of eigenpairs that we want to
% compute.
% K = [10:5:25];
K = [20];
numK = numel(K);

% The entries of MP_index must be in ascending algebraic order.
MP_index = [26];
% MP_index = [1, 13, 19, 26:32];
% There must be a folder names Matrix_Pencil_MaxMPI in 
% FunctionDir/'Matrix Pencils'/.
MaxMPI = max(MP_index); 
numMP_index = numel(MP_index);


% Use solver_switch to decide if a particular eigenvalue solver will be a
% part of the experiment.

solver_switch = [];
solver_switch.eigs_sm = 1;
solver_switch.eigs_sa = 0;
solver_switch.bleigifp = 0;
solver_switch.lobpcg = 0;
solver_switch.geigpen = 1;

% Use matDim_ParamDivide to control the iteration parameters.

matdim_ParamDivide = 1000;
% We will use the following for plotting our results.

structFind = @(x) x > 0; % Find which eigenvalue solvers actually were used.

% For the legend.

entries = {'Shifted eigs-SM','eigs-SA','bleigifp','lobpcg','geigpen'};
entries = entries(find(structfun(structFind,solver_switch)));
% pause;

C = [     0,   0.4470,   0.7410; % C is a color matrix
     0.8500,   0.3250,   0.0980;
     0.9290,   0.6940,   0.1250;
     0.4940,   0.1840,   0.5560;
     0.4660,   0.6740,   0.1880];

 

% common_tol is a vector of different tolerances. We will iterate through
% common_tol, the different matrix pencils, and the different values for
% k order to run our experiments.


for common_tol = [1e-4];
    
    % j iterates through the desired matrix pencils that are specified in
    % MP_index.
    
    J = 1;
    j = MP_index(J);
    
    
    while (j<=MaxMPI)&&(J<=numMP_index)
        %% Here we load the appropriate matrix pencil.
        
        while (j<=MaxMPI)&&(J<=numMP_index)
            
            % Loading the matrix pencils and making sure they satisfy our
            % conditions.
            BaseRoot = ['/Users/frankiecamacho/Desktop/',...
                'Final_Test_Code_GEIGPEN/Matrix_Pencils/Matrix_Pencil_',...
                num2str(j),'/'];
%             BaseRoot = ['/home/fc16/Desktop/Final_Test_Code_GEIGPEN/',...
%                 'Matrix_Pencils/Matrix_Pencil_',num2str(j),'/'];
%             %             cd(BaseRoot);
            load([BaseRoot,'Matrix_Pencil_',num2str(j)]);
            eval(['BMS = Matrix_Pencil_',num2str(j),';']);
            % BMS stands for BaseMatrixStruct.
            clear(['Matrix_Pencil_',num2str(j)]);
            fprintf('MAIN_SCRIPT loaded %s.\nA = %s --- B = %s.\n',...
                BMS.Name,BMS.A_name,BMS.B_name);
            
            % Making sure all of our criteria is met, if not load the next
            % matrix. If yes, save the program output in an ASCII file.
            
            n = BMS.A_size(2);
            if all([BMS.A_size, BMS.B_size] == n)&&(n < size_limit)
                if exist([BaseRoot,'Matrix_Pencil_',num2str(j),' Diary',...
                        ', K = [',num2str(K),'], tol = ',...
                        num2str(common_tol)],'file');
                    delete([BaseRoot,'Matrix_Pencil_',num2str(j),...
                        ' Diary, K = [',num2str(K),'], tol = ',...
                        num2str(common_tol)]);
                end
                diary([BaseRoot,'Matrix_Pencil_',num2str(j),' Diary, K',...
                    ' = [',num2str(K),'], tol = ',...
                    num2str(common_tol)]);
                fprintf('Matrix pencil load was successful.\n\n');
                
                break;
            else
                if n>size_limit
                    fprintf(['The matrix size n = %i is\nabove the',...
                        ' prescribed limit M = %i.\n'],...
                        n,size_limit);
                else
                    fprintf(['Matrix dimensions do not match.\nEither',...
                        ' one or both are square or\nunequal to each ',...
                        'other, or both.\n']);
                end
                
                J = J + 1;
                if J>numMP_index
                    break;
                else
                    j = MP_index(J);
                    if j >MaxMPI;
                        break;
                    end
                end
                fprintf('Will attempt to load the next matrix.\n\n');
            end
            
        end
        
        if j>MaxMPI || J>numMP_index
            break;
        end
        
        %% Choosing the iteration parameters for our solvers.
        
        switch (n>=matdim_ParamDivide)
            case 1
                common_maxit = 1000;
                geigpen_maxitO = 10;
                geigpen_maxitI = 100;
            case 0
                common_maxit = 500;
                geigpen_maxitO = 5;
                geigpen_maxitI = 100;
            otherwise
                warning('Common iteration default parameters assigned.');
        end
        
%% Preallocating vectors and constants to store historical data for graphs.
        
        index = []; T = []; MRE = []; trXAX = []; pnrm2 = []; ERE = [];
        EREplot = []; iter = [];
        
        if solver_switch.eigs_sm
            index.eigs_sm = 0;
            T.eigs_sm = zeros(numK,1);
            MRE.eigs_sm = zeros(numK,1);
            trXAX.eigs_sm = zeros(numK,1);
            pnrm2.eigs_sm = zeros(numK,1);
            out_eigs_sm = [];
            out_eigs_sm.success = 1;
        end
        
        if solver_switch.eigs_sa
            index.eigs_sa = 0;
            T.eigs_sa = zeros(numK,1);
            MRE.eigs_sa = zeros(numK,1);
            trXAX.eigs_sa = zeros(numK,1);
            pnrm2.eigs_sa = zeros(numK,1);
            ERE.eigs_sa = zeros(numK,1);
            EREplot.eigs_sa = 0;
            out_eigs_sa = [];
            out_eigs_sa.success = 1;
        end
        
        if solver_switch.bleigifp
            index.bleigifp = 0;
            T.bleigifp = zeros(numK,1);
            MRE.bleigifp = zeros(numK,1);
            trXAX.bleigifp = zeros(numK,1);
            pnrm2.bleigifp = zeros(numK,1);
            ERE.bleigifp = zeros(numK,1);
            EREplot.bleigifp = 0;
            out_bleigifp = [];
            out_bleigifp.success = 1;
        end
        
        if solver_switch.lobpcg
            index.lobpcg = 0;
            T.lobpcg = zeros(numK,1);
            MRE.lobpcg = zeros(numK,1);
            trXAX.lobpcg = zeros(numK,1);
            pnrm2.lobpcg = zeros(numK,1);
            ERE.lobpcg = zeros(numK,1);
            EREplot.lobpcg = 0;
            iter.lobpcg = zeros(numK,1);
            out_lobpcg = [];
            out_lobpcg.success = 1;
            
        end
        
        if solver_switch.geigpen
            index.geigpen = 0;
            T.geigpen = zeros(numK,1);
            MRE.geigpen = zeros(numK,1);
            trXAX.geigpen = zeros(numK,1);
            pnrm2.geigpen = zeros(numK,1);
            ERE.geigpen = zeros(numK,1);
            EREplot.geigpen = 0;
            iter.geigpen = zeros(numK,3);
            out_geigpen = [];
            out_geigpen.success = 1;
        end
        
        
        for m = 1:numK
            %% Running the experiments for all m in K.
            
            k = K(m);
            
            %%%%%%%%%%% Assigning the options for eigs_sm %%%%%%%%%%%
            opts_eigs_sm = [];
            if solver_switch.eigs_sm
                if out_eigs_sm.success
                    sigma = 'SM';
                    opts_eigs_sm.proceed = 1;
                    opts_eigs_sm.issym = true;
                    opts_eigs_sm.isreal = true;
                    opts_eigs_sm.disp = 0;
                    opts_eigs_sm.tol = common_tol;
                    opts_eigs_sm.maxit = common_maxit;
                    opts_eigs_sm.p = 2*k;
                else
                    opts_eigs_sm.proceed = 0;
                end
            end
            
            %%%%%%%%%%% Assigning the options for eigs_sa %%%%%%%%%%%
            opts_eigs_sa = [];
            if solver_switch.eigs_sa
                if out_eigs_sa.success
                    sigma = 'SA';
                    opts_eigs_sa.proceed = 1;
                    opts_eigs_sa.issym = true;
                    opts_eigs_sa.isreal = true;
                    opts_eigs_sa.disp = 0;
                    opts_eigs_sa.tol = common_tol;
                    opts_eigs_sa.maxit = common_maxit;
                    opts_eigs_sa.p = 2*k;
                else
                    opts_eigs_sa.proceed = 0;
                end
            end
            
            %%%%%%%%%%% Assigning the options for BLEIGIFP %%%%%%%%%%%
            opts_bleigifp = [];
            if solver_switch.bleigifp
                if out_bleigifp.success
                    opts_bleigifp.proceed = 1;
                    opts_bleigifp.DISP = 0;
                    %opts_bleigifp.ADAPTTOL = 0.1;
                    %opts_bleigifp.BLOCKSIZE = 2;
                    opts_bleigifp.SIGMA = 'SA';
                    %opts_bleigifp.MAXITERATION = [];
                    %opts_bleigifp.INNERITERATION = [];
                    opts_bleigifp.TOLERANCE = common_tol;
                    %opts_bleigifp.NORMA = norm(A,1);
                    %opts_bleigifp.NORMB = norm(B,1);
                    opts_bleigifp.SIZE = BMS.A_size(2);
                    %opts_bleigifp.ILUTHRESH = 1e-4;
                    %opts_bleigifp.MAXBLOCKSIZE = 5*opts_bleigifp.BLOCKSIZE;
                    %opts_bleigifp.INITIALVEC = randn(n,k);
                    %opts_bleigifp.UPDATEM = [];
                    %opts_bleigifp.UPDATEP = [];
                else
                    opts_bleigifp.proceed = 0;
                end
            end
            
            %%%%%%%%%%% Assigning the options for LOBPCG %%%%%%%%%%%
            opts_lobpcg = [];
            if solver_switch.lobpcg
                if out_lobpcg.success
                    opts_lobpcg.proceed = 1;
                    opts_lobpcg.residualTolerance = common_tol;
                    opts_lobpcg.maxIterations = common_maxit;
                    opts_lobpcg.verbosityLevel = 0;
                else
                    opts_lobpcg.proceed = 0;
                end
            end
            
            %%%%%%%%%%% Assigning the options for GEIGPEN %%%%%%%%%%%
            opts_geigpen = [];
            if solver_switch.geigpen
                if out_geigpen.success == 1
                    opts_geigpen.proceed = 1;
                    opts_geigpen.arrp = 0;
                    opts_geigpen.proj = 1;
                    opts_geigpen.print = 1;
                    opts_geigpen.kp = ceil(.2*k);
                    opts_geigpen.tol1 = 0.1;
                    opts_geigpen.tol2 = common_tol;
                    opts_geigpen.alf_bool = 1;
                    opts_geigpen.maxit = geigpen_maxitI;
                    opts_geigpen.maxouter = geigpen_maxitO;
                    opts_geigpen.dols = 2;
                    opts_geigpen.istop = false;
                    opts_geigpen.freq = 5;
                    opts_geigpen.rho = 0.2;
                    opts_geigpen.tol1sp = 5;
                else
                    opts_geigpen.proceed = 0;
                end
            end
            
            % Need to make sure that our solver actually runs by making all
            % the inputs are defined.
            if exist('opts_eigs_sm','var')
                if ~isfield(opts_eigs_sm,'proceed')
                    opts_eigs_sm.proceed = 0;
                elseif isempty(opts_eigs_sm.proceed)
                    opts_eigs_sm.proceed = 0;
                end
            else
                opts_eigs_sm = [];
                opts_eigs_sm.proceed = 0;
            end
            
            
            if exist('opts_eigs_sa','var')
                if ~isfield(opts_eigs_sa,'proceed')
                    opts_eigs_sa.proceed = 0;
                elseif isempty(opts_eigs_sa.proceed)
                    opts_eigs_sa.proceed = 0;
                end
            else
                opts_eigs_sa = [];
                opts_eigs_sa.proceed = 0;
            end
            
            
            if exist('opts_bleigifp','var')
                if ~isfield(opts_bleigifp,'proceed')
                    opts_bleigifp.proceed = 0;
                elseif isempty(opts_bleigifp.proceed)
                    opts_bleigifp.proceed = 0;
                end
            else
                opts_bleigifp = [];
                opts_bleigifp.proceed = 0;
            end
            
            
            if exist('opts_lobpcg','var')
                if ~isfield(opts_lobpcg,'proceed')
                    opts_lobpcg.proceed = 0;
                elseif isempty(opts_lobpcg.proceed)
                    opts_lobpcg.proceed = 0;
                end
            else
                opts_lobpcg = [];
                opts_lobpcg.proceed = 0;
            end
            
            
            if exist('opts_geigpen','var')
                if ~isfield(opts_geigpen,'proceed')
                    opts_geigpen.proceed = 0;
                elseif isempty(opts_geigpen.proceed)
                    opts_geigpen.proceed = 0;
                end
            else
                opts_geigpen = [];
                opts_geigpen.proceed = 0;
            end
          
            
            %%%%%%% Calling RUN_COMPARISON to perform experiments  %%%%%%%
            
%             cd(FunctionDir);
            [out_eigs_sm, out_eigs_sa, out_bleigifp, out_lobpcg,...
                out_geigpen] = RUN_COMPARISON(BMS,k,opts_eigs_sm,...
                opts_eigs_sa,opts_bleigifp,opts_lobpcg,opts_geigpen,...
                solver_switch,BaseRoot,FunctionDir,entries,structFind);
            
            eigenValues = figure;
            
            if solver_switch.eigs_sm
                if out_eigs_sm.success
                    index.eigs_sm = index.eigs_sm+1;
                    T.eigs_sm(index.eigs_sm) = out_eigs_sm.T(3);
                    pnrm2.eigs_sm(index.eigs_sm) = out_eigs_sm.pnrm2;
                    MRE.eigs_sm(index.eigs_sm) = out_eigs_sm.maxres;                    
                    trXAX.eigs_sm(index.eigs_sm) = out_eigs_sm.trXAX;
                    ev = out_eigs_sm.ev;
                    plot(real(ev),'color',C(1,:),'LineStyle','none',...
                        'LineWidth',1.2,'marker','d');
                    hold on
                end
            end
            
            if solver_switch.eigs_sa
                if out_eigs_sa.success
                    index.eigs_sa = index.eigs_sa+1;
                    T.eigs_sa(index.eigs_sa) = out_eigs_sa.T;
                    pnrm2.eigs_sa(index.eigs_sa) = out_eigs_sa.pnrm2;
                    MRE.eigs_sa(index.eigs_sa) = out_eigs_sa.maxres;
                    trXAX.eigs_sa(index.eigs_sa) = out_eigs_sa.trXAX;
                    if isfield(out_eigs_sa,'evrelerr')
                        EREplot.eigs_sa = EREplot.eigs_sa + 1;
                        ERE.eigs_sa(EREplot.eigs_sa)=out_eigs_sa.evrelerr;
                    end
                    ev = out_eigs_sa.ev;
                    figure(eigenValues);
                    plot(real(ev),'color',C(2,:),'LineStyle','none',...
                        'LineWidth',1.2,'marker','o');
                    hold on
                    
                end
            end
                        
            if solver_switch.bleigifp
                if out_bleigifp.success
                    index.bleigifp = index.bleigifp+1;
                    T.bleigifp(index.bleigifp) = out_bleigifp.T;
                    pnrm2.bleigifp(index.bleigifp) = out_bleigifp.pnrm2;
                    MRE.bleigifp(index.bleigifp) = out_bleigifp.maxres;
                    trXAX.bleigifp(index.bleigifp) = out_bleigifp.trXAX;
                    if isfield(out_bleigifp,'evrelerr')
                        EREplot.bleigifp = EREplot.bleigifp+1;
                        ERE.bleigifp(index.bleigifp)=out_bleigifp.evrelerr;
                    end
                    figure(eigenValues);
                    ev = out_bleigifp.ev;
                    plot(real(ev),'color',C(3,:),'LineStyle','none',...
                        'LineWidth',1.2,'marker','s');
                    hold on
                end
            end
            
            if solver_switch.lobpcg
                if out_lobpcg.success
                    index.lobpcg = index.lobpcg+1;
                    T.lobpcg(index.lobpcg) = out_lobpcg.T;
                    pnrm2.lobpcg(index.lobpcg) = out_lobpcg.pnrm2;
                    MRE.lobpcg(index.lobpcg) = out_lobpcg.maxres;
                    trXAX.lobpcg(index.lobpcg) = out_lobpcg.trXAX;
                    if isfield(out_lobpcg,'evrelerr')
                        EREplot.lobpcg = EREplot.lobpcg + 1;
                        ERE.lobpcg(index.lobpcg) = out_lobpcg.evrelerr;
                    end
                    iter.lobpcg(index.lobpcg) = out_lobpcg.iter;
                    figure(eigenValues);
                    ev = out_lobpcg.ev;
                    plot(real(ev),'color',C(4,:),'LineStyle','none',...
                        'LineWidth',1.2,'marker','^');
                    hold on
                end
            end
            
            if solver_switch.geigpen
                if out_geigpen.success
                    index.geigpen = index.geigpen+1;
                    T.geigpen(index.geigpen) = out_geigpen.T;
                    pnrm2.geigpen(index.geigpen) = out_geigpen.pnrm2;
                    MRE.geigpen(index.geigpen) = out_geigpen.maxres;
                    trXAX.geigpen(index.geigpen) = out_geigpen.trXAX;
                    if isfield(out_geigpen,'evrelerr')
                        EREplot.geigpen = EREplot.geigpen + 1;
                        ERE.geigpen(index.geigpen) = out_geigpen.evrelerr;
                    end
                    iter.geigpen(index.geigpen,1) = out_geigpen.iter.inner;
                    iter.geigpen(index.geigpen,2) = out_geigpen.iter.outer;
                    figure(eigenValues);
                    ev = out_geigpen.ev;
                    plot(real(ev),'color',C(5,:),'LineStyle','none',...
                        'LineWidth',1.2,'marker','p');
                end
            end
            
            figure(eigenValues)
            
            legEntries = entries(find(structfun(structFind,index)));
            legend(legEntries,'location','best');
            legend('boxoff');
            xlabel('Index');
            ylabel('Value');
            title(['Eigenvalue Estimate Plot, Matrix Pencil ',...
                num2str(j),', tol = ',num2str(common_tol)]);
            set(gca,'fontsize',12); hold on; grid on;
            hold off; drawnow;
            if exist([BaseRoot,'Matrix_Pencil_',num2str(j),...
                    ', k = ',num2str(k),' tol = ',...
                    num2str(common_tol),...
                    ' ev estimate plot.fig'], 'file');
                delete([BaseRoot,'Matrix_Pencil_',num2str(j),...
                    ', k = ',num2str(k),' tol = ',...
                    num2str(common_tol),...
                    ' ev estimate plot.fig']);
            end
            savefig(eigenValues, [BaseRoot,'Matrix_Pencil_',num2str(j),...
                ', k = ',num2str(k),' tol = ',num2str(common_tol),...
                ' ev estimate plot.fig']);
            close(eigenValues);
            
        end
        
        if ~exist([BaseRoot,BMS.A_name,' Spy Plot.fig'],'file')
            A_spy = figure;
        spy(BMS.A);
        title(['Spy Plot of ', strrep(BMS.A_name,'_',' ')]);
        
%         if exist([BaseRoot,BMS.A_name,' Spy Plot.fig'],'file')
%             delete([BaseRoot,BMS.A_name,' Spy Plot.fig']);
%         end
        
        savefig(A_spy,[BaseRoot,BMS.A_name,' Spy Plot.fig']);
        close(A_spy);
        end
        
        
        if ~exist([BaseRoot,BMS.B_name,' Spy Plot.fig'],'file')
            
        B_spy = figure;
        spy(BMS.B);
        title(['Spy Plot of ', strrep(BMS.B_name,'_',' ')]);
%         if exist([BaseRoot,BMS.B_name,' Spy Plot.fig'],'file')
%             delete([BaseRoot,BMS.B_name,' Spy Plot.fig']);
%         end
        savefig(B_spy,[BaseRoot,BMS.B_name,' Spy Plot.fig']);
        close(B_spy);
        end
        if any(structfun(structFind,index))
            
            Time = figure;
            SupNorm = figure;
            PNRM2 = figure;
            projTR = figure;
            evRELERR = figure;
            
            if solver_switch.lobpcg
                LiterFIG = figure;
            end
            
            if solver_switch.geigpen
                GiterFIG = figure;
            end
            
            if solver_switch.eigs_sm
                if index.eigs_sm
                    figure(Time);
                    hold on
                    plot(K(1:index.eigs_sm),T.eigs_sm(1:index.eigs_sm),...
                        'color',C(1,:),'marker', 'd','LineWidth',1.2);
                    
                    figure(SupNorm);
                    hold on
                    plot(K(1:index.eigs_sm),...
                        MRE.eigs_sm(1:index.eigs_sm),'color', C(1,:),...
                        'marker','d','LineWidth',1.2);
                    
                    figure(PNRM2);
                    hold on
                    plot(K(1:index.eigs_sm),...
                        pnrm2.eigs_sm(1:index.eigs_sm),'color', C(1,:),...
                        'marker', 'd', 'LineWidth',1.2);
                    
                    figure(projTR);
                    hold on
                    plot(K(1:index.eigs_sm),...
                        trXAX.eigs_sm(1:index.eigs_sm),'color', C(1,:),...
                        'marker', 'd','LineWidth',1.2);
                end
            end
            
            if solver_switch.eigs_sa
                if index.eigs_sa
                    figure(Time);
                    hold on
                    plot(K(1:index.eigs_sa),T.eigs_sa(1:index.eigs_sa),...
                        'color', C(2,:), 'marker', 'o','LineWidth',1.2);
                    
                    figure(SupNorm);
                    hold on
                    plot(K(1:index.eigs_sa),...
                        MRE.eigs_sa(1:index.eigs_sa),'color', C(2,:),...
                        'marker', 'o','LineWidth',1.2);
                    
                    figure(PNRM2);
                    hold on
                    plot(K(1:index.eigs_sa),...
                        pnrm2.eigs_sa(1:index.eigs_sa),'color', C(2,:),...
                        'marker', 'o','LineWidth',1.2);
                    
                    figure(projTR);
                    hold on
                    plot(K(1:index.eigs_sa),...
                        trXAX.eigs_sa(1:index.eigs_sa),'color',C(2,:),...
                        'marker','o','LineWidth',1.2);
                    
                    figure(evRELERR)
                    hold on
                    plot(K(1:EREplot.eigs_sa),...
                        ERE.eigs_sa(1:EREplot.eigs_sa),'color',C(2,:),...
                        'marker','o','LineWidth',1.2);
                       
                end
            end
            % dont want to plot once 
            if solver_switch.bleigifp
                if index.bleigifp
                    figure(Time);
                    hold on
                    plot(K(1:index.bleigifp),...
                        T.bleigifp(1:index.bleigifp),'color',C(3,:),...
                        'marker','s','LineWidth',1.2);
                    
                    figure(SupNorm);
                    hold on
                    plot(K(1:index.bleigifp),...
                        MRE.bleigifp(1:index.bleigifp),'color',C(3,:),...
                        'marker','s','LineWidth',1.2);
                    
                    figure(PNRM2);
                    plot(K(1:index.bleigifp),...
                        pnrm2.bleigifp(1:index.bleigifp),'color',C(3,:),...
                        'marker','s','LineWidth',1.2);
                    
                    figure(projTR);
                    plot(K(1:index.bleigifp),...
                        trXAX.bleigifp(1:index.bleigifp),'color',C(3,:),...
                        'marker','s','LineWidth',1.2);
                    
                    figure(evRELERR)
                    hold on
                    plot(K(1:EREplot.bleigifp),...
                        ERE.bleigifp(1:EREplot.bleigifp),'color',C(3,:),...
                        'marker','s','LineWidth',1.2);
                end
            end
            
            if solver_switch.lobpcg
                if index.lobpcg
                    figure(Time);
                    hold on
                    plot(K(1:index.lobpcg),T.lobpcg(1:index.lobpcg),...
                        'color',C(4,:),'marker','^','LineWidth',1.2);
                    
                    figure(SupNorm);
                    hold on
                    plot(K(1:index.lobpcg),...
                        MRE.lobpcg(1:index.lobpcg),'color',C(4,:),...
                        'marker','^','LineWidth',1.2);
                    
                    figure(PNRM2);
                    hold on
                    plot(K(1:index.lobpcg),...
                        pnrm2.lobpcg(1:index.lobpcg),'color',C(4,:),...
                        'marker','^','LineWidth',1.2);
                    
                    figure(projTR);
                    hold on
                    plot(K(1:index.lobpcg),...
                        trXAX.lobpcg(1:index.lobpcg),'color',C(4,:),...
                        'marker','^','LineWidth',1.2);
                    
                    figure(evRELERR);
                    hold on
                    plot(K(1:EREplot.lobpcg),...
                        ERE.lobpcg(1:EREplot.lobpcg),'color',C(4,:),...
                        'marker','^','LineWidth',1.2);                   

                end
            end
            
            if solver_switch.geigpen
                if index.geigpen
                    figure(Time);
                    hold on
                    plot(K(1:index.geigpen),T.geigpen(1:index.geigpen),...
                        'color',C(5,:),'marker','p','LineWidth',1.2);
                    
                    figure(SupNorm);
                    hold on
                    plot(K(1:index.geigpen),...
                        MRE.geigpen(1:index.geigpen),'color',C(5,:),...
                        'marker','p','LineWidth',1.2);
                    
                    figure(PNRM2);
                    hold on
                    plot(K(1:index.geigpen),...
                        pnrm2.geigpen(1:index.geigpen),'color',C(5,:),...
                        'marker','p','LineWidth',1.2);
                    
                    figure(projTR);
                    hold on
                    plot(K(1:index.geigpen),...
                        trXAX.geigpen(1:index.geigpen),'color',C(5,:),...
                        'marker','p','LineWidth',1.2);
                    
                    figure(evRELERR);
                    hold on
                    plot(K(1:EREplot.geigpen),...
                        ERE.geigpen(1:EREplot.geigpen),'color',C(5,:),...
                        'marker','p','LineWidth',1.2);                    
                end
            end
            
            if solver_switch.lobpcg
                figure(LiterFIG);
                plot(K(1:index.lobpcg),iter.lobpcg(1:index.lobpcg),...
                    'color','r');
                legend('lobpcg iterCount');
                legend('boxoff');
                xlabel('Number of Eigenpairs');
                ylabel('Iteration Count');
                title(['lobpcg iteration count, Matrix Pencil ',...
                    num2str(j),', tol = ',num2str(common_tol)]);
                set(gca,'fontsize',12); hold on; grid on;
                hold off; drawnow;
                if exist([BaseRoot,'Matrix_Pencil_',num2str(j),...
                        ', K = [',num2str(K),'], tol = ',...
                        num2str(common_tol),...
                        ' lobpcg iteration record.fig'], 'file');
                    delete([BaseRoot,'Matrix_Pencil_',num2str(j),...
                        ', K = [',num2str(K),'], tol = ',...
                        num2str(common_tol),...
                        ' lobpcg iteration record.fig']);
                end
                savefig(LiterFIG, [BaseRoot,'Matrix_Pencil_',num2str(j),...
                    ', K = [',num2str(K),'], tol = ',...
                    num2str(common_tol),' lobpcg iteration record.fig']);
                close(LiterFIG);
            end
            
            if solver_switch.geigpen
                figure(GiterFIG);
                plotyy(K(1:index.geigpen),...
                    iter.geigpen(1:index.geigpen,1),...
                    K(1:index.geigpen),iter.geigpen(1:index.geigpen,2));
                legend('geigpen inner iterCount',...
                    'geigpen outer iterCount');
                legend('boxoff');
                xlabel('Number of Eigenpairs');
                ylabel('Iteration Count');
                title(['geigpen iteration count, Matrix Pencil ',...
                    num2str(j),', tol = ',num2str(common_tol)]);
                set(gca,'fontsize',12); hold on; grid on;
                hold off; drawnow;
                if exist([BaseRoot,'Matrix_Pencil_',num2str(j),...
                        ', K = [',num2str(K),'], tol = ',...
                        num2str(common_tol),...
                        ' geigpen iteration record.fig'], 'file');
                    delete([BaseRoot,'Matrix_Pencil_',num2str(j),...
                        ', K = [',num2str(K),'], tol = ',...
                        num2str(common_tol),...
                        ' geigpen iteration record.fig']);
                end
                savefig(GiterFIG, [BaseRoot,'Matrix_Pencil_',...
                    num2str(j),', K = [',num2str(K),'], tol = ',...
                    num2str(common_tol),' geigpen iteration record.fig']);
                
                close(GiterFIG);
            end

            
            legEntries = entries(find(structfun(structFind,index)));

            figure(Time);

            legend(legEntries,'location','best');
            legend('boxoff');
            xlabel('Number of Eigenpairs');
            ylabel('Running Time');
            title(['Comparison of Times to Solution, Matrix Pencil ',...
                num2str(j),', tol = ',num2str(common_tol)]);
            set(gca,'fontsize',12); hold on; grid on;
            hold off; drawnow;
            if exist([BaseRoot,'Matrix_Pencil_',num2str(j),...
                    ', K = [',num2str(K),'], tol = ',...
                    num2str(common_tol),...
                    ' Eigenpair vs Running Time Comparison.fig'], 'file');
                delete([BaseRoot,'Matrix_Pencil_',num2str(j),...
                    ', K = [',num2str(K),'], tol = ',...
                    num2str(common_tol),...
                    ' Eigenpair vs Running Time Comparison.fig']);
            end
            savefig(Time, [BaseRoot,'Matrix_Pencil_',num2str(j),...
                ', K = [',num2str(K),'], tol = ',num2str(common_tol),...
                ' Eigenpair vs Running Time Comparison.fig']);
            
            close(Time);
            
            
            figure(SupNorm);           
            legend(legEntries,'location','best');
            legend('boxoff');
            xlabel('Number of Eigenpairs');
            ylabel('Semilog of Sup Norm Relative Residual Error')
            title(['Comparison of Solver Accuracies, Matrix Pencil ',...
                num2str(j),', tol = ',num2str(common_tol)]);
            set(gca,'yscale','log','fontsize',12); hold on; grid on;
            hold off; drawnow;
            if exist([BaseRoot,'Matrix_Pencil_',num2str(j),...
                    ', K = [',num2str(K),'], tol = ',...
                    num2str(common_tol),...
                    ' Eigenpair vs MRE Comparison.fig'], 'file');
                delete([BaseRoot,'Matrix_Pencil_',num2str(j),...
                    ', K = [',num2str(K),'], tol = ',...
                    num2str(common_tol),...
                    ' Eigenpair vs MRE Comparison.fig']);
            end
            savefig(SupNorm, [BaseRoot,'Matrix_Pencil_',num2str(j),...
                ', K = [',num2str(K),'], tol = ',num2str(common_tol),...
                ' Eigenpair vs MRE Comparison.fig']);
            
            close(SupNorm);
            
            
            figure(PNRM2);
%             legEntries
%             pause;
            legend(legEntries,'location','best');
            legend('boxoff');
            xlabel('Number of Eigenpairs');
            ylabel('Semilog of B-Orthogonality');
            title(['Comparison of Eigenvalue Solvers, Matrix Pencil ',...
                num2str(j),', tol = ',num2str(common_tol)]);
            set(gca,'yscale','log','fontsize',12); hold on; grid on;
            hold off; drawnow;
            if exist([BaseRoot,'Matrix_Pencil_',num2str(j),...
                    ', K = [',num2str(K),'], tol = ',...
                    num2str(common_tol),...
                    ' Eigenpair vs pnrm2 Comparison.fig'], 'file');
                delete([BaseRoot,'Matrix_Pencil_',num2str(j),...
                    ', K = [',num2str(K),'], tol = ',...
                    num2str(common_tol),...
                    ' Eigenpair vs pnrm2 Comparison.fig']);
            end
            savefig(PNRM2, [BaseRoot,'Matrix_Pencil_',num2str(j),...
                ', K = [',num2str(K),'], tol = ',num2str(common_tol),...
                ' Eigenpair vs pnrm2 Comparison.fig']);
            
            close(PNRM2);
            
            
            figure(projTR);
%             legEntries
%             pause;
            legend(legEntries,'location','best')
            legend('boxoff');
            xlabel('Number of Eigenpairs');
            ylabel('Semilog of Projection Trace: tr(X''AX)');
            title(['Comparison of Eigenvalue Solvers, Matrix Pencil, ',...
                num2str(j),', tol = ',num2str(common_tol)]);
            set(gca,'yscale','log','fontsize',12); hold on; grid on;
            hold off; drawnow;
            if exist([BaseRoot,'Matrix_Pencil_',num2str(j),...
                    ', K = [',num2str(K),'], tol = ',...
                    num2str(common_tol),...
                    ' Eigenpair vs projTR Comparison.fig'], 'file');
                delete([BaseRoot,'Matrix_Pencil_',num2str(j),...
                    ', K = [',num2str(K),'], tol = ',...
                    num2str(common_tol),...
                    ' Eigenpair vs projTR Comparison.fig']);
            end
            savefig(projTR, [BaseRoot,'Matrix_Pencil_',num2str(j),...
                ', K = [',num2str(K),'], tol = ',...
                num2str(common_tol),...
                ' Eigenpair vs projTR Comparison.fig']);
            
            close(projTR);
            
            figure(evRELERR);
%                         legEntries(2:end)
%                         pause;
            legend(legEntries(2:end),'location','best')
            legend('boxoff');
            xlabel('Number of Eigenpairs');
            ylabel('Semilog of ev Relative Error');
            title(['Comparison of Eigenvalue Solvers, Matrix Pencil, ',...
                num2str(j),', tol = ',num2str(common_tol)]);
            set(gca,'yscale','log','fontsize',12); hold on; grid on;
            hold off; drawnow;
            if exist([BaseRoot,'Matrix_Pencil_',num2str(j),...
                    ', K = [',num2str(K),'], tol = ',...
                    num2str(common_tol),...
                    ' Eigenpair vs evRELERR Comparison.fig'], 'file');
                delete([BaseRoot,'Matrix_Pencil_',num2str(j),...
                    ', K = [',num2str(K),'], tol = ',...
                    num2str(common_tol),...
                    ' Eigenpair vs evRELERR Comparison.fig']);
            end
            savefig(evRELERR, [BaseRoot,'Matrix_Pencil_',num2str(j),...
                ', K = [',num2str(K),'], tol = ',...
                num2str(common_tol),...
                ' Eigenpair vs evRELERR Comparison.fig']);
            
            close(evRELERR);
            clear legEntries;
            diary off
            save([BaseRoot,'Matrix_Pencil_',num2str(j),' Overall V, K',...
                    ' = [',num2str(K),'], tol = ',...
                    num2str(common_tol),'.mat'],'solver_switch','K',...
                    'common_tol','index','T', 'MRE', 'trXAX', 'pnrm2',...
                    'ERE','EREplot');
        end
        
        J = J + 1;
        
        if J>numMP_index
            break;
        else
            j = MP_index(J);
            if j>MaxMPI
                break;
            end
        end
              
    end
end

