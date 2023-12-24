function [out_eigs_sm, out_eigs_sa, out_bleigifp, out_lobpcg,...
                out_geigpen] = RUN_COMPARISON(BMS,k,opts_eigs_sm,...
                opts_eigs_sa,opts_bleigifp,opts_lobpcg,opts_geigpen,...
                solver_switch,BaseRoot,FunctionDir,entries,structFind)
            
frres = @(R,r)sqrt(sum(R.^2))./max(1,sqrt(sum(r.^2)));
rng('default');

% cd(BaseRoot);

if exist([BaseRoot,'k = ',num2str(k),' tol = ',...
        num2str(opts_geigpen.tol2),'.mat'],'file');
    delete([BaseRoot,'k = ',num2str(k),' tol = ',...
        num2str(opts_geigpen.tol2),'.mat']);
end

save([BaseRoot,'k = ',num2str(k),' tol = ',num2str(opts_geigpen.tol2),...
    '.mat'],'opts_eigs_sm','opts_eigs_sa','opts_bleigifp','opts_lobpcg',...
    'opts_geigpen');

% cd(FunctionDir);

fprintf('======================================================== \n');
% 
% fprintf(['\n*** Solving Generalized Eigenvalue Problem ***\n',...
%     '%s : (%s,%s)\n{n = %i | sprsty A,B = %f, %f | k = %i | ',...
%     'tol = %12.8e}\n'],BMS.Name,BMS.A_name,BMS.B_name,BMS.A_size(2),...
%     BMS.A_sprsty, BMS.B_sprsty,k,opts_eigs_sm.tol);

fprintf(['\n*** Solving Generalized Eigenvalue Problem ***\n',...
    '%s : (%s,%s)\n{n = %i | sprsty A,B = %f, %f | k = %i | ',...
    'tol = %12.8e}\n'],BMS.Name,BMS.A_name,BMS.B_name,BMS.A_size(2),...
    nnz(BMS.A)/(BMS.A_size(2)^2), nnz(BMS.B)/(BMS.B_size(2)^2),k,opts_eigs_sm.tol);
if solver_switch.eigs_sm
    
    fprintf(['\n----- Running EIGS_SM *****',...
        ' EIGS_SM Display Output ----- \n\n']);
    
    if opts_eigs_sm.proceed == 1
        try
            RUN_EIGS_SM;
            out_eigs_sm.success = 1;
        catch ME_eigs_sm;
            warning('Error for RUN_EIGS_SA thrown.');
            disp(ME_eigs_sm);
            out_eigs_sm = [];
            out_eigs_sm.success = 0;
        end
    else
        warning('Forbidden to run EIGS_SA by MAIN_SCRIPT.');
        out_eigs_sm = [];
        out_eigs_sm.success = 0;
    end
else
    out_eigs_sm = [];
    out_eigs_sm.success = 0;
    
end

if solver_switch.eigs_sa
    
    fprintf(['\n----- Running EIGS_SA *****',...
        ' EIGS_SA Display Output ----- \n\n']);
    
    if opts_eigs_sa.proceed == 1
        try
            RUN_EIGS_SA;
            out_eigs_sa.success = 1;
        catch ME_eigs_sa;            
            warning('Error for RUN_EIGS_SA thrown.');
            disp(ME_eigs_sa);
            out_eigs_sa = [];
            out_eigs_sa.success = 0;
        end
    else
        warning('Forbidden to run EIGS_SA by MAIN_SCRIPT.');
        out_eigs_sa = [];
        out_eigs_sa.success = 0;
    end
else    
    out_eigs_sa = [];
    out_eigs_sa.success = 0;
    
end

if solver_switch.bleigifp
    
    fprintf(['\n----- Running BLEIGIFP ***',...
        ' BLEIGIFP Display Output ----- \n\n']);
    if opts_bleigifp.proceed == 1
        try
            RUN_BLEIGIFP;
            out_bleigifp.success = 1;
        catch ME_bleigifp;            
            warning('Error for BLEIGIFP thrown.');
            disp(ME_bleigifp);
            out_bleigifp = [];
            out_bleigifp.success = 0;
        end
    else
        warning('Forbidden to run BLEIGPIFP by MAIN_SCRIPT.');
        out_bleigifp = [];
        out_bleigifp.success = 0;
    end
else
    out_bleigifp = [];
    out_bleigifp.success = 0;    
end

if solver_switch.lobpcg
    
    fprintf(['\n----- Running LOBPCG *******',...
        ' LOBPCG Display Output ----- \n\n']);
    
    if opts_lobpcg.proceed == 1
        try
            RUN_LOGPCG;
            out_lobpcg.success = 1;
        catch ME_lobpcg;
            warning('Error for LOBPCG thrown.');
            disp(ME_lobpcg);
            out_lobpcg = [];
            out_lobpcg.success = 0;
        end
    else
        warning('Forbidden to run LOBPCG by MAIN_SCRIPT.');
        out_lobpcg = [];
        out_lobpcg.success = 0;
    end
else
    out_lobpcg = [];
    out_lobpcg.success = 0;    
end

if solver_switch.geigpen
    
    fprintf(['\n----- Running GEIGPEN *****',...
        ' GEIGPEN Display Output ----- \n\n']);
    
    if opts_geigpen.proceed == 1
        try
            RUN_GEIGPEN;
            out_geigpen.success = 1;
        catch ME_geigpen;
            warning('Error for GEIGPEN thrown.');
            disp(ME_geigpen);
            out_geigpen = [];
            out_geigpen.success = 0;
        end
    else
        warning('Forbidden to run GEIGPEN by MAIN_SCRIPT.');
        out_geigpen = [];
        out_geigpen.success = 0;
    end
else
    out_geigpen = [];
    out_geigpen.success = 0;    
end


clear opts_eigs_sm;
clear opts_eigs_sa;
clear opts_bleigifp;
clear opts_lobpcg;
clear BMS;

solverCount = sum(structfun(structFind,solver_switch));
check = zeros(solverCount,1);

check(1) = out_eigs_sm.success;
check(2) = out_eigs_sa.success;
check(3) = out_bleigifp.success;
check(4) = out_lobpcg.success;
check(5) = out_geigpen.success;

if any(check)
    
   fprintf('\nTTS:\n');
   if check(1)
       fprintf('eigs-SM:  %12.8e\n',out_eigs_sm.T(3));
   end   
   if check(2)
       fprintf('eigs-SA:  %12.8e\n',out_eigs_sa.T);
   end
   if check(3)
       fprintf('bleigifp: %12.8e\n',out_bleigifp.T);
   end
   if check(4)
       fprintf('lobpcg:   %12.8e\n',out_lobpcg.T);
   end
   if check(5)
       fprintf('geigpen:  %12.8e\n',out_geigpen.T);
   end
 
    
   fprintf('\ntrXAX:\n');
   if check(1)
       fprintf('eigs-SM:  %12.8e\n',out_eigs_sm.trXAX);
   end
   if check(2)
       fprintf('eigs-SA:  %12.8e\n',out_eigs_sa.trXAX);
   end
   if check(3)
       fprintf('bleigifp: %12.8e\n',out_bleigifp.trXAX);
   end
   if check(4)
       fprintf('lobpcg:   %12.8e\n',out_lobpcg.trXAX);
   end
   if check(5)
       fprintf('geigpen:  %12.8e\n',out_geigpen.trXAX);
   end
   
   fprintf('\nB-Orthogonality:\n');
   
   if check(1)
       fprintf('eigs-SM:  %12.8e\n',out_eigs_sm.pnrm2);
   end
   if check(2)
       fprintf('eigs-SA:  %12.8e\n',out_eigs_sa.pnrm2);
   end
   if check(3)
       fprintf('bleigifp: %12.8e\n',out_bleigifp.pnrm2);
   end
   if check(4)
       fprintf('lobpcg:   %12.8e\n',out_lobpcg.pnrm2);
   end
   if check(5)
       fprintf('geigpen:  %12.8e\n',out_geigpen.pnrm2);
   end
   
   fprintf('\nSup Norm Relative Residual:\n');

   if check(1)
       fprintf('eigs-SM:  %12.8e\n',out_eigs_sm.maxres);
   end   
   if check(2)
       fprintf('eigs-SA:  %12.8e\n',out_eigs_sa.maxres);
   end
   if check(3)
       fprintf('bleigifp: %12.8e\n',out_bleigifp.maxres);
   end
   if check(4)
       fprintf('lobpcg:   %12.8e\n',out_lobpcg.maxres);
   end
   if check(5)
       fprintf('geigpen:  %12.8e\n',out_geigpen.maxres);
   end
   
end

check(1) = 0;

if any(check)
    fprintf('\nIteration Count:\n');
    
    if check(4)
        if isfield(out_lobpcg,'iter')
            fprintf('lobpcg:   %d\n',out_lobpcg.iter);
                
        end                
    end

    
    if check(5)
        if isfield(out_geigpen,'iter')
            fprintf('geigpen inner/outer:  %d/%d\n',...
                out_geigpen.iter.inner,out_geigpen.iter.outer);
        end                
    end
end

check(2:3) = 0;

%cd(BaseRoot);
save([BaseRoot,'k = ',num2str(k),' tol = ',num2str(opts_geigpen.tol2),...
    '.mat'],'out_eigs_sm','out_eigs_sa','out_bleigifp', 'out_lobpcg',...
    'out_geigpen','-append');
fprintf('\n');

clear k;
clear opts_geigpen;

%% The function definining the EIGS_SM computation

    function RUN_EIGS_SM
        
        out_eigs_sm = [];
        

        if isfield(BMS,'s')
            s = -min(0,BMS.s);
            t1 = BMS.t;
        else
            s = 0;
            t1 = 0;
        end
        
        tic;
        [X_eigs_sm,D_eigs_sm] = eigs(BMS.A+s*BMS.B,BMS.B,k,'SM',...
            opts_eigs_sm);
        D_eigs_sm = D_eigs_sm - s*speye(k);
        t2 = toc;
        
        t = t1+t2;
        
        AX_eigs_sm = BMS.A*X_eigs_sm; BX_eigs_sm = BMS.B*X_eigs_sm;
        trXAX_eigs_sm = dot(X_eigs_sm(:),AX_eigs_sm(:));
        pnrm2_eigs_sm = norm(X_eigs_sm'*BX_eigs_sm-speye(k),'fro');
        RES_eigs_sm = frres(AX_eigs_sm-BX_eigs_sm*D_eigs_sm,AX_eigs_sm);
        maxres_eigs_sm = max(RES_eigs_sm);
        ev_eigs_sm = sort(diag(D_eigs_sm));
        

        out_eigs_sm.s = s;
        clear s;
        out_eigs_sm.T = [t1, t2, t]';
        clear t1 t2 t;
%         out_eigs_sm.X = X_eigs_sm;
        clear X_eigs_sm;
%         out_eigs_sm.D = sparse(D_eigs_sm);
        clear D_eigs_sm;
%         out_eigs_sm.AX = AX_eigs_sm;
        clear AX_eigs_sm;
%         out_eigs_sm.BX = BX_eigs_sm;
        clear BX_eigs_sm;
        out_eigs_sm.trXAX = trXAX_eigs_sm;
        clear trXAX_eigs_sm;
        out_eigs_sm.pnrm2 = pnrm2_eigs_sm;
        clear pnrm2_eigs_sm;
%         out_eigs_sm.RES = RES_eigs_sm;
        clear RES_eigs_sm;
        out_eigs_sm.maxres = maxres_eigs_sm;
        clear maxres_eigs_sm;        
        out_eigs_sm.ev = ev_eigs_sm;
        clear ev_eigs_sm;
        
    end
%% The function definining the EIGS_SA computation

    function RUN_EIGS_SA
        
        out_eigs_sa = [];
        
        tic;
        [X_eigs_sa,D_eigs_sa] = eigs(BMS.A,BMS.B,k,'SA',...
            opts_eigs_sa);
        t = toc;
        
        AX_eigs_sa = BMS.A*X_eigs_sa; BX_eigs_sa = BMS.B*X_eigs_sa;
        trXAX_eigs_sa = dot(X_eigs_sa(:),AX_eigs_sa(:));
        pnrm2_eigs_sa = norm(X_eigs_sa'*BX_eigs_sa-speye(k),'fro');
        RES_eigs_sa = frres(AX_eigs_sa-BX_eigs_sa*D_eigs_sa,AX_eigs_sa);
        maxres_eigs_sa = max(RES_eigs_sa);
        ev_eigs_sa = sort(diag(D_eigs_sa));
        if exist('out_eigs_sm','var')
            if isfield(out_eigs_sm,'ev')
                if any(out_eigs_sm.ev)&&~any(isnan(out_eigs_sm.ev))
            evrelerr_eigs_sa = norm(out_eigs_sm.ev-ev_eigs_sa)...
                /norm(out_eigs_sm.ev);
                end                
            end
        end
        out_eigs_sa.T = t;
        clear t
%         out_eigs_sa.X = X_eigs_sa;
        clear X_eigs_sa;
%         out_eigs_sa.D = sparse(D_eigs_sa);
        clear D_eigs_sa;
%         out_eigs_sa.AX = AX_eigs_sa;
        clear AX_eigs_sa;
%         out_eigs_sa.BX = BX_eigs_sa;
        clear BX_eigs_sa;
         out_eigs_sa.trXAX = trXAX_eigs_sa;
        clear trXAX_eigs_sa;
        out_eigs_sa.pnrm2 = pnrm2_eigs_sa;
        clear pnrm2_eigs_sa;
%         out_eigs_sa.RES = RES_eigs_sa;
        clear RES_eigs_sa;
        out_eigs_sa.maxres = maxres_eigs_sa;
        clear maxres_eigs_sa;        
        out_eigs_sa.ev = ev_eigs_sa;
        clear ev_eigs_sa;
        out_eigs_sa.evrelerr = evrelerr_eigs_sa;
        clear evrelerr_eigs_sa;
        
    end

%%

    function RUN_BLEIGIFP
        
        tic;
        [d_bleigifp,X_bleigifp] = bleigifp(BMS.A,BMS.B,k,opts_bleigifp);
        k1_bleigifp = min(k,numel(d_bleigifp));
        d_bleigifp = d_bleigifp(1:k1_bleigifp);
        X_bleigifp = X_bleigifp(:,1:k1_bleigifp);
        t = toc;
        
        if k1_bleigifp < k, NoMss_bleigifp = k-k1_bleigifp; end
        
        fprintf('Found %d/%d eigenvalues.\n',k1_bleigifp,k); 
        if exist('NoMss_bleigifp','var')
            fprintf('%d eigenvalues are missing.\n',NoMss_bleigifp);
        else
            fprintf('All eigenvalues converged.\n');
        end
        
        AX_bleigifp = BMS.A*X_bleigifp; BX_bleigifp = BMS.B*X_bleigifp;
        trXAX_bleigifp = dot(X_bleigifp(:),AX_bleigifp(:));
        pnrm2_bleigifp = ...
            norm(X_bleigifp'*BX_bleigifp-speye(k1_bleigifp),'fro');
        RES_bleigifp = ...
            frres(AX_bleigifp-BX_bleigifp*diag(d_bleigifp),AX_bleigifp);
        maxres_bleigifp = max(RES_bleigifp);
        ev_bleigifp = sort(d_bleigifp);
        if exist('out_eigs_sm','var')
            if isfield(out_eigs_sm,'ev')
                if any(out_eigs_sm.ev)&&~any(isnan(out_eigs_sm.ev))
                    evrelerr_bleigifp =...
                        norm(out_eigs_sm.ev(1:k1_bleigifp)...
                        -ev_bleigifp(1:k1_bleigifp))...
                        /norm(out_eigs_sm.ev(1:k1_bleigifp));
                end
            end
        end
                   
        out_bleigifp.T = t;
        clear t;
        %         out_bleigifp.d = d_bleigifp;
        clear d_bleigifp;
        %         out_bleigifp.X = X_bleigifp;
        clear X_bleigifp;
        %         out_bleigifp.AX = AX_bleigifp;
        clear AX_bleigifp;
        %         out_bleigifp.BX = BX_bleigifp;
        clear BX_bleigifp;
        out_bleigifp.k1 = k1_bleigifp;
        clear k1_bleigifp;
        if exist('NoMss_bleigifp','var');
            out_bleigifp.NoMss = NoMss_bleigifp;
            clear NoMss_bleigifp;
        end
        out_bleigifp.trXAX = trXAX_bleigifp;
        clear trXAX_bleigifp;
        out_bleigifp.pnrm2 = pnrm2_bleigifp;
        clear pnrm2_bleigifp;
        %         out_bleigifp.RES = RES_bleigifp;
        clear RES_bleigifp;
        out_bleigifp.maxres = maxres_bleigifp;
        clear maxres_bleigifp;
        out_bleigifp.ev = ev_bleigifp;
        clear ev_bleigifp;
        out_bleigifp.evrelerr = evrelerr_bleigifp;
        clear evrelerr_bleigifp;
        
        
    end

%%

    function RUN_LOGPCG
        
        tic;
        opts_lobpcg.blockVectorX = randn(BMS.A_size(2),k);
        [X_lobpcg,d_lobpcg, failureFlag_lobpcg,lambdaHistory_lobpcg,...
            residualNormsHistory_lobpcg] = ...
            lobpcg(opts_lobpcg.blockVectorX,BMS.A,BMS.B,...
            opts_lobpcg.residualTolerance,opts_lobpcg.maxIterations,...
            opts_lobpcg.verbosityLevel);
        t = toc;
        
        X_lobpcg = X_lobpcg(:,1:k); d_lobpcg = d_lobpcg(1:k);
        AX_lobpcg = BMS.A*X_lobpcg; BX_lobpcg = BMS.B*X_lobpcg;
        trXAX_lobpcg = dot(X_lobpcg(:),AX_lobpcg(:));
        pnrm2_lobpcg = norm(X_lobpcg'*BX_lobpcg-speye(k),'fro');
        RES_lobpcg = frres(AX_lobpcg-BX_lobpcg*diag(d_lobpcg),AX_lobpcg);
        maxres_lobpcg = max(RES_lobpcg);
        ev_lobpcg = sort(d_lobpcg);       
        if exist('out_eigs_sm','var')
            if isfield(out_eigs_sm,'ev')
                if any(out_eigs_sm.ev)&&~any(isnan(out_eigs_sm.ev))
                    evrelerr_lobpcg = norm(out_eigs_sm.ev-ev_lobpcg)...
                        /norm(out_eigs_sm.ev);
                end
            end
        end
        
        out_lobpcg.iter = size(lambdaHistory_lobpcg,2);
        %         out_lobpcg.X = X_lobpcg;
        clear X_lobpcg;
        %         out_lobpcg.d = d_lobpcg;
        clear d_lobpcg;
        out_lobpcg.failureFlag_lobpcg = ...
            failureFlag_lobpcg;
        clear failureFlag_lobpcg;
                out_lobpcg.lambdaHistory_lobpcg = ...
                    lambdaHistory_lobpcg;
        clear lambdaHistory_lobpcg;
                out_lobpcg.residualNormsHistory_lobpcg = ...
                    residualNormsHistory_lobpcg;
        clear residualNormsHistory_lobpcg;
        out_lobpcg.T = t;
        clear t;
        %         out_lobpcg.AX = AX_lobpcg;
        clear AX_lobpcg;
        %         out_lobpcg.BX = BX_lobpcg;
        clear BX_lobpcg;
        out_lobpcg.trXAX = trXAX_lobpcg;
        clear trXAX_lobpcg;
        out_lobpcg.pnrm2 = pnrm2_lobpcg;
        clear pnrm2_lobpcg;
        %         out_lobpcg.RES = RES_lobpcg;
        clear RES_lobpcg;
        out_lobpcg.maxres = maxres_lobpcg;
        clear maxres_lobpcg;
        out_lobpcg.ev = ev_lobpcg;
        clear ev_lobpcg;
        out_lobpcg.evrelerr = evrelerr_lobpcg;
        clear evrelerr_lobpcg;
        
    end

%%

    function RUN_GEIGPEN
%         tic;
%         s = eigs(BMS.A,BMS.B,1,'SA',opts_eigs_sm); s = -min(0,s);
%         t1 = toc;
%         
%         tic;
%         [X_eigs_sm,D_eigs_sm] = eigs(BMS.A+s*BMS.B,BMS.B,k,'SM',...
%             opts_eigs_sm);
%         D_eigs_sm = D_eigs_sm - s*speye(k);
%         t2 = toc;
%         
%         t = t1+t2;
        
%         if isfield(BMS,'s')
%             s = -min(0,BMS.s);
%         else
%             s = 0;
%             BMS.t = 0;
%         end
%         
        tic;
        [X_geigpen,D_geigpen,out_geigpen] = geigpen(BMS.A,BMS.B,...
            BMS.A_size(2),k,[],opts_geigpen);
                t = toc;
%         [X_geigpen,D_geigpen,out_geigpen] = geigpen(BMS.A + s*BMS.B,...
%             BMS.B,BMS.A_size(2),k,[],opts_geigpen);
%         D_geigpen = D_geigpen - s*speye(k);

%         t = BMS.t + t;
        
        AX_geigpen = BMS.A*X_geigpen; BX_geigpen = BMS.B*X_geigpen;
        trXAX_geigpen = dot(X_geigpen(:),AX_geigpen(:));
        pnrm2_geigpen = norm(X_geigpen'*BX_geigpen-speye(k),'fro');
        out_geigpen.ev = sort(diag(D_geigpen));
        ev_geigpen  = out_geigpen.ev;
        RES_geigpen = frres(AX_geigpen-BX_geigpen*D_geigpen,AX_geigpen);
        maxres_geigpen = max(RES_geigpen);
        if exist('out_eigs_sm','var')
            if isfield(out_eigs_sm,'ev')
                if any(out_eigs_sm.ev)&&~any(isnan(out_eigs_sm.ev))
                    evrelerr_geigpen = norm(out_eigs_sm.ev-ev_geigpen)...
                        /norm(out_eigs_sm.ev);
                end
            end
        end
        geigpen_times = out_geigpen.geigpen_times;
        TimesTable = struct2table(geigpen_times);
        disp(TimesTable)
        %         out_geigpen.X = X_geigpen;
        clear X_geigpen;
        %         out_geigpen.D = D_geigpen;
        clear D_geigpen;
        out_geigpen.TimesTable = TimesTable;
        out_geigpen.T = t;
        clear t;
        %         out_geigpen.AX = AX_geigpen;
        clear AX_geigpen;
        %         out_geigpen.BX = BX_geigpen;
        clear BX_geigpen;
        out_geigpen.trXAX = trXAX_geigpen;
        clear trXAX_geigpen;
        out_geigpen.pnrm2 = pnrm2_geigpen;
        clear pnrm2_geigpen;
        out_geigpen.ev = ev_geigpen;
        clear ev_geigpen;
        %         out_geigpen.RES = RES_geigpen;
        clear RES_geigpen;
        out_geigpen.maxres = maxres_geigpen;
        clear maxres_geigpen;
        out_geigpen.evrelerr = evrelerr_geigpen;
        clear evrelerr_geigpen;
    end

end





% if out_eigs_sa.success
%     fprintf(['\neigs  tr(XtAX): %12.8e  B-orth: %5.2e  ',...
%         'maxres: %5.2e\n'],...
%         out_eigs_sa.trXAX,out_eigs_sa.pnrm2,out_eigs_sa.maxres);
% end

% if out_bleigifp.success
%     fprintf(['\nbleigifp tr(XtAX): %12.8e  B-orth: %5.2e  ',...
%         'maxres: %5.2e\n'],...
%         out_bleigifp.trXAX,out_bleigifp.pnrm2,out_bleigifp.maxres);
%     
%     if out_eigs_sa.success
%         if ~any(isnan(out_eigs_sa.ev))
%             fprintf('e-value rel-err: %9.4e\n\n',...
%                 norm(out_eigs_sa.ev(1:out_bleigifp.k1)...
%                 -out_bleigifp.ev(1:out_bleigifp.k1))...
%                 /norm(out_eigs_sa.ev(1:out_bleigifp.k1)));
%         end
%         
%     end
% end

% if out_lobpcg.success
%     fprintf('\nlobpcg tr(XtAX): %12.8e  B-orth: %5.2e  maxres: %5.2e\n',...
%         out_lobpcg.trXAX,out_lobpcg.pnrm2,out_lobpcg.maxres);
%     if out_eigs_sa.success
%         if ~any(isnan(out_eigs_sa.ev))
%             fprintf('e-value rel-err: %9.4e\n\n',...
%                 norm(out_eigs_sa.ev-out_lobpcg.ev)/norm(out_eigs_sa.ev));
%         end
%         
%     end
% end
% 
% if out_geigpen.success
%     fprintf('\ngeigpen tr(XtAX): %12.8e  B-orth: %5.2e  maxres: %5.2e\n',...
%         out_geigpen.trXAX,out_geigpen.pnrm2,out_geigpen.maxres);
%     if out_eigs_sa.success
%         if ~any(isnan(out_eigs_sa.ev))
%             fprintf('e-value rel-err: %9.4e\n\n',...
%                 norm(out_eigs_sa.ev-out_geigpen.ev)/norm(out_eigs_sa.ev));
%         end
%         
%     end
% end
