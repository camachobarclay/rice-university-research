function [X, D, out] = geigpen(A, B, n, k, X, opts)

% Compute the "smallest" k-dimensional eigenspace
% of a symmetric matrix A to a moderate accuracy
%
% Input:
%        A -- (n x n) symmetric matrix
%        B -- (n x n) symmetric matrix
%        n -- (1 x 1) matrix dimension
%        k -- (1 x 1) dimension of invariant subspace
%        X -- (n x k) initial matrix; could be empty
% Output:
%        X -- (n x k) orthonormal, approximate solution
%        D -- (k x k) sparse diagonal, holding eigenvalues
%
% Copyright@ Yin Zhang and co-authors, 2016
% Last modified: Aug. 15, 2016

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 3;    k = round(0.01*size(A,1));           end % Check
if nargin < 4;    X = [];                              end % Check
if nargin < 5;    opts = [];                           end % Check

% setup
symm = @(Z) 0.5*(Z + Z'); % Check
%Ik = eye(kw);

tic;
initialize; % kw, X, AX, BX, ... % Check
geigen_times.initialize = toc;

p = opts.arrp; % Check
proj_def = opts.proj; % Check
tol1sp = opts.tol1sp; % Check
maxouter = opts.maxouter; % Check

iterCount = [];
iterCount.inner = 0;
iterCount.outer = 0;
if opts.print, fprintf('ARR order p = %i\n',p); end % Check
if p>0, ip=symamd(B); L=chol(B(ip,ip),'lower'); end % Check
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
geigpen_times.dot = 0;
geigpen_times.XtBX_I = 0;
geigpen_times.G = 0;
geigpen_times.alfSpike = 0;
geigpen_times.backtrack = 0;
geigpen_times.rGnrm = 0;
geigpen_times.gpen = 0;
% geigpen_times.AX = 0;
% geigpen_times.BX = 0;
% geigpen_times.qr = 0;
% geigpen_times.XtAX = 0;
% geigpen_times.eig = 0;
% geigpen_times.XY = 0;
% geigpen_times.AXY = 0;
% geigpen_times.BXY = 0;
% geigpen_times.Res = 0; 

for outer = 1:maxouter % Check
    t0 = tic; % Check
    opts.outer = outer; % Check
    opts.proj = proj_def; % Check
    if outer <= 1, opts.proj = 0; end % Check
    
    
    %% =========== call solver ===========
    tic;
    [X,G,out] = gpen_solve(A,B,X,AX,BX,opts); % Check
    geigpen_times.gpen = geigpen_times.gpen + toc;
    
    gpen_times = out.gpen_times;
    geigpen_times.AX = geigpen_times.AX + gpen_times.AX;
    geigpen_times.BX = geigpen_times.BX + gpen_times.BX;
    geigpen_times.dot = geigpen_times.dot + gpen_times.dot;
    geigpen_times.XtBX_I = geigpen_times.XtBX_I + gpen_times.XtBX_I;
    geigpen_times.G = geigpen_times.G + gpen_times.G;
    geigpen_times.alfSpike = geigpen_times.alfSpike + gpen_times.alfSpike;
    geigpen_times.backtrack = geigpen_times.backtrack +...
        gpen_times.backtrack;
    geigpen_times.rGnrm = geigpen_times.rGnrm + gpen_times.rGnrm;
    
    %% =========== call ARR ==============
    G = [X G]; % Check
    for j = 1:p % Check
        X = A*X;  X(ip,:) = L'\(L\X(ip,:)); % Check
        G = [G X]; %#ok<*AGROW>  % Check
    end % Check

    tic;
    [X,D,ev,Res,AX,BX,gRR_times] = gRR(A,B,G);
    geigpen_times.gRR = geigpen_times.gRR + toc; % Check
    
    geigpen_times.AX = geigpen_times.AX + gRR_times.AX;
    geigpen_times.BX = geigpen_times.BX + gRR_times.BX;
    geigpen_times.qr = geigpen_times.qr + gRR_times.qr;
    geigpen_times.XtAX = geigpen_times.XtAX + gRR_times.XtAX;
    geigpen_times.XtBX = geigpen_times.XtBX + gRR_times.XtBX;
    geigpen_times.eig = geigpen_times.eig + gRR_times.eig;
    geigpen_times.XY = geigpen_times.XY + gRR_times.XY;
    geigpen_times.AXY = geigpen_times.AXY + gRR_times.AXY;
    geigpen_times.BXY = geigpen_times.BXY + gRR_times.BXY;
    geigpen_times.Res = geigpen_times.Res + gRR_times.Res;

    %% ==================================
    maxres = norm(Res(1:k),inf); % Check
    X = X(:,1:kw);  %#ok<SHVAI> % Check
    AX = AX(:,1:kw); % Check
    BX = BX(:,1:kw); % Check
    if opts.print % Check
        %XtBX_I = X'*BX-Ik; % Check
        %pnrm2 = norm(XtBX_I(:))^2; % Check
        fprintf(['Outer %2i: maxres %8.2e (%3i) obj %15.8e ninc (%i) ',...
            'Proj%i  %6.2fs\n'],outer,maxres,out.iter,out.obj(end),...
            out.ninc,opts.proj,toc(t0)); % Check
        if  opts.print > 1, fprintf('\n'); end % Check
    end % Check    
    
    iterCount.inner = iterCount.inner + out.iter;
    iterCount.outer = iterCount.outer + 1;
    
    if maxres < opts.tol2 && outer > 1, break; end % Check
    if outer >= maxouter; break; end % Check
    %opts.mu = max(1,ev(kw));
    opts.tol1 = opts.tol1/tol1sp; % Check


end

% output
X = X(:,1:k); % Check
D = D(1:k,1:k); % Check
out.ev = ev(1:k); % Check
out.Res = Res(1:k); % Check
maxres = max(Res(1:k)); % Check
out.maxres = maxres; % Check
out.iter = iterCount;
out.geigpen_times = geigpen_times; 
%%%%%%%% end of main body %%%%%%%%%


%%%%%%%%% embedded function %%%%%%%%%%%%%%%%%%%%%%%%
    function initialize % Check
        
        if ~isfield(opts,'tol1');       opts.tol1 = 1e-1;       end % Check
        if ~isfield(opts,'tol2');       opts.tol2 = 1e-4;       end % Check
        if ~isfield(opts,'maxit');      opts.maxit = 100;       end % Check
        if ~isfield(opts,'print');      opts.print = 0;         end % Check
        if ~isfield(opts,'proj');       opts.proj = 1;          end % Check
        if ~isfield(opts,'dols');       opts.dols = 2;          end % Check
        if ~isfield(opts,'arrp');       opts.arrp = 0;          end % Check
        if ~isfield(opts,'kp');         opts.kp = round(.2*k);  end % Check
        if ~isfield(opts,'maxouter');   opts.maxouter = 20;     end % Check
        if ~isfield(opts,'tol1sp');     opts.tol1sp = 5;        end % Check
        if ~isfield(opts,'istop');      opts.istop = false;     end % Check
        if ~isfield(opts,'freq');       opts.freq = 5;          end % Check
        if ~isfield(opts,'rho');        opts.rho = 0.2;         end % Check
        if ~isfield(opts,'alf_bool');   opts.alf_bool = 1;      end % Check
        if ~isfield(opts,'outer');      opts.outer = 0;         end % Check
        % add extra kp guard vectors
        kw = k + opts.kp; % Check
                
        % generate random X if necessary
        if ~exist('X','var') || isempty(X) % Check
            X = randn(n,kw); % Check
        end % Check
        % scale X columns to unit size
        d = 1./sqrt(sum(X.^2)); % Check
        X = bsxfun(@times,X,d); % Check
        
        % set parameter pmu
        tic;
        [X,~,ev,Res,AX,BX, gRR_times] = gRR(A,B,X); % Check
        geigpen_times.gRR = toc; 
     
        geigpen_times.AX = gRR_times.AX;
        geigpen_times.BX = gRR_times.BX;
        geigpen_times.qr = gRR_times.qr;
        geigpen_times.XtAX = gRR_times.XtAX;
        geigpen_times.XtBX = gRR_times.XtBX;
        geigpen_times.eig = gRR_times.eig;
        geigpen_times.XY = gRR_times.XY;
        geigpen_times.AXY = gRR_times.AXY;
        geigpen_times.BXY = gRR_times.BXY;
        geigpen_times.Res = gRR_times.Res;
        
        pmu = max(ev); % Check

        if ~isfield(opts,'mu') % Check
            pmu = abs(pmu); % Check
            opts.mu = pmu; % Check
        end % Check
        
    end %initialize % Check

%%%%%%%%% embedded function %%%%%%%%%%%%%%%%%%%%%%%%
    function [X,D,ev,Res,AX,BX,gRR_times] = gRR(A,B,X) % Check
% compute Ritz pairs
        gRR_times = [];

        tic; [X,~] = qr(X,0); gRR_times.qr = toc; % Check
        tic; AX = A*X; gRR_times.AX = toc; % Check
        tic; XtAX = X'*AX; gRR_times.XtAX = toc; XtAX = symm(XtAX); % Check
        tic; BX = B*X; gRR_times.BX = toc; 
        tic; XtBX = X'*BX; gRR_times.XtBX = toc; XtBX = symm(XtBX); % Check
        
        tic; [V,D] = eig(XtAX,XtBX); gRR_times.eig = toc; % Check
        
        ev = real(diag(D));
        [ev,io] = sort(ev,'ascend'); % Check
        
        tic; X  =  X*V(:,io); gRR_times.XY = toc; % Check
        tic; AX = AX*V(:,io); gRR_times.AXY = toc; % Check
        tic; BX = BX*V(:,io); gRR_times.BXY = toc; % Check
        
        D = D(io,io); % Check
        
        tic;
        Res = sqrt(sum((AX-BX*D).^2))./max(1,sqrt(sum(AX.^2))); % Check       
        Res = Res'; % Check
        gRR_times.Res = toc;
        
    end %RR % Check

end % Check
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    %% solver
    function [X, G, out] = gpen_solve(A,B,X,AX,BX,opts) % Check
    %
    %function [X, G, out] = gpen_solve(A,B,X,AX,BX,opts)
    %
    % Solving: min tr(X'AX) + 0.5*mu*||X'BX-I||^2_F
    %
    % Input:
    %        A -- (n x n) symmetric matrix
    %        B -- (n x n) symmetric matrix
    %        k -- (1 x 1) dimension of invariant subspace
    %        X -- (n x k) initial matrix; could be empty
    %       AX -- (n x k) A*X
    %       BX -- (n x k) B*X
    % Output:
    %        X -- (n x k) computed solution

    %% set parameters
    pmu = opts.mu; % Check
    tol = opts.tol1; % Check
    maxit = opts.maxit; % Check
    do_proj = opts.proj; % Check
    alf_bool = opts.alf_bool; % Check
    dols = opts.dols; % Check
    istop = opts.istop; % Check
    freq = opts.freq; % Check
    rho = opts.rho; % Check
    
    ninc = 0; % Check
    
    gpen_times = [];
    
    Ik = speye(size(X,2)); % Check


    
    tic; trXAX = dot(X(:),AX(:)); gpen_times.dot = toc; % Check
    tic; XtBX_I = X'*BX; gpen_times.XtBX_I = toc; 
    XtBX_I = XtBX_I - Ik;% Check
    
    pnrm2 = norm(XtBX_I(:))^2; % Check
    obj  = trXAX/2 + pmu/4 * pnrm2; % Check
    if opts.print > 2 % Check
        fprintf('objs:\t[%9.2e %9.2e %9.2e]\n',trXAX,pnrm2,obj); % Check
    end % Check
    out.alf  = zeros(1,maxit); % Check
    out.obj  = zeros(1,maxit); % Check

    % gradient or projected gradient 
        if do_proj % Check 
            tic;
            R = (BX'*BX)\(BX'*AX); % Check
            G = (AX - BX*R); % Check
            gpen_times.G = toc; % Check
        else
            tic; G = AX + BX*(pmu*XtBX_I); gpen_times.G = toc; % Check
        end % Check

    alf = 1e-4*norm(X(:))/norm(G(:)); % Check
    % Preallocate Memory
    
    gpen_times.alfSpike = 0;
    gpen_times.AX = 0;
    gpen_times.BX = 0;
    gpen_times.backtrack = 0;
    gpen_times.rGnrm = 0;
    
    %% ---------- loop -----------
    for iter = 1:maxit % Check
        
        objp = obj; Xd = X; Gd = G; % save obj, X and G % Check
        
        %% control alf not to spike too much
        tic;        
        if iter > 2 && alf_bool % Check
            alf = min(alf,25*mean(out.alf(iter-2:iter)));  % Check
        end % Check
        gpen_times.alfSpike = gpen_times.alfSpike + toc;
        
        %% update X with a backtracking line search
        tic;
        for backtrack = 0:dols % Check
            X = Xd - alf*G; % Check %Should it be + alf*G
            tic; 
            AX = A*X;
            gpen_times.AX = gpen_times.AX + toc; % Check
            
            tic; 
            BX = B*X; 
            gpen_times.BX = gpen_times.BX + toc; % Check
            
            tic; 
            trXAX = dot(X(:),AX(:)); % Check
            gpen_times.dot = gpen_times.dot + toc;
            
            tic;
            XtBX_I = X'*BX;
            gpen_times.XtBX_I = gpen_times.XtBX_I + toc;
            XtBX_I = XtBX_I - Ik; % Check
            
            pnrm2 = norm(XtBX_I(:))^2; % Check
            obj = trXAX/2 + pmu/4 * pnrm2; % Check
            ok = obj < objp; % Check
            if ok; break; end % Check
            alf = rho*alf; % Check
        end % Check
        
        gpen_times.backtrack = gpen_times.backtrack + toc;
        
        if obj >= objp, ninc = ninc + 1; end % Check
        
        %% update (projected) gradient
        %do_proj = mod(iter,2) == 0;
        if do_proj % Check
            tic;
            R = (BX'*BX)\(BX'*AX); % Check
            G = (AX - BX*R); % Check
            gpen_times.G = gpen_times.G + toc;
        else
            tic;
            G = AX + BX*(pmu*XtBX_I); % Check;
            gpen_times.G = gpen_times.G + toc;
        end % Check
        
        %% record info
        out.obj(iter) = obj; % Check
        out.alf(iter) = alf; % Check

        %% check stopping rules
        tfimf = mod(iter,freq); % Check
        if tfimf == 0 % Check
            tic;
            rGnrm = max(sqrt(sum(G.^2))./max(1,sqrt(sum(AX.^2)))); % Check
    %        rGnrm = norm(G,'fro')/norm(AX,'fro'); % Check
            gpen_times.rGnrm = gpen_times.rGnrm + toc;
            
            istop = rGnrm < tol; % Check
        end % Check
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        % print iter information
        if opts.print > 1 && (~tfimf||iter == maxit||istop) % Check
            fprintf('iter: %3i  rGnrm: %6.2e  pnrm2: %6.2e  obj: ',...
                '%9.6e\n',iter,rGnrm,pnrm2,obj); % Check
        end     % Check
        
        % exit the loop
        if istop; break; end     % Check
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%% BB step length
        Xd = X - Xd; % Check
        Gd = G - Gd; % Check
        
        tic;
        Gdnrm2 = dot(Gd(:),Gd(:)); % Check
        dXdG = dot(Xd(:),Gd(:)); % Check
        gpen_times.dot = gpen_times.dot + toc;
        
        alf = abs(dXdG)/Gdnrm2; % Check
        %if dXdG <= 0, warning('dXdG <= 0'); end
            
    end % ----- loop ------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    out.opts = opts; % Check
    out.alf = out.alf(1:iter); % Check
    out.obj = out.obj(1:iter); % Check
    out.iter = iter; % Check
    out.ninc = ninc; % Check
%     out.Xd = Xd;
    out.gpen_times = gpen_times;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end % Check