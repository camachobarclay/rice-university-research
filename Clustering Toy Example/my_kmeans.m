function [C, mu2, MU, superrval,timeval] = my_kmeans(X,k,maxiter, tol)
   
    % This takes note of the size of the data X, which we will assume is an
    % N x p matrix.
    tic;
    [N, p] = size(X);
    
    % The k-means algorithm depends on having at least k distinct points to
    % put into the k distinct bins. Hence, we perform some auxilary error
    % checking to make sure that our clustering makes sense for the data X.
    
    UX = unique(X, 'rows');
    [UN, Up] = size(UX);
    
    % This tests the dimensions N,p and k to make sure that they make
    % sense.
    
    
    if N<1||p<1
        error('Matrix dimensions don''t make sense.\n');
    elseif UN<1||p<1
        error('UX matrix dimensions don''t make sense.\n');
    elseif UN>N||Up<p
        error(['Function ''''unique(X,''rows'')'''' isn''t doing what',...
            'you think it''s doing.\n'])
    elseif ~(mod(k,1)==0);
        error('Variable ''k'' is not an integer.\n');
    elseif (UN==N)&&(Up==p)
        if norm(X)==NaN||norm(X)==Inf
            error('Matrix ''X'' is poorly conditioned.\n');
        elseif k>N
            error('There are more bins than there are data points.\n');
        elseif k==N
            C = [1:N].';
            mu2 = X;
            MU = X;
        elseif 0<k
            [C, mu2,MU,superrval] = employkmeans(X,k,N,p,maxiter,tol);
        else
            error('Something went wrong, note that UN==N, and Up==p.\n');
        end
    elseif (0<UN)&&(UN<N)&&(Up==p)
        fprintf('Number of unique data points is strictly less than N.\n');
        if (norm(UX)==NaN)||(norm(UX)==Inf)
            error('Matrix ''UX'' is poorly conditioned.\n');
        elseif k>UN
            error(['Variable k is less than or equal to N but strictly',... 
                'greater than the number of distinct data points.']);
        elseif k==UN
            C = [1:UN].';
            mu2 = UX;
            MU = UX;
            fprintf(['Number of bins is equal to the number of distinct'... 
                'data points. MU does not hold k initial guesses',... 
                'for the centers. Instead, MU==UX.\n']);
        elseif 0<k
            [C, mu2,MU,superrval] = employkmeans(X,k,N,p,maxiter, tol);
            MU = UX;
            fprintf(['MU no longer holds initial k centers. Instead,',...
                'MU==UX. the unique set of data points.\n']);
        else
            error('Something went wrong, note that UN<N.\n');
        end
    else
        error('Neither conditional evaluated. Something went wrong.\n');
    end
    timeval = toc;
    
    
end


