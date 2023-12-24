function [C2, mu2, MU, superrval]=employkmeans(X,k,N,p,maxiter,tol)

    % Here we randomly choose k data points to be the initial cluster
    % centers. We require that they be distinct.     
    Z = X(randperm(N).',:);
    mu1 = Z(randi(N,1,k),:);  
    
    % This makes sure that the k chosen centers are all distinct.    
    Y = zeros(k,p);    
    while sum(size(unique(mu1,'rows'))==[k,p])-2
            mu1 = mu1-Y;
            Y = randn(k,p);
            mu1 = mu1+Y;
            fprintf('Still in first loop.\n');   
    end
% mu1 = [24, 84; 14,56; 17,35; 42,79; 34,56; 40,40; 32,17; 67,86; 61,56;
%     63,40; 51,18; 83,73; 86,55; 80,33; 85,15];

    % MU is the initial set of starting centers.    
    MU  = mu1;
    % C1 and C2 will store the old and the new clusters of the data.    
    C1 = zeros(N,1);
    C2 = C1;
    % These values will be used to decide which bin the ith data point
    % belongs to.    
    distvals = zeros(k,1);
    minval = 0;
    ind = 0;      
    % A counter that counts the number of outer iterations in the following
    % while loop.
    counter = 0;
    superrval = zeros(1,N);
	for i=1:N
      distvals = sqrt(sum(abs(mu1 - repmat(X(i,:),k,1)).^2,2));
%         distvals = sum(abs(mu1 - repmat(X(i,:),k,1)),2);
        [minval, ind] = min(distvals);
        C1(i) = ind;
    end
    % Now we recompute the centers of the k different clusters.
    for j=1:k
      mu2(j,:) = mean(X(C1(:)==j,:));
%         mu2(j,:) = median(X(C1(:)==j,:));
    end    
    Y = zeros(k,p);
    while sum(size(unique(mu2,'rows'))==[k,p])-2
            mu2 = mu2-Y;
            Y = randn(k,p);
            mu2 = mu2+Y;
            fprintf('Still in first loop.\n');            
    end
    for i=1:N
      distvals = sqrt(sum(abs(mu2 - repmat(X(i,:),k,1)).^2,2));
%         distvals = sum(abs(mu2 - repmat(X(i,:),k,1)),2);
        % Then we find the minimim distance and the corresponding index
        % and store it. This places each x into a cluster.
        [minval, ind] = min(distvals);
        C2(i) = ind;
    end
    while 1 
        counter = counter+1;
        % The following is our breaking parameter.
        superrval(counter) = max(max(abs(mu1(C1)-mu2(C2))));
        if (superrval(counter)<tol)||(counter>maxiter)
            fprintf('Breaking now...\n');
            break;
        end
        % Now we are looking at each vector X(i,:) for i=1,...,N
        mu1 = mu2;
        for j=1:k            
          mu2(j,:) = mean(X(C2(:)==j,:));
%             mu2(j,:) = median(X(C2(:)==j,:));
        end
        while sum(size(unique(mu2,'rows'))==[k,p])-2
                mu2 = mu2-Y;
                Y = randn(k,p);
                mu2 = mu2+Y;
                fprintf('Still in third loop.\n');            
        end
        C1 = C2;
        for i=1:N
          distvals = sqrt(sum(abs(mu2 - repmat(X(i,:),k,1)).^2,2));            
%             distvals = sum(abs(mu2 - repmat(X(i,:),k,1)),2);
            % Then we find the minimim distance and the corresponding index
            % and store it. This places each x into a cluster.
            [minval, ind] = min(distvals);
            C2(i) = ind;
        end
        % Now we recompute the centers of the k different clusters.] 
    end
    % This makes sure that we have every bin with at least one data point.
    superrval = superrval(1:counter);
    if length(unique(C2(:)))<k
        fprintf('Certain bins missing data points, will employ'+...
        'the k-means algorithm recursively.\n');
        [C2, mu2, MU,superrval] = employkmeans(X,k,N,p,maxiter,tol);
    end  
end