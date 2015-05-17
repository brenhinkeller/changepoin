%% Make synthetic dataset

data=[randn(30,1); randn(50,1)+6; randn(100,1)+2; randn(47,1)+4];
figure; plot(data)
hold on; plot([1 30 31 80 81 180 181 227],[0 0 6 6 2 2 4 4],'m')


%% Brute force version
npmin=0; % Minimum number of changepoints
npmax=3; % Maximum number of changepoints


% dmin=2; % Minimum number of points needed between any two changepoints (in any partition)

possiblepoints=size(data,1)-1; % Number of possible changepoint locations


nsims=20000;
r2C=NaN(nsims,1);
r2C(1)=Inf; % Current sum squared residual

for i=2:nsims;
    % Propose a random number of changepoints
    np=round(rand(1,1)*(npmax-npmin)+npmin);

    % Propose a random set of boundary points
    bndPntsP=[0; sort(ceil(rand(np,1)*possiblepoints)); size(data,1)];
    
    % means=NaN(npmax+2,1);
    meansP=NaN(length(bndPntsP)-1,1);
    
    
    % Add a new changepoint
    
    % Delete a changepoint
    
    % Move a changepoint
    
    
    % Calculate probability of accepting proposed model
    
    
    % Can evolve the acceptance criteria....
    
    % Fit between changepoints and recalculate residuals
    r2P=0;
    for k=1:length(bndPntsP)-1;
        meansP(k)=mean(data(bndPntsP(k)+1:bndPntsP(k+1)));
        r2P=r2P+sum((data(bndPntsP(k)+1:bndPntsP(k+1))-meansP(k)).^2);
    end
    
    if r2P<r2C(i-1)
        r2C(i)=r2P;
        bndPntsC=bndPntsP;
        meansC=meansP;
    else
        r2C(i)=r2C(i-1);
    end
    
end

%% Plot results
% figure; plot(find(data),data)
model=sortrows([bndPntsC(1:end-1)+1 meansC; bndPntsC(2:end) meansC]);
hold on, plot(model(:,1), model(:,2),'r');

figure; plot(r2C)


