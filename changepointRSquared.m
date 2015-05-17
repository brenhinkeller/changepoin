%% Make synthetic dataset
tic;
clear
% data=[randn(30,1); randn(50,1)+6; randn(100,1)+2; randn(47,1)+4];%; randn(100,1)+2; randn(100,1)-5];
% save data data;
load data
figure; plot(data)
hold on; plot([1 30 31 80 81 180 181 227],[0 0 6 6 2 2 4 4],'m')

%% Markov chain version
npmin=0; % Minimum number of changepoints
npmax=3; % Maximum number of changepoints


% dmin=2; % Minimum number of points needed between any two changepoints (in any partition)

possiblepoints=size(data,1)-1; % Number of possible changepoint locations


nsims=10000;
r2C=NaN(nsims,1);
r2C(1)=Inf; % Current sum squared residual


% Propose a random number of changepoints
% np=round(rand(1,1)*(npmax-npmin)+npmin);
np=1;

% Propose an initial random set of boundary points
bndPntsC=[0; sort(ceil(rand(np,1)*possiblepoints)); size(data,1)];

% means=NaN(npmax+2,1);
meansP=NaN(length(bndPntsC)-1,1);

for i=2:nsims;

    choice=ceil(rand(1,1)*3);
    bndPntsP=bndPntsC;
    if choice==1&&np<npmax;
            % Add a new changepoint
            bndPntsP=unique(sort([bndPntsP; ceil(rand(1,1)*possiblepoints)]));
            meansP=NaN(length(bndPntsP)-1,1); 
    elseif choice==2&&np>npmin
            % Delete a changepoint
            unlucky=ceil(rand(1,1)*(length(bndPntsC)-2)+1);
            bndPntsP(unlucky)=[];
            meansP=NaN(length(bndPntsP)-1,1);

    elseif choice==3
            % Move a changepoint
            lucky=ceil(rand(1,1)*(length(bndPntsP)-2)+1);
            bndPntsP(lucky)=ceil(rand(1,1)*possiblepoints);
            bndPntsP=unique(sort(bndPntsP));
    end
    
    
    % Can evolve the acceptance criteria....
    
    % Fit between changepoints and recalculate residuals
    r2P=0;
    for k=1:length(bndPntsP)-1;
        meansP(k)=mean(data(bndPntsP(k)+1:bndPntsP(k+1)));
        r2P=r2P+sum((data(bndPntsP(k)+1:bndPntsP(k+1))-meansP(k)).^2);
    end
    
    
    % Calculate probability of accepting proposed model
    if r2P<r2C(i-1)
        r2C(i)=r2P;
        bndPntsC=bndPntsP;
        meansC=meansP;
        np=length(bndPntsP)-2;
    else
        r2C(i)=r2C(i-1);
    end
    
end

%% Plot results
% figure; plot(find(data),data)
model=sortrows([bndPntsC(1:end-1)+1 meansC; bndPntsC(2:end) meansC]);
hold on, plot(model(:,1), model(:,2),'r');

figure; plot(r2C); ylabel('R squared')

toc
