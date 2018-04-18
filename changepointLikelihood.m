%% Make synthetic dataset
tic;
clear
data=[randn(300,1); randn(500,1)+6; randn(1000,1)+2; randn(470,1)+4];%; randn(100,1)+2; randn(100,1)-5];
data=repmat(data,2 ,1);
% save data data;
% load data
data=data-mean(data);
data=data/std(data);

figure; plot(data)
% hold on; plot([1 30 31 80 81 180 181 227],[0 0 6 6 2 2 4 4],'m')

%% Markov chain version using data uncertainties
possiblepoints=size(data,1)-1; % Number of possible changepoint locations

npmin=0; % Minimum number of changepoints
npmax=possiblepoints; % Maximum number of changepoints


dmin=10; % Minimum number of points needed between any two changepoints (in any partition)


nsims=10^6;
llC=NaN(nsims,1);
llC(1)=-Inf; % Current sum squared residual


% Propose a random number of changepoints
% np=round(rand(1,1)*(npmax-npmin)+npmin);
np=NaN(nsims,1);
np(1)=1;


r2C=NaN(nsims,1);
r2C(1)=NaN; % Current sum squared residual

% Propose an initial random set of boundary points
bndPntsC=[0; sort(ceil(rand(np(1),1)*possiblepoints)); size(data,1)];

% means=NaN(npmax+2,1);
meansP=NaN(length(bndPntsC)-1,1);

varC=NaN(size(data));
for k=1:length(bndPntsC)-1;
    varC(bndPntsC(k)+1:bndPntsC(k+1))=var(data(bndPntsC(k)+1:bndPntsC(k+1)));
end



for i=2:nsims;
    choice=ceil(rand(1,1)*3);
    bndPntsP=bndPntsC;
    if choice==1&&np(i-1)<npmax;
            % Add a new changepoint
            bndPntsP=unique(sort([bndPntsP; ceil(rand(1,1)*possiblepoints)]));
            meansP=NaN(length(bndPntsP)-1,1); 
    elseif choice==2&&np(i-1)>npmin
            % Delete a changepoint
            unlucky=ceil(rand(1,1)*(length(bndPntsC)-2)+1);
            bndPntsP(unlucky)=[];
            meansP=NaN(length(bndPntsP)-1,1);

    elseif choice==3&&np(i-1)>0
            % Move a changepoint
            lucky=ceil(rand(1,1)*(length(bndPntsP)-2)+1);
            bndPntsP(lucky)=ceil(rand(1,1)*possiblepoints);
            bndPntsP=unique(sort(bndPntsP));
    end
    
    % Remove adjacent changepoints
    test=bndPntsP(2:end-1)-bndPntsP(1:end-2)<=dmin;
    if sum(test)>1;
        bndPntsP([false; test; false])=[];
        meansP=NaN(length(bndPntsP)-1,1);
    end
    
    % Can evolve the acceptance criteria....
    
    % calculate log likelyhood of new model fitting data given last variance
%     llP=-2*log10(possiblepoints^length(bndPntsP))-log10(factorial(length(bndPntsP)));
    llP=-log10(possiblepoints^length(bndPntsP))-log10(nchoosek(possiblepoints,length(bndPntsP)));
    r2P=0;
    for k=1:length(bndPntsP)-1;
        meansP(k)=mean(data(bndPntsP(k)+1:bndPntsP(k+1)));
        llP = llP + sum( log10( 2*normcdf( -abs(data(bndPntsP(k)+1:bndPntsP(k+1)) - meansP(k)), 0, varC(bndPntsP(k)+1:bndPntsP(k+1)) )));
        r2P=r2P+sum((data(bndPntsP(k)+1:bndPntsP(k+1))-meansP(k)).^2);
    end
    
    
    % Calculate random probability of accepting proposed model
    if rand(1,1)<10^(llP-llC(i-1));
        % If accepted
        bndPntsC=bndPntsP;
        meansC=meansP;
        np(i)=length(bndPntsP)-2;
        
        % Update variance
        for k=1:length(bndPntsC)-1;
            varC(bndPntsC(k)+1:bndPntsC(k+1))=var(data(bndPntsC(k)+1:bndPntsC(k+1)));
        end
        
        % Update log likelyhood given new variance
%         llP=-2*log10(possiblepoints^length(bndPntsP))-log10(factorial(length(bndPntsP)));
        llP=-log10(possiblepoints^length(bndPntsP))-log10(nchoosek(possiblepoints,length(bndPntsP)));
        for k=1:length(bndPntsP)-1;
            meansP(k)=mean(data(bndPntsP(k)+1:bndPntsP(k+1)));
            llP = llP + sum( log10( 2*normcdf( -abs(data(bndPntsP(k)+1:bndPntsP(k+1)) - meansP(k)), 0, varC(bndPntsP(k)+1:bndPntsP(k+1)) )));
        end
        llC(i)=llP;
        
%         % Update figure
%         hold off; plot(data)
%         hold on; plot([1 30 31 80 81 180 181 227],[0 0 6 6 2 2 4 4],'m')
%         model=sortrows([bndPntsC(1:end-1)+1 meansC; bndPntsC(2:end) meansC]);
%         hold on, plot(model(:,1), model(:,2),'r');
        

        r2C(i)=r2P;
    
    else
        r2C(i)=r2C(i-1);
        llC(i)=llC(i-1);
        np(i)=np(i-1);

    end 
    
end

%% Plot results
% figure; plot(find(data),data)
model=sortrows([bndPntsC(1:end-1)+1 meansC; bndPntsC(2:end) meansC]);
hold on, plot(model(:,1), model(:,2),'r'); title('changepoint w/ liklihood')

figure; plot(llC); ylabel('log likelihood'); title('changepoint w/ liklihood')
figure; plot(r2C); ylabel('R squared'); title('changepoint w/ liklihood')
figure; plot(np); ylabel('number of changepoints'); title('changepoint w/ liklihood')

toc