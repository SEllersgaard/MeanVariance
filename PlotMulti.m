%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This programme reads empircal data and plots the associated mean-variance
% frontier. It also simulates data into the future, estimates sample
% parameters and plots the mean-variance frontiers again. The point is that
% only be using the true drift (mu) do we get something that resembles the
% original mean-variance frontier.
%
% Simon Ellersgaard Nielsen, 24.11.15.
%
% Parameter specification *************************************************

years = 20;                     % Years of data to look at, ending July 2015
dt    = 1/252;                  % Time step
datapoints = years*(1/dt);      

years2 = 5;                     % Years of data to forecast
dt2    = 1/252;                 % Time step in forecast
datapoints2 = years2*(1/dt2);

simulations = 50;               % No. of simulations to run

ymax = 0.2;                     % Max expected return on figure
topx = 0.2;                     % Max expected varianc on figure

% Extract Data ************************************************************

Index = csvread('DataIndex.CSV');
RateData  = csvread('DataRiskFree.CSV',0,1);
Rate = RateData(1:length(Index),4);

% Shorten data ************************************************************

IndexEs = Index(end-datapoints+1:end,2:end);
for i = 1:1:datapoints   % clear missing data
   for j = 1:1:5
      if IndexEs(i,j) == -99.99 || IndexEs(i,j) == -999
          warning('Data points incomplete: set to zero')
          IndexEs(i,j) = 0;
      end
   end
end
RateEs  = Rate(end-datapoints+1:end,1);
IndexExcess = (IndexEs - RateEs*ones(1,5))/100; % excess return

% Compute estimators ******************************************************

mu = mean(IndexExcess)/dt;      % Drift of each of the five indices

Sigma = zeros(5,5);             % Covariance of each of the five indices
for i = 1:1:5
   for j = 1:1:5 
    
       sum = 0;
       for k = 1:1:datapoints
          sum = sum + (IndexExcess(k,i)-mu(i)*dt)*(IndexExcess(k,j)-mu(j)*dt); 
       end
       
       Sigma(i,j) = (1/dt)*1/(datapoints-1)*sum;
       
   end
end


% Plot true mean-variance frontier ****************************************

% These are the parameters the go into the mean-variance equation
Sigmainv = inv(Sigma);
atrue = mu*Sigmainv*mu';
btrue = mu*Sigmainv*ones(5,1);
ctrue = ones(1,5)*Sigmainv*ones(5,1);
dtrue = atrue*ctrue-btrue^2;

x = 0:0.0001:topx;              % Expected return axis

[y,y2] = sqfun(x,atrue,btrue,ctrue,dtrue);

% Plot the true mean-variance frontier for subgraphs **********************

for i = 1:1:5
figure(i)
plot(x,y,'k','LineWidth',2)
hold on
plot(x,y2,'k','LineWidth',2)
end

% Simulate data points ****************************************************

sig = chol(Sigma,'lower');          %Cholesky decompose covariance

muplot = zeros(simulations,5);      %Encode all the simulated means
varplot = zeros(simulations,15);    %Encode all the simulated covariances

musimmean = zeros(1,5);             %The mean of the simulated mus
Sigmasimmean = zeros(5,5);          %The mean of the simulated covariances

x = 0:0.0001:topx;

ABCD = zeros(simulations,4);        %Store the a,b,c,d coefficients with 
                                    %simulated mu and Sigma


for p = 1:1:simulations

    IndexSim = zeros(datapoints2,5);

    for i = 1:1:datapoints2
       IndexSim(i,:) = mu*dt2 + sqrt(dt2)*randn(1,5)*(sig'); 
    end

    mu2 = mean(IndexSim)/dt2;

    Sigma2 = zeros(5,5);
    for i = 1:1:5
       for j = 1:1:5 

           sum = 0;
           for k = 1:1:datapoints2
              sum = sum + (IndexSim(k,i)-mu2(i)*dt2)*(IndexSim(k,j)-mu2(j)*dt2); 
           end

           Sigma2(i,j) = (1/dt2)*1/(datapoints2-1)*sum;

       end
    end
    
    % Calculate the mean of the estimators
    musimmean = musimmean + mu2;
    Sigmasimmean = Sigmasimmean + Sigma2;
    
    % Store the estimators for plotting purposes
    muplot(p,:) = mu2;
    varplot(p,:) = [Sigma2(1,1),Sigma2(1,2),Sigma2(1,3),Sigma2(1,4),Sigma2(1,5),...
        Sigma2(2,2),Sigma2(2,3),Sigma2(2,4),Sigma2(2,5),Sigma2(3,3),Sigma2(3,4),...
        Sigma2(3,5),Sigma2(4,4),Sigma2(4,5),Sigma2(5,5)];
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Simulated mu, simulated Sigma
    
    figure(1)  
    Sigmainv = inv(Sigma2);
    atrue = mu2*Sigmainv*mu2';
    btrue = mu2*Sigmainv*ones(5,1);
    ctrue = ones(1,5)*Sigmainv*ones(5,1);
    dtrue = atrue*ctrue-btrue^2;
    
    ABCD(p,:) = [atrue,btrue,ctrue,dtrue];

    [y,y2] = sqfun(x,atrue,btrue,ctrue,dtrue);

    plot(x,y,'-.k')
    hold on
    plot(x,y2,'-.k')
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     % Simulated mu, true Sigma
    
    figure(2)
    Sigmainv = inv(Sigma);
    atrue = mu2*Sigmainv*mu2';
    btrue = mu2*Sigmainv*ones(5,1);
    ctrue = ones(1,5)*Sigmainv*ones(5,1);
    dtrue = atrue*ctrue-btrue^2;

    [y,y2] = sqfun(x,atrue,btrue,ctrue,dtrue);

    plot(x,y,'-.k')
    hold on
    plot(x,y2,'-.k')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     % True mu, simulated Sigma
    figure(3)
    Sigmainv = inv(Sigma2);
    atrue = mu*Sigmainv*mu';
    btrue = mu*Sigmainv*ones(5,1);
    ctrue = ones(1,5)*Sigmainv*ones(5,1);
    dtrue = atrue*ctrue-btrue^2;

    [y,y2] = sqfun(x,atrue,btrue,ctrue,dtrue);

    plot(x,y,'-.k')
    hold on
    plot(x,y2,'-.k')

end

% Specify grapg properties*************************************************

figure(1)
xlim([0,0.2])
ylim([0,ymax])
xlabel('Variance','FontSize',14)
ylabel('Expected Excess Return','FontSize',14)
hti1 = title('Simulated \mu, Simulated \Sigma  ');

figure(2)
xlim([0,0.2])
ylim([0,ymax])
xlabel('Variance','FontSize',14)
ylabel('Expected Excess Return','FontSize',14)
hti2 = title('Simulated \mu, True \Sigma  ');

figure(3)
xlim([0,0.2])
ylim([0,ymax])
xlabel('Variance','FontSize',14)
ylabel('Expected Excess Return','FontSize',14)
hti3 = title('True \mu, Simulated \Sigma  ');

figure(4)
xlim([0,0.2])
ylim([0,ymax])
xlabel('Variance','FontSize',14)
ylabel('Expected Excess Return','FontSize',14)
hti4 = title('Mean simulated \mu, Mean simulated \Sigma  ');

figure(5)
xlim([0,0.2])
ylim([0,ymax])
xlabel('Variance','FontSize',14)
ylabel('Expected Excess Return','FontSize',14)
hti5 = title('Average Mean-Variance Frontier  ');

set(hti1,'FontSize',14)
set(hti2,'FontSize',14)
set(hti3,'FontSize',14)
set(hti4,'FontSize',14)
set(hti5,'FontSize',14)

% Plot mean variance frontier with the average estimator*******************

% Simulated mu, true Sigma
figure(4)
musimmean = musimmean/simulations;
Sigmasimmean = Sigmasimmean/simulations;

Sigmainv = inv(Sigmasimmean);
atrue = musimmean*Sigmainv*musimmean';
btrue = musimmean*Sigmainv*ones(5,1);
ctrue = ones(1,5)*Sigmainv*ones(5,1);
dtrue = atrue*ctrue-btrue^2;

[y,y2] = sqfun(x,atrue,btrue,ctrue,dtrue);

plot(x,y,'-.k')
hold on
plot(x,y2,'-.k')

% Plot the "mean" mean-variance frontier (horizontal average) for estimated
% mu and Sigma. ***********************************************************

figure(5)
varm = zeros(1,length(x));
for i = 1:1:length(x) % We pretend this is the y axis now
   
    sumINV = 0;
    
    for j = 1:1:simulations  % sum over all curves 
       
        % For each curve we take the inverse function to get a variance. We
        % then compute the average variance
        sumINV = sumINV + (ABCD(j,1)-2*ABCD(j,2)*x(i)+ABCD(j,3)*x(i).^2)/ABCD(j,4);
        
    end
    
    varm(1,i) = sumINV/simulations;
    
end
plot(varm,x,'-.k')

% Plot the running estimators mu and Sigma against their true value *******

muplottrue = zeros(simulations,5);
varplottrue = zeros(simulations,15);

for p = 1:1:simulations
muplottrue(p,:) = mu;
varplottrue(p,:) = [Sigma(1,1),Sigma(1,2),Sigma(1,3),Sigma(1,4),Sigma(1,5),...
        Sigma(2,2),Sigma(2,3),Sigma(2,4),Sigma(2,5),Sigma(3,3),Sigma(3,4),...
        Sigma(3,5),Sigma(4,4),Sigma(4,5),Sigma(5,5)];
end

number = 1:1:simulations;

figure(6) % plot mu
for i = 1:1:5
plot(number,muplot(:,i)','k')
hold on
plot(number,muplottrue(:,i)',':k')
hold on
end
xlabel('Simulation #','FontSize',14)
ylabel('Estimated \mu','FontSize',14)
ax = gca;
ax.XTick = [1,5:5:simulations];
%title('Simulated \mu, True \Sigma  ')
xlim([1,simulations])
text(simulations+1,mu(4)-0.02,'_1')
text(simulations+1,mu(4)-0.01,'_2')
text(simulations+1,mu(4)-0.03,'_3')
text(simulations+1,mu(4),'_4')
text(simulations+1,mu(4)-0.04,'_5')

figure(7) % Plot Sigma
for i = 1:1:15
plot(number,varplot(:,i)','k')
hold on
plot(number,varplottrue(:,i)',':k')
hold on
end
xlabel('Simulation #','FontSize',14)
ylabel('Estimated \Sigma','FontSize',14)
ax = gca;
ax.XTick = [1,5:5:simulations];
xlim([1,simulations])

text(simulations+1,Sigma(1,1),'_{11}')
text(simulations+1,Sigma(1,2),'_{12}')
text(simulations+1,Sigma(1,3),'^{13}')
text(simulations+1,Sigma(1,4),'_{14}')
text(simulations+1,Sigma(1,5),'_{15}')
text(simulations+1,Sigma(2,2),'^{22}')
text(simulations+1,Sigma(2,3),'_{23}')
text(simulations+1,Sigma(2,4),'_{24}')
text(simulations+1,Sigma(2,5),'_{25}')
text(simulations+1,Sigma(3,3),'_{33}')
text(simulations+1,Sigma(3,4),'_{34}')
text(simulations+1,Sigma(3,5),'_{35}')
text(simulations+1,Sigma(4,4),'^{44}')
text(simulations+1,Sigma(4,5),'_{45}')
text(simulations+1,Sigma(5,5),'_{55}')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



