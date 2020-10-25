

years = 20;
dt    = 1/252;
datapoints = years*(1/dt);

years2 = 5;
dt2    = 1/252;
datapoints2 = years2*(1/dt2);

ymax = 0.2;

% Extract Data ************************************************************
Index = csvread('DataIndex.CSV');
RateData  = csvread('DataRiskFree.CSV',0,1);
Rate = RateData(1:length(Index),4);


% Shorten data ************************************************************
IndexEs = Index(end-datapoints+1:end,2:end);
for i = 1:1:datapoints   % clear missing data
   for j = 1:1:5
      if IndexEs(i,j) == -99.99 || IndexEs(i,j) == -999
          i
          j
          IndexEs(i,j) = 0;
      end
   end
end
RateEs  = Rate(end-datapoints+1:end,1);
IndexExcess = (IndexEs - RateEs*ones(1,5))/100; % excess return

% Compute estimators ******************************************************
mu = mean(IndexExcess)/dt;

Sigma = zeros(5,5);
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

Sigmainv = inv(Sigma);
atrue = mu*Sigmainv*mu';
btrue = mu*Sigmainv*ones(5,1);
ctrue = ones(1,5)*Sigmainv*ones(5,1);
dtrue = atrue*ctrue-btrue^2;


x = 0:0.0001:0.2;
y = zeros(length(x),1);
y2 = zeros(length(x),1);
for i=1:1:length(x)
    if (1/ctrue)*(dtrue*x(i)+btrue^2/ctrue-atrue) > 0
        y(i) = (btrue/ctrue) + sqrt((1/ctrue)*(dtrue*x(i)+btrue^2/ctrue-atrue));
    else
        y(i) = -100;
    end
    if (1/ctrue)*(dtrue*x(i)+btrue^2/ctrue-atrue) > 0
        y2(i) = (btrue/ctrue) - sqrt((1/ctrue)*(dtrue*x(i)+btrue^2/ctrue-atrue));
    else
        y2(i) = -100;
    end
end
for i = 1:1:length(x)-1
    if y(i) == -100 && y(i+1) > -100
        jump = i;
        y(jump) = (btrue/ctrue);
    end
    
    if y2(i) == -100 && y2(i+1) > -100
        jump2 = i;
        y2(jump2) = (btrue/ctrue);
    end
    
end
y(y==-100) = NaN;
y2(y2==-100) = NaN;


%[y,y2] = sqfun(x,atrue,btrue,ctrue,dtrue);

figure(1)
plot(x,y,'k','LineWidth',2)
hold on
plot(x,y2,'k','LineWidth',2)


% Simulate data points ****************************************************

sig = chol(Sigma,'lower');


muplot = zeros(simulations,5);
varplot = zeros(simulations,15);


musimmean = zeros(1,5);
Sigmasimmean = zeros(5,5);


simnum = [5,50,1000];

for j = 1:1:length(simnum)
    
    yello = j;
    
    simulations = simnum(j);


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
    
    musimmean = musimmean + mu2;
    Sigmasimmean = Sigmasimmean + Sigma2;
      

end


figure(1)
xlim([0,0.2])
ylim([0,ymax])
xlabel('Variance','FontSize',14)
ylabel('Expected Excess Return','FontSize',14)
title('Mean simulated \mu, Mean simulated \Sigma  ')

% Simulated mu, true Sigma
figure(1)
musimmean = musimmean/simulations;
Sigmasimmean = Sigmasimmean/simulations;

Sigmainv = inv(Sigmasimmean);
atrue = musimmean*Sigmainv*musimmean';
btrue = musimmean*Sigmainv*ones(5,1);
ctrue = ones(1,5)*Sigmainv*ones(5,1);
dtrue = atrue*ctrue-btrue^2;

x = 0:0.0001:0.2;
y = zeros(length(x),1);
y2 = zeros(length(x),1);
for i=1:1:length(x)
    if (1/ctrue)*(dtrue*x(i)+btrue^2/ctrue-atrue) > 0
        y(i) = (btrue/ctrue) + sqrt((1/ctrue)*(dtrue*x(i)+btrue^2/ctrue-atrue));
    else
        y(i) = -100;
    end
    if (1/ctrue)*(dtrue*x(i)+btrue^2/ctrue-atrue) > 0
        y2(i) = (btrue/ctrue) - sqrt((1/ctrue)*(dtrue*x(i)+btrue^2/ctrue-atrue));
    else
        y2(i) = -100;
    end
end
for i = 1:1:length(x)-1
    if y(i) == -100 && y(i+1) > -100
        jump = i;
        y(jump) = (btrue/ctrue);
    end

    if y2(i) == -100 && y2(i+1) > -100
        jump2 = i;
        y2(jump2) = (btrue/ctrue);
    end

end
y(y==-100) = NaN;
y2(y2==-100) = NaN;

if yello == 1
plot(x,y,':k')
hold on
plot(x,y2,':k')
hold on
%legend('a')
elseif yello == 2
plot(x,y,'-.k')
hold on
plot(x,y2,'-.k')
hold on
%legend('b')
elseif yello == 3
plot(x,y,'--k')
hold on
plot(x,y2,'--k')
hold on
else
plot(x,y,':k')
hold on
plot(x,y2,':k')
hold on 
end

end

%h2 = legend('5 simulations','50 simulations','500 simulations','5000 simulations');
%set(h2,'FontSize',14)






