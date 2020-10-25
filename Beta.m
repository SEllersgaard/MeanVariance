years = 20;
dt    = 1/252;
datapoints = years*(1/dt);

% Extract Data ************************************************************
Index49 = csvread('DataIndex49.CSV');
RateData  = csvread('DataRiskFree.CSV',0,1);

RmRf = RateData(1:length(Index49),1);
Rate = RateData(1:length(Index49),4);

% Shorten data ************************************************************
IndexEs = Index49(end-datapoints+1:end,2:end);
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
RmRfExcess = RmRf(end-datapoints+1:end,1)/100;
IndexExcess = (IndexEs - RateEs*ones(1,49))/100; % excess return

% Calculate beta vector ***************************************************

covRiRm = zeros(1,49);

meanRmRf = mean(RmRfExcess)/dt;
meanInEx = mean(IndexExcess)/dt;

%varRmRf = var(RmRfExcess)/d;

varsum = 0;
for j = 1:1:datapoints
    varsum = varsum + (RmRfExcess(j)-meanRmRf*dt)*(RmRfExcess(j)-meanRmRf*dt);
end
varRmRf = varsum/((datapoints-1)*dt);



for t = 1:1:49
    
    covsum = 0;
    for j = 1:1:datapoints
       covsum = covsum + (IndexExcess(j,t)-meanInEx(1,t)*dt)*(RmRfExcess(j)-meanRmRf*dt); 
    end
    covRiRm(1,t) = covsum/((datapoints-1)*dt); 
    
end

beta = covRiRm/varRmRf;

%meanInEx(11) = [];
%beta(11) = [];

OLSbetatop = 0;
OLSbetabot = 0;
betamean = mean(beta);
meanmeanInEx = mean(meanInEx);
for i = 1:1:49
   OLSbetatop = OLSbetatop + (beta(i)-betamean)*(meanInEx(i)-meanmeanInEx);
   OLSbetabot = OLSbetabot + (beta(i)-betamean)*(beta(i)-betamean);
end
OLSbeta = OLSbetatop/OLSbetabot;
OLSalpha = meanmeanInEx - OLSbeta*betamean;

x = 0:0.1:1.6;
y = OLSalpha + OLSbeta*x;


plot(x,y,'k')
xlabel('Beta','FontSize',14)
ylabel('Expected Excess Return','FontSize',14)

mdl = fitlm(x,y)

%lsline
hold on

y2 = meanRmRf*x;
plot(x,y2,'--k')
hold on

h2 = legend('Empirical SML ','Theoretical SML ');
set(h2,'FontSize',14)


scatter(beta,meanInEx,'Xk')







