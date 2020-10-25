

years = 20; %10;
dt    = 1/252;
datapoints = years*(1/dt);

% Extract Data ************************************************************
Index = csvread('DataIndex.CSV');
RateData  = csvread('DataRiskFree.CSV',0,1);
Rate = RateData(1:length(Index),4);

% Shorten data ************************************************************
IndexEs = Index(end-datapoints+1:end,2:end);
RateEs  = Rate(end-datapoints+1:end,1);    
IndexExcess = (IndexEs- RateEs*ones(1,5))/100; % excess return

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
    if (1/ctrue)*(dtrue*x(i)+btrue^2/ctrue-atrue) >= 0
        y(i) = (btrue/ctrue) + sqrt((1/ctrue)*(dtrue*x(i)+btrue^2/ctrue-atrue));
    else
        y(i) = -100;
    end
    if (1/ctrue)*(dtrue*x(i)+btrue^2/ctrue-atrue) >= 0
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
plot(x,y,'k','LineWidth',2)
hold on
plot(x,y2,'k','LineWidth',2)
xlim([0,0.2])
ylim([0,0.2])
xlabel('Variance','FontSize',14)
ylabel('Expected Excess Return','FontSize',14)
hti = title('True \mu, True \Sigma  ');
set(hti,'FontSize',14)

hold on

scatter(diag(Sigma)',mu,'Xr','LineWidth',2) %'MarkerEdgeColor',[0.7 0.7 0.7]

% Plot Security Market Line ***********************************************

%rmean = mean(RateEs)/100;

%piCML = Sigmainv*mu'/(ones(1,5)*Sigmainv*mu');

%muCML = (mu+rmean*ones(1,5))*piCML;

%varCML = piCML'*Sigma*piCML;

[val,index] = max(y./sqrt(x')); %maximsie Sharpe ratio

muCML = y(index);
varCML = x(index);

scatter(varCML,muCML,'MarkerEdgeColor','k','MarkerFaceColor','w','LineWidth',2)

y3 = (muCML/sqrt(varCML))*sqrt(x);
% for i = 1:1:length(x)
%    if y3(i) > muCML
%       y3(i) = NaN; 
%    end
% end
plot(x,y3,'--k')
str1 = 'market portfolio';
text(varCML+0.005,muCML-0.005,str1)

str2 = 'capital market line';
text(varCML-0.015,muCML-0.05,str2)



