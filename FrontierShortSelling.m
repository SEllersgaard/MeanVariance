

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

% Sigmainv = inv(Sigma);
% atrue = mu*Sigmainv*mu';
% btrue = mu*Sigmainv*ones(5,1);
% ctrue = ones(1,5)*Sigmainv*ones(5,1);
% dtrue = atrue*ctrue-btrue^2;
% 
% 
% x = 0:0.0001:0.2;
% y = zeros(length(x),1);
% y2 = zeros(length(x),1);
% for i=1:1:length(x)
%     if (1/ctrue)*(dtrue*x(i)+btrue^2/ctrue-atrue) >= 0
%         y(i) = (btrue/ctrue) + sqrt((1/ctrue)*(dtrue*x(i)+btrue^2/ctrue-atrue));
%     else
%         y(i) = -100;
%     end
%     if (1/ctrue)*(dtrue*x(i)+btrue^2/ctrue-atrue) >= 0
%         y2(i) = (btrue/ctrue) - sqrt((1/ctrue)*(dtrue*x(i)+btrue^2/ctrue-atrue));
%     else
%         y2(i) = -100;
%     end
% end
% for i = 1:1:length(x)-1
%     if y(i) == -100 && y(i+1) > -100
%         jump = i;
%         y(jump) = (btrue/ctrue);
%     end
%     
%     if y2(i) == -100 && y2(i+1) > -100
%         jump2 = i;
%         y2(jump2) = (btrue/ctrue);
%     end
%     
% end
% 
% y(y==-100) = NaN;
% y2(y2==-100) = NaN;
% plot(x,y,'k','LineWidth',2)
% hold on
% plot(x,y2,'k','LineWidth',2)
% xlim([0,0.2])
% ylim([0,0.2])
% xlabel('Variance','FontSize',14)
% ylabel('Expected Excess Return','FontSize',14)
% title('True \mu, True \Sigma  ')
% 
% hold on
% 
% scatter(diag(Sigma)',mu,'Xk','LineWidth',2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dx = (max(mu)-min(mu))/100;

muP = min(mu):dx:max(mu);
%muP = 0.0850:0.001:0.105;
varP = zeros(1,length(muP));
varP2 = zeros(1,length(muP));

H = Sigma;
f = [];
A = [];
b = [];
Aeq = [mu;ones(1,5)];
lb = zeros(5,1);
ub = [];




for i = 1:1:length(muP)
    
    x0 = zeros(5,1);
    beq = [muP(i);1];
    
    %fun = @(x) 0.5*[x(1),x(2),x(3),x(4),x(5)]*Sigma*[x(1);x(2);x(3);x(4);x(5)];
    %pi = fmincon(fun,x0,A,b,Aeq,beq,lb,ub);
    
   
    pi = quadprog(H,f,A,b,Aeq,beq,lb,ub,ones(5,1));
    
    pi2 = quadprog(H,f,A,b,Aeq,beq,[],ub,ones(5,1));
    
    
   varP(1,i) = pi'*Sigma*pi;      
   varP2(1,i) = pi2'*Sigma*pi2;
   
end

plot(varP,muP,'-.k','LineWidth',2);
hold on
plot(varP2,muP,'k','LineWidth',2);
%xlim([0,0.2])
hold on
scatter(diag(Sigma)',mu,'Xr','LineWidth',2)
xlabel('Variance','FontSize',14)
ylabel('Expected Excess Return','FontSize',14)
hti = title('True \mu, True \Sigma  ');
set(hti,'FontSize',14)
h2 = legend('No Short Selling','With Short Selling');
set(h2,'FontSize',14)
%xlim([0,0.2])
%ylim([0,0.14])


[val,index] = max(muP./sqrt(varP));
[val2,index2] = max(muP./sqrt(varP2));

scatter(varP(index),muP(index),'MarkerEdgeColor','k','MarkerFaceColor','w','LineWidth',2)
scatter(varP2(index2),muP(index2),'MarkerEdgeColor','k','MarkerFaceColor','w','LineWidth',2)

Data = [varP(index),muP(index),muP(index)/sqrt(varP(index));...
    varP2(index2),muP(index2),muP(index2)/sqrt(varP2(index2))]






