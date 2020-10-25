function [Y,Y2] = sqfun(x,atrue,btrue,ctrue,dtrue)

y = zeros(length(x),1);         % Variance axis (+ve solution)
y2 = zeros(length(x),1);        % Variance acis (-ve solution)

% Plot the two solutions
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
% Fix the gap in the parabola
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

Y = y;
Y2 = y2;


end

