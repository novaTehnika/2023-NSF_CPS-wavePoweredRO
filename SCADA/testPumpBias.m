
f = linspace(-1,1,101);

for i = 1:numel(f)
    [f1(i),f2(i)] = fcn(f(i),0.5,0.5);
    F(i) = mean([f1(i),f2(i)]);
end

figure
plot(f,f1)
hold on
plot(f,f2)
plot(f,F)
grid on
legend('f1','f2','F')

function [f1,f2] = fcn(f,HIL_H1T_P1fBias,HIL_H1T_P2fBias)
    P1 = -HIL_H1T_P1fBias;
    P2 = HIL_H1T_P2fBias;
    
    f1 = (1/(1+P1))*(f-P1) * (f<=P1) ...
       + (2*P2/(P2-P1))*(f-P1) * ((P1 < f) & (f < P2)) ...
       + ((1-2*P2)/(1-P2)*f + 1 - (1-2*P2)/(1-P2)) * (f>=P2);
    
    f2 = ((1+2*P1)/(1+P1)*f + (1+2*P1)/(1+P1) - 1) * (f<=P1) ...
       + (2*P1/(P1-P2))*(f-P2) * ((P1 < f) & (f < P2)) ...
       + (1/(1-P2))*(f-P2) * (f>=P2);
end