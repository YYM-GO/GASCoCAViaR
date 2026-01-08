function out = trafficlightonly2(p2n)
% test whether the score of A (VaR and CoVaR) is less than B(VaR and CoVaR)
%
% 1 - grey: A is significant superior for VaR
% 2 - green: insignificant score differences of the VaR forecasts.A is superior for CoVaR
% 3 - yellow: all insignificant
% 4 - orange: insignificant score differences of the VaR forecasts.B is superior for CoVaR
% 5 - red: A is significant inferior for VaR
if p2n < 0.05
    out = 2; %'green';
else if p2n > 0.95
        out = 4; % 'orange'
    else
        out = 3; % 'yellow';
    end
end
end

