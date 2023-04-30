%% Function that plots Electrical Power vs True Anomaly
function [f1,f2] = powerplots(P,theta,deg_rate)
% Plot of Electrical Power
f1 = figure;
plot(theta,P{1});
hold on;
plot(theta,P{2});
plot(theta,P{3});
plot(theta,P{4});
plot(theta,P{5});
xlabel('True Anomaly [deg]');
ylabel('Power Generated [W]');
legend('Zenith Face','Nadir Face','Pos $v$ Face','Neg $v$ Face',...
    'N and S Faces','Location','north');
xlim([0,360]);

% Combining the Electrical Power
P_comb = P{1}+P{2}+P{3}+P{4}+P{5}+P{6};
% Plot of Combined Thermal Power
f2 = figure;
plot(theta,P_comb);
hold on;
plot(theta,P_comb.*(1-deg_rate));
xlabel('True Anomaly [deg]');
ylabel('Power Generated [W]');
legend('Combined','With 2.5\% Degradation','Location','north');
xlim([0,360]);
end