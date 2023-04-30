%% Function that plots Thermal Power vs True Anomaly
function [f1,f2] = powerplots(Q,theta)
% Plot of Thermal Power
f1 = figure;
plot(theta,Q{1});
hold on;
plot(theta,Q{2});
plot(theta,Q{3});
plot(theta,Q{4});
plot(theta,Q{5});
xlabel('True Anomaly [deg]');
ylabel('Thermal Power [W]');
legend('Zenith Face','Nadir Face','Pos $v$ Face','Neg $v$ Face',...
    'N and S Faces','Location','north');
xlim([0,360]);

% Combining the Thermal Power
Q_comb = Q{1}+Q{2}+Q{3}+Q{4}+Q{5}+Q{6};
% Plot of Combined Thermal Power
f2 = figure;
plot(theta,Q_comb);
xlabel('True Anomaly [deg]');
ylabel('Thermal Power [W]');
legend('Combined','Location','north');
xlim([0,360]);
end