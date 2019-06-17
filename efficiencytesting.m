% use (slightly) more accurate averaging.
% 300-500 eV, 2.481nm-4.136nm

% Find out how one order changes with wavelength and number of strata in a
% single plot, then compare it to existing values
n = 100:20:1200;
waves = 1239.84193 ./ n;
x0 = zeros(1,length(n));
x1 = zeros(1,length(n));
x2 = zeros(1,length(n));
x3 = zeros(1,length(n));
x4 = zeros(1,length(n));
x5 = zeros(1,length(n));
xtot = zeros(1,length(n));
for i = 1:length(n)
    output = calceff(125,10,waves(i));
    for j=1:length(output)
        if output(j).m1 == 0
            x0(i) = output(j).eff4;
        elseif output(j).m1 == -1
            x1(i) = output(j).eff4;
        elseif output(j).m1 == -2
            x2(i) = output(j).eff4;
        elseif output(j).m1 == -3
            x3(i) = output(j).eff4;
        elseif output(j).m1 == -4
            x4(i) = output(j).eff4;
        elseif output(j).m1 == -5
            x5(i) = output(j).eff4;
        end
        
    end
    xtot(i) = x1(i)+x2(i)+x3(i)+x4(i)+x5(i);
end
plot(n,xtot,n,x0,n,x1,n,x2,n,x3,n,x4,n,x5)
legend('Total','0th Order','1st Order','2nd Order','3rd Order','4th Order','5th Order')
xlabel('Energy (eV)')
ylabel('4th Efficiency')
xlim([100 1200])

% Need to plot efficiency as a function of energy on the defined range
% Plot multiple orders (0,1,2,3)


% n = 0;
% ns = [];
% eff1 = [];
% eff2 = [];
% eff3 = [];
% eff4 = [];
% while n < 150
%     output = calceff(n,0,2.47);
%     ns = [ns n];
%     eff1 = [eff1 output.eff1];
%     eff2 = [eff2 output.eff2];
%     eff3 = [eff3 output.eff3];
%     eff4 = [eff4 output.eff4];
%     fprintf('n = %d\n',n);
%     n = n+1;
% end
%    
% plot(ns,eff1,ns,eff2,ns,eff3,ns,eff4);
% xlabel('Number of Blocks')
% ylabel('eff4 Efficiency')
% legend('Eff1','Eff2','Eff3','Eff4')
