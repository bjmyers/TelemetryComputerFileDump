% use (slightly) more accurate averaging.
% 300-500 eV, 2.481nm-4.136nm

% Find out how one order changes with wavelength and number of strata in a
% single plot, then compare it to existing values
n = 2:.1:10;
energies = zeros(length(n));
x1 = zeros(1,length(n));
x2 = zeros(1,length(n));
x3 = zeros(1,length(n));
for i = 1:length(n)
    output = calceff(25,3,n(i));
    try
        x1(i) = output(2).eff4;
    catch
        x1(i) = 0;
    end
    try
        x2(i) = output(3).eff4;
    catch
        x2(i) = 0;
    end
    try
        x3(i) = output(4).eff4;
    catch
        x3(i) = 0;
    end
    energies(i) = n(i);
end
% sum = x1 + x2 + x3;
plot(energies,x1,energies,x2,energies,x3)
legend('1st Order','2nd Order','3rd Order')
xlabel('Wavelength(nm)')
ylabel('4th Efficiency')
xlim([2.481 10])

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
