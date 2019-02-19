% use (slightly) more accurate averaging.
% 300-500 eV, 2.481nm-4.136nm

n = 2.481:.01:4.136;
energies = [];
x1 = [];
x2 = [];
x3 = [];
for i = n
    output = calceff(125,3,i);
    x1 = [x1 output(2).eff4];
    x2 = [x2 output(3).eff4];
    x3 = [x3 output(4).eff4];
    energies = [energies i];
end

plot(energies,x1,energies,x2,energies,x3)
legend('1st Order','2nd Order','3rd Order')
xlabel('Wavelength(nm)')
ylabel('4th Efficiency')

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
