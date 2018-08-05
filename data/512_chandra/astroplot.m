addpath(genpath('C:\Program Files\MATLAB\R2018a\toolbox\Wavelab850'))
addpath(genpath('D:\Philippe\SwitchDrive\AstroQUT_matlab\ASTROQUT v2.0'))

fid = fopen("MATLAB/1000/sol.data", "r");
sol = fread(fid, 263168, 'double');
fclose(fid);
fid = fopen("MATLAB/1000/divx.data", "r");
divx = fread(fid, 263168, 'double');
fclose(fid);
fid = fopen("AstroQUT/5000/computed_sol.data", "r");
sol_c = fread(fid, 263168, 'double');
fclose(fid);
fid = fopen("AstroQUT/5000/computed_divx.data", "r");
divx_c = fread(fid, 263168, 'double');
fclose(fid);

M = 512;
MM = 1024;
N = 262144;
R = 256;

sol_c = sol_c ./ divx_c;
flassow = sol(1:M);
flassow_c = sol_c(1:M);
flassos = sol(M+1:MM);
flassos_c = sol_c(M+1:MM);

[~,~,qmf,S]=setUpOperatorsWS(M,N,R,0,0,0,'Daubechies',8,true,1);

flasso = IWT_PO(flassow, 0, qmf) + S*flassos;
flasso_c = IWT_PO(flassow_c, 0, qmf) + S*flassos_c;

x=linspace(-M/2,M/2,M)*R/(M/2);
x=x((M/2+1):M);

figure
loglog(x,flasso((M/2+1):M), 'color', 'blue', 'LineWidth', 2, 'DisplayName', 'MATLAB')
hold on
loglog(x,flasso_c((M/2+1):M), 'color', 'red','LineStyle', ':', 'LineWidth', 2, 'DisplayName', 'C++')
hold off

figure
loglog(x,flasso((M/2):-1:1), 'color', 'blue', 'LineWidth', 2, 'DisplayName', 'MATLAB')
hold on
loglog(x,flasso_c((M/2):-1:1), 'color', 'red', 'LineStyle', ':', 'LineWidth', 2, 'DisplayName', 'C++')
hold off
