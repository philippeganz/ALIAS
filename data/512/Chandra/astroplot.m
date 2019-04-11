load 'D:\Philippe\SwitchDrive\AstroQUT_matlab\ASTROQUT\data_and_results\FEOchandra.mat'
load 'D:\Philippe\GoogleDrive\Documents\Cours\Thesis\AstroQUT\data\512_chandra\MATLAB\2000\ChandraFEO.mat'

M = 512;
MM = 1024;
N = 262144;
R = 256;

fhat = sol.fhat;
fid = fopen("computed_fhat.data", "r");
fhat_c = fread(fid, M, 'double');
fclose(fid);

figure
semilogy(fhat((R+1):round(M-R*(1-1/sqrt(2)))), 'color', 'blue', 'LineWidth', 2, 'DisplayName', 'MATLAB')
hold on
semilogy(fhat_c((R+1):round(M-R*(1-1/sqrt(2)))), 'color', 'red','LineStyle', ':', 'LineWidth', 2, 'DisplayName', 'C++')
hold off

figure
semilogy(fhat(R:-1:round(R*(1-1/sqrt(2)))), 'color', 'blue', 'LineWidth', 2, 'DisplayName', 'MATLAB')
hold on
semilogy(fhat_c(R:-1:round(R*(1-1/sqrt(2)))), 'color', 'red', 'LineStyle', ':', 'LineWidth', 2, 'DisplayName', 'C++')
hold off

figure
semilogy(fhat(round(R*(1-1/sqrt(2)):R+round(R*(1/sqrt(2))))), 'color', 'blue', 'LineWidth', 2, 'DisplayName', 'MATLAB')
hold on
semilogy(fhat_c(round(R*(1-1/sqrt(2)):R+round(R*(1/sqrt(2))))), 'color', 'red', 'LineStyle', ':', 'LineWidth', 2, 'DisplayName', 'C++')
hold off

dhat = sol.dhatpos;
fid = fopen("computed_dhat.data", "r");
dhat_c = fread(fid, N, 'double');
fclose(fid);
fid = fopen("computed_dhat_CC.data", "r");
dhat_CC = fread(fid, N, 'double');
fclose(fid);

figure
final = (Freal-Oreal)./Ereal;
imshow(100*final/max(max(final)))

figure
[X,Y] = meshgrid(1:M,1:M);
imshow(100*final/max(max(final)))
hold on
plot3(X(dhat>0),Y(dhat>0),dhat(dhat>0),'ro')
hold off

figure
imshow(100*final/max(max(final)))
hold on
plot3(X(dhat_c>0),Y(dhat_c>0),dhat_c(dhat_c>0),'ro')
hold off

figure
imshow(100*final/max(max(final)))
hold on
plot3(X(dhat_CC>0),Y(dhat_CC>0),dhat_CC(dhat_CC>0),'ro')
hold off
