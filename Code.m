close all
clear
clc

%% Read image

addpath(genpath('images'));
Im = im2double(imread('images/Image2'));

if ndims(Im) == 3
    Im = rgb2gray(Im);
    
end
    
figure(1), imshow(Im),title('original image');

%% Edge detection

BW1=edge(Im,'prewitt');
BW2=edge(Im,'canny');
BW3=edge(Im,'sobel');
figure(2),
subplot(2,2,1),imshow(Im), title('Original','FontSize',20); 
subplot(2,2,2),imshow(BW1), title('prewitt','FontSize',20); 
subplot(2,2,3),imshow(BW2), title('canny','FontSize',20);
subplot(2,2,4),imshow(BW3), title('sobel','FontSize',20);

%% Cornor detection

dx = [-1 0 1; -1 0 1; -1 0 1];   % Derivative masks
dy = dx';

Ix = conv2(Im, dx, 'same');      % Image derivatives
Iy = conv2(Im, dy, 'same');

% set the parameter for Gaussian convolution used in Harris Corner Detector

SIGMA_gaussian=1.7;
g = fspecial('gaussian',max(1,fix(3*SIGMA_gaussian)+1), SIGMA_gaussian);

Ix2 = conv2(Ix.^2, g, 'same'); % Smoothed squared image derivatives
Iy2 = conv2(Iy.^2, g, 'same');
Ixy = conv2(Ix.*Iy, g, 'same');

% cim = det(M) - k trace(M)^2.

k = 0.03;
cim = (Ix2.*Iy2 - Ixy.^2) - k * (Ix2 + Iy2);

% remove boundaries of cim which is going to have large values because of zero padding of the image

BORDER=100;
cim(1:BORDER,:)=0;
cim(end-BORDER:end,:)=0;
cim(:,end-BORDER:end)=0;
cim(:,1:BORDER)=0;

% Thresholding the cim

T=mean(cim(:));
CIM=cim;
CIM(find(cim<T))=0;

% perform nonlocal maximum suppression on the thresholded measure

support=true(40);

% compute maximum over a square neighbor of size 40 x 40

maxima=ordfilt2(CIM,sum(support(:)),support);

% determine the locations where the max over the neigh or 40 x 40 corresponds to the cim values

[loc_x,loc_y]=find((cim==maxima).*(CIM>0));
indx = find((cim==maxima).*(CIM>0));

% draw a cross on the image in the local maxima
figure(3), imshow(Im,[]), hold on,
plot(loc_y,loc_x,'g+', 'LineWidth', 4)
title('Cornor Detection')

%% Conics

% Manually extract feature points on the two ellipse
xl = [2388;2395;2331;2234;2209];
yl = [1767;2232;2278;2265;2131];

xr = [3697;3703;3563;3578;3611];
yr = [917;1060;1193;1122;1026];


% Calculating conic matrix
Al=[xl.^2 xl.*yl yl.^2 xl yl ones(size(xl))];
Ar=[xr.^2 xr.*yr yr.^2 xr yr ones(size(xr))];
Nl = null(Al);
Nr = null(Ar);
ccl = Nl(:, 1);
ccr = Nr(:, 1);
[al, bl, cl, dl, el, fl]=deal(ccl(1),ccl(2),ccl(3),ccl(4),ccl(5),ccl(6));
[ar, br, cr, dr, er, fr]=deal(ccr(1),ccr(2),ccr(3),ccr(4),ccr(5),ccr(6));

% here is the matrix of the conic
Cl=[al bl/2 dl/2; bl/2 cl el/2; dl/2 el/2 fl];
Cr=[ar br/2 dr/2; br/2 cr er/2; dr/2 er/2 fr];

% Generate the conic on image
iml=zeros(3096,4128);
for i=1:3096
    for j=1:4128
        iml(i,j)=[j i 1]*Cl*[j i 1]';
    end
end

imr=zeros(3096,4128);
for i=1:3096
    for j=1:4128
        imr(i,j)=[j i 1]*Cr*[j i 1]';
    end
end

figure(4)
imshow(Im .* (iml > 0 & imr > 0));
title("Conics Analyse");
hold on;
scatter(xl,yl, 100, 'filled');
scatter(xr,yr, 100, 'filled');
% Calculate center of conics

Cenl = [al bl/2; bl/2 cl] \ [-dl/2; -el/2];
Cenr = [ar br/2; br/2 cr] \ [-dr/2; -er/2];
Cenl = [Cenl; 1];
Cenr = [Cenr; 1];
[xcl, ycl, xcr, ycr] = deal(Cenl(1), Cenl(2), Cenr(1), Cenr(2));

% Connect the two centers, calculate the intersections of the conics

Line = cross(Cenl, Cenr);

res = zeros([1, 2000]);
minl = 100;
minr = 100;
for xx = xcl:xcr
    yy = -(Line(1) * xx + Line(3))/Line(2);
    p = [xx;yy;1];    
    if (abs(p.' * Cl * p) < minl)
        minl = abs(p.' * Cl * p);
        P_1_x = xx;
    end
    yy = -(Line(1) * xx + Line(3))/Line(2);
    p = [xx;yy;1];
    if (abs(p.' * Cr * p) < minr)
        minr = abs(p.' * Cr * p);
        P_2_x = xx;
    end
end

P_1_y = fix(-(Line(1) * P_1_x + Line(3))/Line(2));
P_2_y = fix(-(Line(1) * P_2_x + Line(3))/Line(2));


%% Calculate the ratio between diameter and distance

FNT_SZ = 20;
figure(5)
imshow(Im .* (iml > 0 & imr > 0));
title("Conics Analyse");
hold on;
scatter(xcl,ycl, 100, 'filled');
text(xcl, ycl+50, 'E_c_l', 'FontSize', FNT_SZ, 'Color', 'w')
scatter(xcr,ycr, 100, 'filled');
text(xcr, ycr+50, 'E_c_r', 'FontSize', FNT_SZ, 'Color', 'w')
scatter(P_1_x,P_1_y, 100, 'filled');
text(P_1_x, P_1_y+50, 'P_1', 'FontSize', FNT_SZ, 'Color', 'w')
scatter(P_2_x,P_2_y, 100, 'filled');
text(P_2_x, P_2_y+50, 'P_2', 'FontSize', FNT_SZ, 'Color', 'w')
plot([xcl, xcr], [ycl, ycr], 'LineWidth',2);

Y = [xcl, ycl];
Z = [P_1_x, P_1_y];
X1 = [P_2_x, P_2_y];
X2 = [xcr, ycr];

CR = (norm(X1-Y)/norm(X1-Z))/(norm(X2-Y)/norm(X2-Z));
Ratio = ((2 * CR - 2) + sqrt((2 * CR - 2) ^ 2 - 4 * (1 - CR)));
% Here is the ratio
Ratio = 1 / ((1 / Ratio) + 1)

%% Determine K

% Manually extract some parallel feature point pairs
[X11, Y11, X12, Y12, X13, Y13, X14, Y14, X21, Y21, X22, Y22, X23, Y23, X24, Y24, X31, Y31, X32, Y32, X33, Y33, X34, Y34] = deal(1048, 1942, 734, 1725, 1174, 1530, 641, 1231, 2378, 2019, 3627, 1102, 2415, 1744, 3682, 910, 1048, 1948, 1060, 1838, 732, 1725, 746, 1625);

% Calculate the coordinates of vanishing points
[C1, C2, C3] = deal(cross(cross([X11, Y11, 1], [X12, Y12, 1]), cross([X13, Y13, 1], [X14, Y14, 1])), cross(cross([X21, Y21, 1], [X22, Y22, 1]), cross([X23, Y23, 1], [X24, Y24, 1])), cross(cross([X31, Y31, 1], [X32, Y32, 1]), cross([X33, Y33, 1], [X34, Y34, 1])));
[x1, y1, x2, y2, x3, y3] = deal(C1(1)/C1(3), C1(2)/C1(3), C2(1)/C2(3), C2(2)/C2(3), C3(1)/C3(3), C3(2)/C3(3));

figure(6),
imshow(Im),
hold on;
scatter(x1,y1, 100, 'filled');
scatter(x2,y2, 100, 'filled');
scatter(x3,y3, 100, 'filled');

plot([X11,x1],[Y11,y1],'LineWidth',2);
plot([X13,x1],[Y13,y1],'LineWidth',2);
plot([X21,x2],[Y21,y2],'LineWidth',2);
plot([X23,x2],[Y23,y2],'LineWidth',2);
plot([X31,x3],[Y31,y3],'LineWidth',2);
plot([X33,x3],[Y33,y3],'LineWidth',2);
text(x1, y1+250, 'V_x', 'FontSize', FNT_SZ, 'Color', 'b')
text(x2, y2+250, 'V_y', 'FontSize', FNT_SZ, 'Color', 'b')
text(x3, y3+250, 'V_z', 'FontSize', FNT_SZ, 'Color', 'b')
hold off; 
%%
% Calculate K by vanishing points
[vx, vy, vz] = deal([x1; y1; 1], [x2; y2; 1], [x3; y3; 1]);

syms fa fb px py
K = [fa 0 px; 0 fb py; 0 0 1];
w = inv(K*K');
eq1 = fa / fb == 4128 / 3096;
eq2 = vx' * w * vy == 0;
eq3 = vy' * w * vz == 0;
eq4 = vz' * w * vx == 0;

rs = solve([eq1, eq2, eq3, eq4]);

% Here is the matrix K

ans_K = double([rs.fa 0 rs.px; 0 rs.fb rs.py; 0 0 1]);
ans_K

%% 3D reconstruction

% Calculating the image coordinate of the world coordinate's origin
[y, w1, w2] = deal([X11; Y11; 1], [X12; Y12; 1], [x1; y1; 1]);
k = 2 * norm(w2 - y) / norm(w1 - y);
[xo, yo] = deal((w2(1) - k * w1(1)) / (1-k), (w2(2) - k * w1(2)) / (1-k));

figure(7),
imshow(Im),
hold on;
scatter(x1,y1, 100, 'filled');
scatter(x2,y2, 100, 'filled');
scatter(x3,y3, 100, 'filled');
scatter(xo,yo, 20, 'filled');
plot([xo,x1],[yo,y1]);
plot([xo,x2],[yo,y2]);
plot([xo,x3],[yo,y3]);
text(x1, y1+250, 'V_x', 'FontSize', FNT_SZ, 'Color', 'b')
text(x2, y2+250, 'V_y', 'FontSize', FNT_SZ, 'Color', 'b')
text(x3, y3+250, 'V_z', 'FontSize', FNT_SZ, 'Color', 'b')
text(xo, yo+250, 'O_I', 'FontSize', FNT_SZ, 'Color', 'w')
hold off; 

% Construct P matrix

P = [x1 x2 x3 xo; y1 y2 y3 yo; 1 1 1 1];
[M, m] = deal(P(:, 1:3), P(:, 4));


% Calculate the 3 axes-units image coordinate

UP = P * [1 0 0; 0 1 0; 0 0 1; 1 1 1];
[Ux, Uy, Uz] = deal(UP(:, 1), UP(:, 2), UP(:, 3));
Ux = Ux / Ux(3);
Uy = Uy / Uy(3);
Uz = Uz / Uz(3);

figure(8),
imshow(Im),
hold on;
scatter(x1,y1, 50, 'filled', 'MarkerEdgeColor', 'k');
scatter(xo,yo, 50, 'filled', 'MarkerEdgeColor', 'w');
scatter(X11,Y11, 50, 'filled', 'MarkerEdgeColor', 'w');
scatter(X12,Y12, 50, 'filled', 'MarkerEdgeColor', 'w');
scatter(Ux(1),Ux(2), 50, 'filled', 'MarkerEdgeColor', 'k');
plot([xo,x1],[yo,y1]);
text(x1, y1+250, 'V_x', 'FontSize', FNT_SZ, 'Color', 'b')
text(Ux(1)-150, Ux(2)+150, 'U_x', 'FontSize', FNT_SZ, 'Color', 'b')
text(xo-150, yo+150, 'O_I', 'FontSize', FNT_SZ, 'Color', 'w')
text(X11-50, Y11+150, 'X_1', 'FontSize', FNT_SZ, 'Color', 'w')
text(X12-250, Y12+50, 'X_2', 'FontSize', FNT_SZ, 'Color', 'w')
hold off; 

% Calculate world coordinate of X1 and X2

[p1, p2] = deal([X11; Y11; 1], [X12; Y12; 1]);
[CR1, CR2] = deal((norm(Ux - p1)/norm(Ux - m)) / (norm(vx - p1) / norm(vx - m)), (norm(Ux - m)/norm(Ux - p2)) / (norm(vx - m) / norm(vx - p2)));

% Here is X1's world coordinate
X_1_w = [1-CR1; 0; 0; 1]
% Here is X2's world coordinate
X_2_w = [1-1/CR2; 0; 0; 1]

% Evaluate
check_p1 = P*X_1_w;
check_p1 = check_p1/check_p1(3);
check_p2 = P*X_2_w;
check_p2 = check_p2/check_p2(3);


%% Determine the camera position
O = -M\m