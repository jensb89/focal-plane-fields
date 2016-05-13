%% Focal fields - TEM00 Gaussian Mode
% This script calculates the field of a TEM00 Gaussian Mode in the focus of
% a microscope objective with a given NA. The formula used here was 
% taken from "Novotny, L. and Hecht, B. (2006); Principles of Nano-Optics".
%
% © Jens Brauer, 2016

%% Constants
lambda = 800;      % Wavelength [nm]
k = 2*pi/lambda;   % Wave vector [1/nm]
w0 = 5E6;          % Beam waist [nm]
n1 = 1;            % refractive index air
n2 = 1.518;        % refractive index glass (MO)
E_0 = 1;           % E_0
z = 0;             % deviation from focal plane
NA = 1.4;          % NA
f0 = 1;            % filling factor (1 means the MO is fully illuminated)
numpoints = 1000;  % number of points for the integral


%% theta_max
theta_max = asin(NA/n2);
theta_max_deg = rad2deg(theta_max);
theta = linspace(0,theta_max,numpoints);
f = w0/(f0*sin(theta_max));

%% Calculations
x = linspace(-lambda,lambda,200);
y = linspace(-lambda,lambda,200);
Ex = zeros(1,length(y)*length(x));
Ey = Ex;
Ez = Ex;
xy = meshgrid(x,y);
y_vals = xy(1,:);
tic;
parfor i=1:numel(xy)
        yp = xy(1,mod(i-1,200)+1);
        xp = xy(i);
        % rho
        rho = sqrt(xp.^2 + yp.^2);
        %phi
        if xp > 0
            phi = atan(yp/xp);
        elseif xp < 0 && yp>=0
            phi = atan(yp/xp) + pi;
        elseif xp < 0 && yp<0
            phi = atan(yp/xp) + pi;
        elseif xp==0 && yp>0
            phi = pi/2;
        elseif xp==0 && yp<0
            phi = -pi/2;
        end

        % fw 
        fw = exp(-1/f0^2 * sin(theta).^2/sin(theta_max)^2);

        % integrals
        I00 = sum(fw.*sqrt(cos(theta)).*sin(theta).*(1+cos(theta)).*besselj(0,k*rho*sin(theta)).*exp(1i*k*z*cos(theta)))/numpoints;
        I01 = sum(fw.*sqrt(cos(theta)).*sin(theta).^2.*besselj(1,k*rho*sin(theta)).*exp(1i*k*z*cos(theta)))/numpoints;
        I02 = sum(fw.*sqrt(cos(theta)).*sin(theta).*(1-cos(theta)).*besselj(2,k*rho*sin(theta)).*exp(1i*k*z*cos(theta)))/numpoints;
        
        % fields
        Ex(i) = 1i*k*f/2 * sqrt(n1/n2) * E_0 * exp(1i*k*f) * (I00 + I02*cos(2*phi));
        Ey(i) = 1i*k*f/2 * sqrt(n1/n2) * E_0 * exp(1i*k*f) * (I02*sin(2*phi));
        Ez(i) = 1i*k*f/2 * sqrt(n1/n2) * E_0 * exp(1i*k*f) * (-2*1i*I01*cos(phi));

end

Ex = reshape(Ex,200,200);
Ey = reshape(Ey,200,200);
Ez = reshape(Ez,200,200);

time = toc;

%% PLotting
fig=figure;
set(fig,'Position',[191 532 1639 417])

subplot(1,3,1);
im1 = imagesc(x,y,abs(Ex).^2);
title('|E_x|²')
max_int1 = max(max(get(im1,'CData')));


subplot(1,3,2);
im2 = imagesc(x,y,abs(Ey).^2);
title('|E_y|²')
max_int2 = max(max(get(im2,'CData')));
text(500,700,sprintf('x%5.0f',max_int1/max_int2),'Color','w','FontSize',12)


subplot(1,3,3);
im3 = imagesc(x,y,abs(Ez).^2);
title('|E_z|²')
max_int3 = max(max(get(im3,'CData')));
text(500,700,sprintf('x%5.0f',max_int1/max_int3),'Color','w','FontSize',12)

%export_fig(['Ex_Ey_Ez_NA_',num2str(NA),'_z',num2str(z),'.png'],'-dpng','-r300','-transparent')