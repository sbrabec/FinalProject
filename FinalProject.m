function[] = FinalProject()
clc

global vendor;
vendor.W1 = 0;
vendor.V1 = 0;
vendor.L1 = 0;
vendor.U1 = 0;
vendor.G1 = 0;
vendor.m1 = 0;
vendor.m2 = 0;
vendor.m3 = 0;
vendor.R1 = 0;
vendor.R2 = 0;
vendor.R3 = 0;


vendor.fig = figure('numbertitle','off','name','Force Analysis');

vendor.LengthWDisplay = uicontrol('style','text','units','normalized',...
    'position', [.034 .80 .09 .095],'string','Input length W','horizontalalignment','right');

vendor.W1 = uicontrol('style','edit','units','normalized','position',...
    [.15 .80 .09 .05],'horizontalalignment','right');

vendor.LengthVDisplay = uicontrol('style','text','units','normalized',...
    'position', [.034 .70 .09 .095],'string','Input length V','horizontalalignment','right');

vendor.V1 = uicontrol('style','edit','units','normalized','position',...
    [.15 .70 .09 .05],'horizontalalignment','right');

vendor.LengthUDisplay = uicontrol('style','text','units','normalized',...
    'position', [.034 .60 .09 .095],'string','Input length U','horizontalalignment','right');

vendor.U1 = uicontrol('style','edit','units','normalized','position',...
    [.15 .60 .09 .05],'horizontalalignment','right');

vendor.LengthLDisplay = uicontrol('style','text','units','normalized',...
    'position', [.034 .50 .09 .095],'string','Input length L','horizontalalignment','right');

vendor.L1 = uicontrol('style','edit','units','normalized','position',...
    [.15 .50 .09 .05],'horizontalalignment','right');

vendor.LengthGDisplay = uicontrol('style','text','units','normalized',...
    'position', [.034 .40 .09 .095],'string','Input length G','horizontalalignment','right');

vendor.G1 = uicontrol('style','edit','units','normalized','position',...
    [.15 .40 .09 .05],'horizontalalignment','right');

vendor.Lengthm1Display = uicontrol('style','text','units','normalized',...
    'position', [.25 .80 .09 .095],'string','Input m1','horizontalalignment','right');

vendor.m1 = uicontrol('style','edit','units','normalized','position',...
    [.35 .80 .09 .05],'horizontalalignment','right');

vendor.Lengthm2Display = uicontrol('style','text','units','normalized',...
    'position', [.25 .70 .09 .095],'string','Input m2','horizontalalignment','right');

vendor.m2 = uicontrol('style','edit','units','normalized','position',...
    [.35 .70 .09 .05],'horizontalalignment','right');

vendor.Lengthm3Display = uicontrol('style','text','units','normalized',...
    'position', [.25 .60 .09 .095],'string','Input m3','horizontalalignment','right');

vendor.m3 = uicontrol('style','edit','units','normalized','position',...
    [.35 .60 .09 .05],'horizontalalignment','right');

vendor.LengthR1Display = uicontrol('style','text','units','normalized',...
    'position', [.25 .50 .09 .095],'string','Input length R1','horizontalalignment','right');

vendor.R1 = uicontrol('style','edit','units','normalized','position',...
    [.35 .50 .09 .05],'horizontalalignment','right');

vendor.LengthR2Display = uicontrol('style','text','units','normalized',...
    'position', [.25 .40 .09 .095],'string','Input length R2','horizontalalignment','right');

vendor.R2 = uicontrol('style','edit','units','normalized','position',...
    [.35 .40 .09 .05],'horizontalalignment','right');

vendor.LengthR3Display = uicontrol('style','text','units','normalized',...
    'position', [.25 .30 .09 .095],'string','Input length R3','horizontalalignment','right');

vendor.R3 = uicontrol('style','edit','units','normalized','position',...
    [.35 .30 .09 .05],'horizontalalignment','right');

vendor.answer = uicontrol('style','text','units','normalized',...
    'position', [.40 .90 .09 .095],'string','Final Answer','horizontalalignment','right');

vendor.solve = uicontrol('style','pushbutton','units','normalized',...
    'position',[.20 .007 .14 .05],'string','solve','callback',...
    {@solveforForce});

vendor.FinalAnswerDisplay = uicontrol('style','text','units','normalized','position',...
    [.50 0 .50 1],'horizontalalignment','right');

vendor.ExcelButton = uicontrol('style','pushbutton','units','normalized',...
    'position',[.034 .007 .14 .05],'string','Convert to Excel','callback', {@converttoExcel});

end

function [] = solveforForce(~,~)

global vendor;

W1 = str2double(vendor.W1.String);
V1 = str2double(vendor.V1.String);
U1 = str2double(vendor.U1.String);
L1 = str2double(vendor.L1.String);
G1 = str2double(vendor.G1.String);
m1 = str2double(vendor.m1.String);
m2 = str2double(vendor.m2.String);
m3 = str2double(vendor.m3.String);
R1 = str2double(vendor.R1.String);
R2 = str2double(vendor.R2.String);
R3 = str2double(vendor.R3.String);

g = -9.81;
i = sqrt(-1);

start_angle = 0;
step_angle = -0.1;
stop_angle = -35;

angular_disp = [start_angle:step_angle:stop_angle]';
N = length(angular_disp);

Fpx = 45;

theta = angle(W1);
rho = angle(V1);
sigma = angle(U1);
delta = angle(L1);

Gx = real(G1);
Gy = imag(G1);

W = abs(W1);
V = abs(V1);
U = abs(U1);
L = abs(L1);
G = abs(G1);

f_alpha_gamma = @(x,B) [...
    W*cos(theta+B)+V*cos(rho+x(1))-U*cos(sigma+x(2))-Gx;
    W*sin(theta+B)+V*sin(rho+x(1))-U*sin(sigma+x(2))-Gy];


R1y = imag(R1);
R1x = real(R1);
R2y = imag(R2);
R2x = real(R2);
R3y = imag(R3);
R3x = real(R3); 

Ta0 = ones(N,1)*Inf;

Fa0 = ones(N,2)*Inf;
Fa1 = ones(N,2)*Inf;
Fb0 = ones(N,2)*Inf;
Fb1 = ones(N,2)*Inf;

tau = ones(N,1)*Inf;

options = optimoptions('fsolve','Display','off');

x0 = [0,0];

for j=1:N
    B = angular_disp(j) * pi / 180;
    [x] = fsolve(f_alpha_gamma, x0, options, B);
   
    
    A = x(1);
    G = x(2);
    
    x0 = x;
   
    tau(j) = abs((sigma + G - (rho + A)) * 180 / pi);
    
    if tau(j) > 180
        tau(j) = 360 - tau(j);
    end
    
Wy=W*sin(theta+B);
Wx=W*cos(theta+B);
Vy=V*sin(rho+A);
Vx=V*cos(rho+A);
Uy=U*sin(sigma+G);
Ux=U*cos(sigma+G);
Ly=L*sin(delta+A);
Lx=L*cos(delta+A);


R1j = R1x*cos(B)-R1y*sin(B);
R2j = R2x*cos(A)-R2y*sin(A);
R3j = R3x*cos(G)-R3y*sin(G);  
    
a = [0,1,0,1,0,0,0,0,0;
     0,0,1,0,1,0,0,0,0;
     1,0,0,-Wy,Wx,0,0,0,0;
     0,0,0,-1,0,1,0,0,0;
     0,0,0,0,-1,0,1,0,0;
     0,0,0,0,0,-Vy,Vx,0,0;
     0,0,0,0,0,-1,0,1,0;
     0,0,0,0,0,0,-1,0,1;
     0,0,0,0,0,-Uy,Ux,0,0;];
 b = [0,-m1*g,-m1*g*R1j,...
     0,-m2*g-Fpx,-m2*g*R2j-Lx*Fpx,...
     0,-m3*g,-m3*g*R3j]';
 x = a\b;
 
 vendor.Ta0 = (x(1));
 vendor.Fa0x = (x(2));
 vendor.Fa0y = (x(3));
 vendor.Fa1x = (x(4));
 vendor.Fa1y = (x(5));
 vendor.Fb1x = (x(6));
 vendor.Fb1y = (x(7));
 vendor.Fb0x = (x(8));
 vendor.Fb0y = (x(9)); 
 
 
end

output = sprintf('Ta0 = %s \nFa0x = %s \nFa0y = %s \nFa1x = %s \nFa1y = %s \nFb1x = %s \nFb1y = %s\nFb0x = %s \nFb0y = %s', ...
        vendor.Ta0, vendor.Fa0x, vendor.Fa0y, vendor.Fa1x, vendor.Fa1y, vendor.Fb1x, vendor.Fb1y, vendor.Fb0x, vendor.Fb0y);
    
displayOutput(output)

end

function[] = displayOutput(output)

global vendor
vendor.FinalAnswerDisplay.String = output;

end

function[] = converttoExcel(~,~)

solveforForce()
global vendor

filename = 'ForceAnalysis.csv';
fid = fopen(filename,'w');
fprintf(fid,'%s\n',[num2str(vendor.Ta0),',', num2str(vendor.Fa0x),',', num2str(vendor.Fa0y),',', ...
    num2str(vendor.Fa1x),',', num2str(vendor.Fa1y),',', num2str(vendor.Fb1x),',', ...
    num2str(vendor.Fb1y),',', num2str(vendor.Fb0x),',' num2str(vendor.Fb0y),',']);
fclose(fid);
end