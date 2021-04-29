function[] = FinalProject()

%this is a test

global vendor;
vendor.length = 0;
vendor.mass = 0;
vendor.centerofmass = 0;
vendor.Finalanswer = 0;

vendor.fig = figure('numbertitle','off','name','Force Analysis');

vendor.LengthDisplayMessage = uicontrol('style','text','units','normalized',...
    'position', [.034 .90 .09 .095],'string','Input Length and angle in radians','horizontalalignment','right');

vendor.LengthWDisplay = uicontrol('style','text','units','normalized',...
    'position', [.034 .80 .09 .095],'string','Input length W','horizontalalignment','right');

vendor.W1 = uicontrol('style','edit','units','normalized','position',...
    [.10 .80 .09 .05],'string', 'edit','horizontalalignment','right');

vendor.LengthVDisplay = uicontrol('style','text','units','normalized',...
    'position', [.034 .70 .09 .095],'string','Input length V','horizontalalignment','right');

vendor.V1 = uicontrol('style','edit','units','normalized','position',...
    [.15 .78 .09 .05],'string', 'edit','horizontalalignment','right');

vendor.LengthUDisplay = uicontrol('style','text','units','normalized',...
    'position', [.034 .60 .09 .095],'string','Input length U','horizontalalignment','right');

vendor.U1 = uicontrol('style','edit','units','normalized','position',...
    [.15 .78 .09 .05],'string', 'edit','horizontalalignment','right');

vendor.LengthLDisplay = uicontrol('style','text','units','normalized',...
    'position', [.034 .50 .09 .095],'string','Input length L','horizontalalignment','right');

vendor.L1 = uicontrol('style','edit','units','normalized','position',...
    [.15 .78 .09 .05],'string', 'edit','horizontalalignment','right');

vendor.LengthGDisplay = uicontrol('style','text','units','normalized',...
    'position', [.034 .40 .09 .095],'string','Input length G','horizontalalignment','right');

vendor.G1 = uicontrol('style','edit','units','normalized','position',...
    [.15 .78 .09 .05],'string', 'edit','horizontalalignment','right');

vendor.Lengthm1Display = uicontrol('style','text','units','normalized',...
    'position', [.034 .80 .09 .095],'string','Input m1','horizontalalignment','right');

vendor.m1 = uicontrol('style','edit','units','normalized','position',...
    [.15 .78 .09 .05],'string', 'edit','horizontalalignment','right');

vendor.Lengthm2Display = uicontrol('style','text','units','normalized',...
    'position', [.034 .90 .09 .095],'string','Input m2','horizontalalignment','right');

vendor.m2 = uicontrol('style','edit','units','normalized','position',...
    [.15 .78 .09 .05],'string', 'edit','horizontalalignment','right');

vendor.Lengthm3Display = uicontrol('style','text','units','normalized',...
    'position', [.034 .90 .09 .095],'string','Input m3','horizontalalignment','right');

vendor.m3 = uicontrol('style','edit','units','normalized','position',...
    [.15 .78 .09 .05],'string', 'edit','horizontalalignment','right');

vendor.LengthR1Display = uicontrol('style','text','units','normalized',...
    'position', [.034 .90 .09 .095],'string','Input length R1','horizontalalignment','right');

vendor.R1 = uicontrol('style','edit','units','normalized','position',...
    [.15 .78 .09 .05],'string', 'edit','horizontalalignment','right');

vendor.LengthR2Display = uicontrol('style','text','units','normalized',...
    'position', [.034 .90 .09 .095],'string','Input length R2','horizontalalignment','right');

vendor.R2 = uicontrol('style','edit','units','normalized','position',...
    [.15 .78 .09 .05],'string', 'edit','horizontalalignment','right');

vendor.LengthR3Display = uicontrol('style','text','units','normalized',...
    'position', [.034 .90 .09 .095],'string','Input length R3','horizontalalignment','right');

vendor.R3 = uicontrol('style','edit','units','normalized','position',...
    [.15 .78 .09 .05],'string', 'edit','horizontalalignment','right');

vendor.CenterDisplayMessage = uicontrol('style','text','units','normalized',...
    'position', [.034 .90 .09 .095],'string','Input location of center of mass','horizontalalignment','right');


vendor.MassDisplayMessage = uicontrol('style','text','units','normalized',...
    'position', [.034 .90 .09 .095],'string','Input link mass','horizontalalignment','right');


vendor.solve = uicontrol('style','pushbutton','units','normalized',...
    'position',[.20 .007 .14 .05],'string','solve','callback', {@solve});

vendor.answer = uicontrol('style','text','units','normalized',...
    'position', [.034 .90 .09 .095],'string','Final Answer','horizontalalignment','right');

vendor.FinalAnswerDisplay = uicontrol('style','text','units','normalized','position',...
    [.15 .78 .09 .05],'string', num2str(vendor.Finalanswer),'horizontalalignment','right');

vendor.solve = uicontrol('style','pushbutton','units','normalized',...
    'position',[.034 .007 .14 .05],'string','Convert to Excel','callback', {@converttoExcel});

end

function [] = solve()

global vendor;

W1 = str2num(vendor.LengthWDisplay.string);
V1 = str2num(vendor.LengthVDisplay.string);
U1 = str2num(vendor.LengthUDisplay.string);
L1 = str2num(vendor.LengthLDisplay.string);
G1 = str2num(vendor.LengthGDisplay.string);
m1 = str2num(vendor.Lengthm1Display.string);
m2 = str2num(vendor.Lengthm2Display.string);
m3 = str2num(vendor.Lengthm3Display.string);
R1 = str2num(vendor.LengthR1Display.string);
R2 = str2num(vendor.LengthR2Display.string);
R3 = str2num(vendor.LengthR3Display.string);

g = -9.81;

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
     0,-m2*g-Fpy,-m2*g*R2j-Lx*Fpy,...
     0,-m3*g,-m3*g*R3j]';
 x = a\b;
 
 Ta0 = x(1);
 Fa0x = x(2);
 Fa0y = x(3);
 Fa1x = x(4);
 Fa1y = x(5);
 Fb1x = x(6);
 Fb1y = x(7);
 Fb0x = x(8);
 Fb0y = x(9); 
 
 
 
end

vendor.answer.string = vendor.Finalanswer;

end

function[] = converttoExcel()

end