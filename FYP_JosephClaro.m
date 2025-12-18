% FYP - Joseph Claro - 21317623

clf
%% Variables
U = 1;                       % Linear flow velocity
a = pi/8;                    % Linear flow angle
gam = 6;                     % Vortex circulation
b = 1;                       % Joukowsky parameter/blade lengthscale

C = 10;                      % Streamfunction = constant (+/- range)

center = [-0.15;0.4];         % Pre-transform circle center - circle intersects (1,0)
R = sqrt((center(1)-b)^2+(center(2)^2));

N = 20;                      % no. of grid lines
n = 200;                     % no. of points per grid line
circ_n = 500;                % no. of points on circle/airfoil

gridcolor = [0.5 0.5 0.9];   % Grid line color
gridWidth = 4;               % Grid line Thickness

titleSize = 16;
subtitleSize = 13;

plt = 1;                    % Plot indexing variable
%% Circle/Joukowsky functions

jouk = @(x,y)[x+(b^2.*x)./(x.^2+y.^2);      % Joukowski Transform
              y-(b^2.*y)./(x.^2+y.^2)]; 

joukX = @(x,y)x+(b^2.*x)./(x.^2+y.^2);      % Transform, split into a function per axis
joukY = @(x,y)y-(b^2.*y)./(x.^2+y.^2);

% defines pre-image circle
circFunc = @(t)[R*cos(t) + center(1);       % Parametric func. for circle
                R*sin(t) + center(2)];

circMesh = linspace(0,2*pi,circ_n);
circ = circFunc(circMesh);                  % points on circle

airf = jouk(circ(1,:),circ(2,:));           % transformed circle

plotLeft = center(1)-3*R;                   % plotting boundaries
plotRight = center(1)+3*R;
plotBottom = center(2)-3*R;
plotTop = center(2)+3*R;

%% Flow Streamfunction

streamf_linear = @(x,y) U.*(y.*cos(a)-x.*sin(a)-((R^2).*((y-center(2)).*cos(a)-(x-center(1)).*sin(a)))./((x-center(1)).^2+(y-center(2)).^2));
streamf_vortex = @(x,y) (gam/(2.*pi)).*log((x-center(1)).^2+(y-center(2)).^2);

streamfunc = @(x,y,C) streamf_linear(x,y)+streamf_vortex(x,y)-C;

isOutside = @(x,y) ((x-center(1)).^2 + (y-center(2)).^2 > 1.5.*R); % Boolean to check if a point is inside the input circle

plotStreamf = @(x,y,C) isOutside(x,y).*streamfunc(x,y,C);

%% Grid lines
% VERTICAL LINES
vertLines = linspace(plotLeft, plotRight, N);

vertGridPts = ones(n,N,2);

vertDiscr = linspace(plotTop, plotBottom, n)';
horizDiscr = linspace(plotLeft, plotRight, n);

for i = 1:n
    vertGridPts(i,:,1) = vertGridPts(i,:,1).*vertLines;
end
for i = 1:N
    vertGridPts(:,i,2) = vertGridPts(:,i,2).*vertDiscr;
end


% HORIZONTAL LINES
horizLines = linspace(plotTop, plotBottom, N)';

horizGridPts = ones(N,n,2);

for i = 1:N
    horizGridPts(i,:,1) = horizGridPts(i,:,1).*horizDiscr;
end
for i = 1:n
    horizGridPts(:,i,2) = horizGridPts(:,i,2).*horizLines;
end

% Remove points located inside circle (-> inside airfoil)
horizGridPts_ = removePointsInside(horizGridPts,center,R);
vertGridPts_ = removePointsInside(vertGridPts,center,R);
%% Joukowsky Grid Lines
% Initialise transformed grid lines
joukHoriz = zeros(size(horizGridPts_));
joukVert = zeros(size(vertGridPts_));

joukHoriz_ = joukHoriz;
joukVert_ = joukVert;

% Apply Joukowsky Transform to entire grid
joukHoriz(:,:,1) = joukX(horizGridPts(:,:,1),horizGridPts(:,:,2));
joukHoriz(:,:,2) = joukY(horizGridPts(:,:,1),horizGridPts(:,:,2));
joukVert(:,:,1) = joukX(vertGridPts(:,:,1),vertGridPts(:,:,2));
joukVert(:,:,2) = joukY(vertGridPts(:,:,1),vertGridPts(:,:,2));

% Apply Joukowsky Transform to grid - inner points removed
joukHoriz_(:,:,1) = joukX(horizGridPts_(:,:,1),horizGridPts_(:,:,2));
joukHoriz_(:,:,2) = joukY(horizGridPts_(:,:,1),horizGridPts_(:,:,2));
joukVert_(:,:,1) = joukX(vertGridPts_(:,:,1),vertGridPts_(:,:,2));
joukVert_(:,:,2) = joukY(vertGridPts_(:,:,1),vertGridPts_(:,:,2));

%% Plot Entire grid - View Singularity


figure(plt)
clf

% Pre-Transform, grid only
sgtitle(sprintf("Joukowsky Transform, $b = %.2f$",b),"fontsize",titleSize, "interpreter","latex");
subplot(1,2,1);
hold on;

plot(vertGridPts(:,:,1), vertGridPts(:,:,2),".","MarkerSize",gridWidth,"Color",gridcolor);
plot(horizGridPts(:,:,1), horizGridPts(:,:,2),".","MarkerSize",gridWidth,"Color",gridcolor);

xlim([2*plotLeft 2*plotRight]);
ylim([2*plotBottom 2*plotTop]);
xline(0,"linewidth",1); yline(0,"linewidth",1);

subtitle(sprintf("Pre-Transform"), "fontsize",subtitleSize,"interpreter","latex");

% Post-Transform, grid only
subplot(1,2,2)
hold on;

plot(joukVert(:,:,1), joukVert(:,:,2),".","MarkerSize",gridWidth,"Color",gridcolor);
plot(joukHoriz(:,:,1)', joukHoriz(:,:,2)',".","MarkerSize",gridWidth,"Color",gridcolor);

xlim([2*plotLeft 2*plotRight]);
ylim([2*plotBottom 2*plotTop]);
xline(0,"linewidth",1); yline(0,"linewidth",1);

subtitle(sprintf("Post-Transform"), "fontsize",subtitleSize,"interpreter","latex");
shg
plt = plt + 1;





%% plot Circle/Airfoil on Grid


figure(plt);
clf

sgtitle(sprintf("Joukowsky Airfoil, $(x_0,y_0) = (%.2f, %.2f), b = %.2f$",center(1),center(2),b),"fontsize",titleSize, "interpreter","latex");
subplot(1,2,1);
hold on;   

% Pre-Transform Circle/Grid
plot(vertGridPts_(:,:,1), vertGridPts_(:,:,2),".","MarkerSize",gridWidth,"Color",gridcolor);
plot(horizGridPts_(:,:,1)', horizGridPts_(:,:,2)',".","MarkerSize",gridWidth,"Color",gridcolor);
plot(circ(1,:),circ(2,:),"r-","linewidth",1.3);

xlim([plotLeft plotRight]);
ylim([plotBottom plotTop]);
xline(0,"linewidth",1); yline(0,"linewidth",1);

subtitle(sprintf("Pre-Transform"), "fontsize",subtitleSize,"interpreter","latex");

% Transformed circle/grid
subplot(1,2,2)
hold on;

plot(joukVert_(:,:,1), joukVert_(:,:,2),".","markersize",gridWidth, "Color",gridcolor);
plot(joukHoriz_(:,:,1)', joukHoriz_(:,:,2)',".","markersize",gridWidth, "Color",gridcolor);
plot(airf(1,:),airf(2,:),"r-","linewidth",1.5);

xlim([plotLeft plotRight]);
ylim([plotBottom plotTop]);
xline(0,"linewidth",1); yline(0,"linewidth",1);

subtitle(sprintf("Post-Transform"),"fontsize",subtitleSize, "interpreter","latex");

shg
plt = plt + 1;      % Next plot

%% Plot Streamlines - Circle/Airfoil
figure(plt)
clf

streamLeft = 2*plotLeft;        % Figure boundaries
streamRight = 2*plotRight;
streamTop = 2*plotTop;
streamBottom = 2*plotBottom;

sgtitle(sprintf('Flow Streamlines, $(x_0,y_0)=(%.2f, %.2f),b = %.2f,\\Gamma =%.2f, \\alpha=%.2f$', ...
    center(1),center(2),b,gam,a),'FontSize',titleSize,'Interpreter','latex');

subplot(1,2,1);

streamCs = -C:(U/2):C;
fpCirc = zeros(length(streamCs),1);
circXVals = [];     % streamline x-values
circYVals = [];     % streamline y-values

 for i = streamCs
     fpCirc = fimplicit(@(x,y) plotStreamf(x,y,i), [streamLeft streamRight streamBottom streamTop],'b-');
     circXVals = [circXVals; fpCirc.XData'];    % add another streamline's data to vector
     circYVals = [circYVals; fpCirc.YData'];    % " "
 end

hold on;
% CIRCLE STREAMLINES
plot(circ(1,:),circ(2,:),"r-","linewidth",1.3);
plot(circXVals, circYVals, "b.","MarkerSize",0.9*gridWidth)
xlim([streamLeft streamRight]);
ylim([streamBottom streamTop]);
xline(0,"linewidth",1); yline(0,"linewidth",1);
subtitle(sprintf("Pre-Transform"), "fontsize",subtitleSize,"interpreter","latex");

subplot(1,2,2)
hold on;

joukXVals = joukX(circXVals,circYVals);       % apply transform to streamline data
joukYVals = joukY(circXVals,circYVals);       % " "

% AIRFOIL STREAMLINES
plot(airf(1,:),airf(2,:),"r-","linewidth",1.3); 
plot(joukXVals, joukYVals, "b.","MarkerSize",0.9*gridWidth)
xlim([2*plotLeft 2*plotRight]);
ylim([2*plotBottom 2*plotTop]);
xline(0,"linewidth",1); yline(0,"linewidth",1);
subtitle(sprintf("Post-Transform"), "fontsize",subtitleSize,"interpreter","latex");

plt = plt + 1;
%-----------------------------------------------------
%% 3D blade - b functions (runnable as section)

plt = 4;
figure(plt);
clf
hold on;

plotHeight = 10;
plotWidth = 70;

sgtitle('$B(x_l)$ Examples','interpreter','latex');

subplot(3,1,1); % FIRST

r0 = 1.5;       % radius of blade at base
W = 4.5;          % measure of max blade width
K = 9;          % measure of initial widening sharpness
L = 30;         % Length of turbine blade
T = 0.1;          % slope of blade as it narrows

title(sprintf('i) $R_0 = %.2f, W=%.2f, K=%.2f, L=%.2f, T=%.2f$',r0,W,K,L,T),'Interpreter','latex','FontSize',6);

% points along turbine blade length
xL = @(L) linspace(0,L,80);

% func. varying 'b' parameter 
B = @(x, r0,W,K,L,T) r0 + ((W*K/L).*x)./(1+((K/L).*x).^2) - (T*K/L).*x;

plot(xL(L), B(xL(L), r0,W,K,L,T), '-r');
xlim([0 plotWidth]); ylim([0 plotHeight]);
xlabel('$x_l$','Interpreter','latex');
ylabel('$B(x_l)$','Interpreter','latex','Rotation',-pi/2);

subplot(3,1,2); % SECOND

r0 = 3;       % radius of blade at base
W = 7.5;          % measure of max blade width
K = 4.5;          % measure of initial widening sharpness
L = 60;         % Length of turbine blade
T = 0.75;          % slope of blade as it narrows

plot(xL(L), B(xL(L), r0,W,K,L,T), '-b');
xlim([0 plotWidth]); ylim([0 plotHeight]);
xlabel('$x_l$','Interpreter','latex');
ylabel('$B(x_l)$','Interpreter','latex','Rotation',-pi/2);


subplot(3,1,3); % THIRD

r0 = 2;       % radius of blade at base
W = 5;          % measure of max blade width
K = 3;          % measure of initial widening sharpness
L = 45;         % Length of turbine blade
T = 0.80;          % slope of blade as it narrows

xL = @(L) linspace(0,L,80);     % points along turbine blade
B = @(x, r0,W,K,L,T) r0 + ((W*K/L).*x)./(1+((K/L).*x).^2) - (T*K/L).*x;  % varying 'b' param. 

plot(xL(L), B(xL(L), r0,W,K,L,T), '-m');
xlim([0 plotWidth]); ylim([0 plotHeight]);
xlabel('$x_l$','Interpreter','latex');
ylabel('$B(x_l)$','Interpreter','latex','Rotation',-pi/2);

plt = plt + 1;
%% Functions
function grid_ = removePointsInside(grid,center,Radius)
    for i = 1:size(grid,1)
        for j = 1:size(grid,2)
            if norm([grid(i,j,1)-center(1);grid(i,j,2)-center(2)]) < Radius
                grid(i,j,1) = missing;
                grid(i,j,2) = missing;
            end
        end    
    end
    grid_ = grid;
end

