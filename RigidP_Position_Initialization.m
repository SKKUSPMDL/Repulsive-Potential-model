function [Pcenter, Seta, Lp, nump] = RigidP_Position_Initialization(nump, Hp)

% Plengthlist = ( [5, 10, 15, 20, 25, 30, 35, 40, 45, 60]  )*1e-6;
% 
% Prob = [3.125, 9.375, 16.56, 9.375, 17.1875, 22.5, 7.8125, 6.25, 1.5625, 6.25] * 1e-2;
% 
% Pnum = round( Prob*nump, 0);
% 
% ParticleLength = [];
% for i = 1:length(Prob)
%     ParticleLength = [ParticleLength, Plengthlist(i) * ones(1,Pnum(i) ) ];
% end
% % ParticleLength = [20*1e-6*ones(1,50), 30*1e-6*ones(1,50), 40*1e-6*ones(1,50)];
% 
% h = histfit(ParticleLength, 10);
% pd = fitdist(ParticleLength', 'Kernel');
% 
% 
% 
% 
% Plengthlist = [Plengthlist-2.5*1e-6, 65*1e-6, 70*1e-6, 75*1e-6];
% Prob = pdf(pd, Plengthlist ) * 5*1e-6;
% 
% Pnum = round( Prob*(nump+8), 0);
% nump = sum(Pnum);
% 
% ParticleLength = [];
% for i = 1:length(Prob)
%     ParticleLength = [ParticleLength, Plengthlist(i) * ones(1,Pnum(i) ) ];
% end



Plengthlist = 1.5 * ( [25, 35, 45, 55, 65, 75] )*1e-6;

Prob = [23.6, 23.6, 27, 10.9, 14.5, 0.4] * 1e-2;


Pnum = round( Prob*nump, 0);
nump = sum(Pnum);

ParticleLength = [];
for i = 1:length(Prob)
    ParticleLength = [ParticleLength, Plengthlist(i) * ones(1,Pnum(i) ) ];
end

% % y = cdf(pd, Plengthlist);

PsizeIdx = randperm( length(ParticleLength) );
% PsizeIdx = length(ParticleLength):-1:1;


% numpx = 5;
% numpy = 31;

numpx = 60;
numpy = 250;

% numpx = 26;
% numpy = 120;


xgrid = linspace( -250*1e-6, 0, numpx+3);
dxp = xgrid(2) - xgrid(1);
ygrid = linspace( -150*1e-6 , 250*1e-6 , numpy+3);
dyp = ygrid(2) - ygrid(1);

% check = 1;
Lp = ParticleLength( PsizeIdx );


Pcenter = zeros(nump,2);
Pcorner = zeros(nump,8);
% Seta = pi/2 * ones(nump,1);
Seta =  zeros(nump,1);

[X, Y] = meshgrid( xgrid(2:numpx+2), ygrid(2:numpy+2) );
ngrid = size(X,1) * size(X,2);

% nmix = randperm(ngrid);
% X = X(nmix);
% Y = Y(nmix);

for i = 1:nump

check = 1;
iter = 0;
l = Lp(i);
pcorn_sub = reshape( Pcorner(1:i-1,:), [i-1,2,4]);
pcent_sub = Pcenter(1:i-1,:);

% x1 = pcorn_sub(:,1,:);
% y1 = pcorn_sub(:,2,:);
% x2 = cat(3, pcorn_sub(:,1,2:end), pcorn_sub(:,1,1));
% y2 = cat(3, pcorn_sub(:,2,2:end), pcorn_sub(:,2,1));
% 
% q2q1 = zeros(i-1,2,4);
% 
% q2q1(:,:,1) = pcorn_sub(:,:,2) - pcorn_sub(:,:,1);
% q2q1(:,:,2) = pcorn_sub(:,:,3) - pcorn_sub(:,:,2);
% q2q1(:,:,3) = pcorn_sub(:,:,4) - pcorn_sub(:,:,3);
% q2q1(:,:,4) = pcorn_sub(:,:,1) - pcorn_sub(:,:,4);
% nq2q1 = vecnorm( q2q1, 2, 2);

x3 = repmat( pcorn_sub(:,1,:), [1,1,4]);
y3 = repmat( pcorn_sub(:,2,:), [1,1,4]);
x4 = repmat( cat(3, pcorn_sub(:,1,2:end), pcorn_sub(:,1,1)), [1,1,4]);
y4 = repmat( cat(3, pcorn_sub(:,2,2:end), pcorn_sub(:,2,1)), [1,1,4]);

% eta = Seta(i);
while check > 0 
    seta = pi/2 - 0.1*(rand-0.5);
    pcent = [ xgrid( randi([2, numpx+1]) ) , ygrid( randi( [2, numpy+1] ) ) ];
    % pcent = [ X(iter+1) , Y(iter+1) ];
    pv = 0.5*l*[cos(seta), sin(seta)];
    hv = 0.5*Hp*[cos(seta+pi/2), sin(seta+pi/2)];

    pcorn = [pcent+pv+hv, pcent+pv-hv, pcent-pv-hv, pcent-pv+hv];
    
    x1 = repmat( reshape( repmat( pcorn(1,[1,3,5,7]), [4,1] ), [1,1,16]), [i-1,1,1])  ;
    y1 = repmat( reshape( repmat( pcorn(1,[2,4,6,8]), [4,1] ), [1,1,16]), [i-1,1,1])  ;

    x2 = repmat( reshape( repmat( pcorn(1,[3,5,7,1]), [4,1] ), [1,1,16]), [i-1,1,1])  ;
    y2 = repmat( reshape( repmat( pcorn(1,[4,6,8,2]), [4,1] ), [1,1,16]), [i-1,1,1])  ;
    
    %------------------------------------------------------------------------------------------------------------------------
    % Line segments intersect parameters
    u = ((x1-x3).*(y1-y2) - (y1-y3).*(x1-x2)) ./ ((x1-x2).*(y3-y4)-(y1-y2).*(x3-x4));
    t = ((x1-x3).*(y3-y4) - (y1-y3).*(x3-x4)) ./ ((x1-x2).*(y3-y4)-(y1-y2).*(x3-x4));
    %------------------------------------------------------------------------------------------------------------------------

    check = sum( (u >= -0.1 & u <= 1.1) & (t >= -0.1 & t <= 1.1), 'all') + ( max(pcorn([1,3,5,7]) ) >= 0 | min(pcorn([1,3,5,7]) )  <= -250*1e-6 ) | ( max(pcorn([2,4,6,8]) ) >= 250*1e-6 | min(pcorn([2,4,6,8]) ) <= -150*1e-6) ;
    


    iter = iter+1;

    if iter> ngrid-1
        1
        break
    end
end
Pcenter(i,:) = pcent;
Pcorner(i,:) = pcorn;
Seta(i,1) = seta;

end


end