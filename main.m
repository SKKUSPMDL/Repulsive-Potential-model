clc; clear all; close all

numIter = 1;
errCount = 0;
sucCount = 0;
tryNum = 0;

errIdx = [];



c = 4;

path = "./Case"+num2str(c) ;

if c == 1
    mt = 650; %%
    % siterNum = 89000;
elseif c == 2
    mt = 300; %%
    % siterNum = 89000;
elseif c == 3
    mt = 200; %%
    % siterNum = 89000;
elseif c == 4
    mt = 40; %?
end

% siterNum = 1000;


for i = 1:numIter
    
    while true
        try
            

            [ uelem, velem, pelem, dx, dy, xlimit, ylimit,  L, H, Lp, Hp, density, gamma, dt, Node, Element, U, Wall, viscosity] = RigidP_Fluid_Case(c);
            dt = dt * mt;

            aidx = min(xlimit, [], 2) >= 0;
        
            Umax = max(U(aidx));
            siterNum = 81e-6 / (0.5*Umax*dt);

            nump = 313;
            [Pcenter, Seta, Lp, nump] = RigidP_Position_Initialization(nump, Hp);
            d_L = 0.1e-8;
            d_col = 1e-6;f_col = 1.0;
            d_wcol = 2.5e-6;f_wcol = 1.0;
            
            Lp = Lp'; 
            
            Pvect = 0.5*Lp.*[cos(Seta), sin(Seta)];
            Hvect = 0.5*Hp.*[cos(pi/2 + Seta), sin(pi/2 + Seta)];
            
            rho_BN = 2100;
            density = 1410;
            
            Mass = rho_BN * Lp * Hp;
            I = (Lp.^2 + Hp.^2) .* Mass / 12;
            V = zeros(nump,2);
            omega = zeros(nump,1);
            
            Pcorner = [Pcenter + Pvect + Hvect, Pcenter + Pvect - Hvect, Pcenter - Pvect - Hvect, Pcenter - Pvect + Hvect];
            
            Snum = 5*ones(1,4);
            Pxgrid = zeros(Snum(1), Snum(1), nump);
            Pygrid = zeros(Snum(1), Snum(1), nump);
            dA = zeros(1,1,nump);
            for ig = 1:1:nump
                [X,Y] = meshgrid( linspace(-Lp(ig)/2,Lp(ig)/2,Snum(1)), linspace(-Hp/2,Hp/2,Snum(1)) );
            
                delx = X(1,2) - X(1,1);
                dely = Y(2,1) - Y(1,1);
            
                dA(ig) = delx*dely;
                Pxgrid(:,:,ig) = X;
                Pygrid(:,:,ig) = Y;
            
            end
            
            Slength = [Hp*ones(nump,1), Lp, Hp*ones(nump,1), Lp];


            Pcorner = gpuArray(Pcorner);Pcenter = gpuArray(Pcenter);Pvect = gpuArray(Pvect);Hvect = gpuArray(Hvect);Seta = gpuArray(Seta);Node = gpuArray(Node);Element = gpuArray(Element);
            uelem = gpuArray(uelem);velem = gpuArray(velem);xlimit = gpuArray(xlimit);ylimit = gpuArray(ylimit);Lp = gpuArray(Lp);Mass = gpuArray(Mass);I = gpuArray(I);V = gpuArray(V);omega = gpuArray(omega);Slength = gpuArray(Slength);dx = gpuArray(dx);dy = gpuArray(dy);
            Pxgrid = gpuArray(Pxgrid);Pygrid = gpuArray(Pygrid);dA = gpuArray(dA);Wall = gpuArray(Wall);Snum = gpuArray(Snum);
            
            f1 = figure;
            PlotFieldonMesh(Node(:,2:3), Element, U )
            hold on
            xp =gather( Pcorner(:,[1,3,5,7,1]) );
            yp = gather( Pcorner(:,[2,4,6,8,2]) );
            for k = 1:nump
                fplot(k) = plot( xp( k, : ), yp( k, : ), 'LineWidth', 1);
            end
            for k = 1:nump
                set(fplot(k), 'XDataSource', sprintf('xp(%d,:)', k ) ,...
                              'YDataSource', sprintf('yp(%d,:)', k ))
            end

            v = VideoWriter( path + "/" + sprintf("Case%d_T%d_v[dL%.d_dcol%.d_fcol%.1f_dwcol%.d_fwcol%.1f]",c,tryNum,d_L,d_col,f_col,d_wcol,f_wcol) + ".avi");
            open(v);

            uptidx = gpuArray( upperTriangularIndexing(nump) );
            lotidx = gpuArray( lowerTriangularIndexing(nump) );

            
            for is = 0:siterNum
            
            V_old = V;
            omega_old = omega;
            
            [F_ibm, rF_ibm, Pu, Pru] = HydronamicForce( d_L, Pcorner, Pcenter, Seta, V_old, omega_old, Snum, Slength, dA, xlimit, ylimit, uelem, velem, dx, dy, Node, Element, Pxgrid, Pygrid);
            
            if is == 0
                Pru_n_1 = Pru;
                Pru_n = Pru;

                Pu_n_1 = Pu;
                Pu_n = Pu;
            else
                Pru_n_1 = Pru_n;
                Pru_n = Pru;
            
                Pu_n_1 = Pu_n;
                Pu_n = Pu;
            end
            
            V = V_old + (-density.*F_ibm + density.*(Pu_n - Pu_n_1) ) ./ Mass;
            omega = omega_old + (-density.*rF_ibm + density.*(Pru_n - Pru_n_1) ) ./ I;
            overIdx =find( max(Pcorner(:,[1,3,5,7]) , [], 2) >= 0.99 * 100e-6 );
            V(overIdx, :) = 0;
            omega(overIdx,:) = 0;
            
            [VFC, omegaFC] = Collision2( d_col, f_col, Pcorner, nump, Pcenter, Mass, I, V, omega, uptidx, lotidx, overIdx);
            

            
            [VFW, omegaFW] = WallCollision( d_wcol, f_wcol, Pcorner, nump, Pcenter, Mass, I, V, omega, Wall );
            
            V = V + (VFC-V) + (VFW-V);
            omega = omega + (omegaFC-omega) + (omegaFW-omega);
            
            V(overIdx, :) = 0;
            omega(overIdx,:) = 0;
            
            Pcenter = Pcenter + V*dt;
            Seta = Seta + omega*dt;
            
            Pvect = 0.5*Lp.*[cos(Seta), sin(Seta)];
            Hvect = 0.5*Hp.*[cos(pi/2 + Seta), sin(pi/2 + Seta)];
            Pcorner = [Pcenter + Pvect + Hvect, Pcenter + Pvect - Hvect, Pcenter - Pvect - Hvect, Pcenter - Pvect + Hvect];
            
            if mod(is, 200) == 0
                is
                xp = gather( Pcorner(:,[1,3,5,7,1]) );
                yp = gather( Pcorner(:,[2,4,6 ,8,2]) );
                refreshdata(fplot)
                drawnow
                frame = getframe(gcf);
                writeVideo(v, frame);
            end
            
            end
            
            save( path + "/Result"+"/" + num2str(i) +".mat" , "Pcenter", "Pcorner", "Pvect", "Hvect", "Seta" )

            exportgraphics(f1, path + "/"+"Case"+num2str(c)+"_"+"T" + num2str(tryNum) + ".png", "Resolution", 300)
            close(v)
            close all
            break

        catch
            xp = gather( Pcorner(:,[1,3,5,7,1]) );
            yp = gather( Pcorner(:,[2,4,6,8,2]) );
            refreshdata(fplot)
            drawnow
            frame = getframe(gcf);
            writeVideo(v, frame);
            close(v)

            save( path + "/Error"+"/T" + num2str(tryNum) +".mat" , "Pcenter", "Pcorner", "Pvect", "Hvect", "Seta" )
            exportgraphics(f1, path + "/"+"Case"+num2str(c)+"_"+"T" + num2str(tryNum) + ".png", "Resolution", 300)

            errIdx = [errIdx,1];
            errCount = errCount + 1;
            tryNum = tryNum + 1;

            close all
            reset(gpuDevice)
        end
        
        % clearvars -except numIter errCount sucCount tryNum errIdx c i path
    end
    
    tryNum = tryNum + 1;
    errIdx = [errIdx,0];
    reset(gpuDevice)
end











%%
function [F_ibm, rF_ibm, Pu, Pru_sum] = HydronamicForce( d_l, Pcorner, Pcenter, seta, V, omega, Snum, Slength, dA, xlimit, ylimit, uelem, velem, dx, dy, Node, Element, Pgridx, Pgridy)
    % d_l = 0.1e-8;
    % d_l = 1.3889e-6;
    Np = size(Pcorner,1);
% 
    % Node = Node(:,2:3)';
    % Element = Element';
    Pcorner = reshape( Pcorner, [Np,2,4]);
    
    line2 = cat(3, Pcorner(:,:,2:end), Pcorner(:,:,1) );
    line1 = Pcorner;
    ds = Slength ./ (Snum-1);
    ds = reshape( ds, [Np,1,4]);

    xmat = zeros( Np, Snum(1), 4, 'gpuArray');
    ymat = zeros( Np, Snum(1), 4, 'gpuArray');

    xmat(:,end,:) = line2(:,1,:);
    xmat(:,1,:) = line1(:,1,:);

    ymat(:,end,:) = line2(:,2,:);
    ymat(:,1,:) = line1(:,2,:);


    xmat(:,2:end-1,:) = line1(:,1,:) + ( ( line2(:,1,:) - line1(:,1,:) ) / ( Snum(1) - 1 )  ) .* (1:Snum(1)-2) ;
    ymat(:,2:end-1,:) = line1(:,2,:) + ( ( line2(:,2,:) - line1(:,2,:) ) / ( Snum(1) - 1 )  ) .* (1:Snum(1)-2) ;

    xvec = reshape( xmat, [], 1);
    yvec = reshape( ymat, [], 1);

    
    u = U_IBG([xvec, yvec], xlimit, ylimit, uelem, velem, dx, dy, Node, Element);
    
    ui = reshape( u(:,1), [Np, Snum(1), 4]);
    vi = reshape( u(:,2), [Np, Snum(1), 4]);

    rx = xmat - Pcenter(:,1);
    ry = ymat - Pcenter(:,2);

    ux_gamma = V(:,1) + (-omega.*ry);
    uy_gamma = V(:,2) + (omega.*rx);

    fx_ibm = ux_gamma - ui;
    fy_ibm = uy_gamma - vi;

    F_ibm = zeros(Np,2, 'gpuArray');
    rF_ibm = zeros(Np,1, 'gpuArray');


    F_ibm(:,1) = F_ibm(:,1) + sum( sum( 0.5 * (fx_ibm(:,2:end,:) + fx_ibm(:,1:end-1,:)) .* ds * d_l, 2), 3);
    F_ibm(:,2) = F_ibm(:,2) + sum( sum( 0.5 * (fy_ibm(:,2:end,:) + fy_ibm(:,1:end-1,:)) .* ds * d_l, 2), 3);

    rf_ibm = rx.*fy_ibm - ry.*fx_ibm;

    rF_ibm(:,1) = rF_ibm(:,1) + sum( sum( 0.5 * (rf_ibm(:,2:end,:) + rf_ibm(:,1:end-1,:)) .* ds * d_l, 2), 3);




    
    cos_seta = reshape( cos(seta), [1,1,Np]);
    sin_seta = reshape( sin(seta), [1,1,Np]);

    Pgridx_ = Pgridx .* cos_seta - Pgridy .* sin_seta;
    Pgridy_ = Pgridx .* sin_seta + Pgridy .* cos_seta;

    Pgridx_ = Pgridx_ + reshape(Pcenter(:,1),[1,1,Np]);
    Pgridy_ = Pgridy_ + reshape(Pcenter(:,2),[1,1,Np]);

    tx = reshape(Pgridx_, [], 1);
    ty = reshape(Pgridy_, [], 1);

    ut = U_IBG([tx, ty], xlimit, ylimit, uelem, velem, dx, dy, Node, Element);
    
    Puxgrid = reshape( ut(:,1), 5,5,Np);
    Puygrid = reshape( ut(:,2), 5,5,Np);

    Pux = sum( ( ( Puxgrid(1:end-1, 1:end-1,:) + Puxgrid(2:end, 1:end-1,:) + Puxgrid(1:end-1, 2:end,:) + Puxgrid(2:end, 2:end,:) ) / 4 ) .* dA, [1, 2]);
    Puy = sum( ( ( Puygrid(1:end-1, 1:end-1,:) + Puygrid(2:end, 1:end-1,:) + Puygrid(1:end-1, 2:end,:) + Puygrid(2:end, 2:end,:) ) / 4 ) .* dA, [1, 2]);

    Pu = [ reshape(Pux,Np,1), reshape(Puy,Np,1)];

    Prx = Pgridx_ - reshape( Pcenter(:,1), [1,1,Np]);
    Pry = Pgridy_ - reshape( Pcenter(:,2), [1,1,Np]);
    
    Pru = Prx.*Puygrid - Pry.*Puxgrid;
    Pru_sum = sum( ( ( Pru(1:end-1, 1:end-1,:) + Pru(2:end, 1:end-1,:) + Pru(1:end-1, 2:end,:) + Pru(2:end, 2:end,:) ) / 4 ) .* dA, [1, 2]);

    Pru_sum = reshape( Pru_sum,[Np,1]);
    

    
    

end

function [VF, omegaF] = WallCollision( d_wcol, f_wcol, Pcorner, nump, Pcenter, Mass, I, V, omega, Wall )
    nw = size(Wall,1);
    Corner1 = Wall(4,1:2);
    Corner2 = Wall(6,1:2);

    

    wQ2Q1 = Wall(:,3:4) - Wall(:,1:2);

    nwQ2Q1 = vecnorm( wQ2Q1, 2, 2) ;

    Pcorner = reshape( Pcorner, [nump,2,4] );
    
    Q2Q1 = cat(3, Pcorner(:,:,2:end), Pcorner(:,:,1) ) - Pcorner;
    nQ2Q1 = vecnorm( Q2Q1, 2, 2);

    PQc1x = Corner1(:,1) - ( Pcorner(:,1,:) );
    PQc1y = Corner1(:,2) - ( Pcorner(:,2,:) );

    PQc2x = Corner2(:,1) - ( Pcorner(:,1,:) );
    PQc2y = Corner2(:,2) - ( Pcorner(:,2,:) );

    wPQx = pagetranspose( Pcorner(:,1,:) ) - Wall(:,1);
    wPQy = pagetranspose( Pcorner(:,2,:) ) - Wall(:,2);
    % 


    dnorm = ( 0 < (wQ2Q1(:,1).*wPQx + wQ2Q1(:,2).*wPQy) ./ (nwQ2Q1(:,1).^2) & (wQ2Q1(:,1).*wPQx + wQ2Q1(:,2).*wPQy) ./ (nwQ2Q1(:,1).^2)  < 1    )  ;

    dmat = abs( wQ2Q1(:,1) .* wPQy - wQ2Q1(:,2) .* wPQx ) ./ nwQ2Q1(:,1) .*  dnorm ;

    dmat(dmat==0) = 1e+4;

    [d, didx] = min(dmat, [], 3);
    [d2, d2idx] = mink(dmat, 2, 3);
    d2 = d2(:,:,end);
    d2idx = d2idx(:,:,end);

    % Paridx = abs(d-d2)< 1e-7 & d<1e+4 & d2<1e+4;

    dcorner = rem(didx,4);
    dcorner(dcorner==0) = 4;

    dcorner2 = rem(d2idx,4);
    dcorner2(dcorner2==0) = 4;
    
    x0 = zeros(nw, nump, 'gpuArray'); 
    y0 = zeros(nw, nump, 'gpuArray');

    % rx = zeros(nump, nump);
    % ry = zeros(nump, nump);


    [~,dcC1] = find( dcorner == 1);
    x0( dcorner == 1) = Pcorner(dcC1, 1,1); 
    y0( dcorner == 1) = Pcorner(dcC1, 2,1); 

    x1 = repmat(Wall(:,1), [1, nump]);
    x2 = repmat(Wall(:,3), [1, nump]);
    y1 = repmat(Wall(:,2), [1, nump]);
    y2 = repmat(Wall(:,4), [1, nump]);

    [~,dcC2] = find( dcorner == 2);
    x0( dcorner == 2) = Pcorner(dcC2, 1,2);
    y0( dcorner == 2) = Pcorner(dcC2, 2,2);

    [~,dcC3] = find( dcorner == 3);
    x0( dcorner == 3) = Pcorner(dcC3, 1,3);
    y0( dcorner == 3) = Pcorner(dcC3, 2,3);

    [~,dcC4] = find( dcorner == 4);
    x0( dcorner == 4 ) = Pcorner(dcC4, 1,4);
    y0( dcorner == 4 ) = Pcorner(dcC4, 2,4);


    x20 = zeros(nw, nump, 'gpuArray'); 
    y20 = zeros(nw, nump, 'gpuArray');


    [~,dcC1] = find( dcorner2 == 1);
    x20( dcorner2 == 1) = Pcorner(dcC1, 1,1); 
    y20( dcorner2 == 1) = Pcorner(dcC1, 2,1); 

    [~,dcC2] = find( dcorner2 == 2);
    x20( dcorner2 == 2) = Pcorner(dcC2, 1,2);
    y20( dcorner2 == 2) = Pcorner(dcC2, 2,2);

    [~,dcC3] = find( dcorner2 == 3);
    x20( dcorner2 == 3) = Pcorner(dcC3, 1,3);
    y20( dcorner2 == 3) = Pcorner(dcC3, 2,3);

    [~,dcC4] = find( dcorner2 == 4);
    x20( dcorner2 == 4 ) = Pcorner(dcC4, 1,4);
    y20( dcorner2 == 4 ) = Pcorner(dcC4, 2,4);


    % x0(Paridx) = (x0(Paridx) + x20(Paridx) ) / 2;
    % y0(Paridx) = (y0(Paridx) + y20(Paridx) ) / 2;

    tx = (x1 - x2) ./ sqrt( (x1 - x2).^2 + (y1 - y2).^2   );
    ty = (y1 - y2) ./ sqrt( (x1 - x2).^2 + (y1 - y2).^2   );

    nx = (x1 - x0) - ( (x1-x0).*tx + (y1-y0).*ty ) .* tx;
    ny = (y1 - y0) - ( (x1-x0).*tx + (y1-y0).*ty ) .* ty;

    d(d<1e+4) = (nx(d<1e+4).^2 + ny(d<1e+4).^2).^0.5;




    dc1mat = abs( Q2Q1(:,1,:) .* PQc1y - Q2Q1(:,2,:) .* PQc1x ) ./ nQ2Q1 .* ( 0 < ( Q2Q1(:,1,:) .* PQc1x + Q2Q1(:,2,:) .* PQc1y ) ./ (nQ2Q1.^2) & ( Q2Q1(:,1,:) .* PQc1x + Q2Q1(:,2,:) .* PQc1y ) ./ (nQ2Q1.^2) < 1 );

    dc2mat = abs( Q2Q1(:,1,:) .* PQc2y - Q2Q1(:,2,:) .* PQc2x ) ./ nQ2Q1 .* ( 0 < ( Q2Q1(:,1,:) .* PQc2x + Q2Q1(:,2,:) .* PQc2y ) ./ (nQ2Q1.^2) & ( Q2Q1(:,1,:) .* PQc2x + Q2Q1(:,2,:) .* PQc2y ) ./ (nQ2Q1.^2) < 1 );
    
    dc1mat(dc1mat==0) = 1e+4;
    dc2mat(dc2mat==0) = 1e+4;

    [dc1, dc1idx] = min(dc1mat, [], 3);
    [dc2, dc2idx] = min(dc2mat, [], 3);
    
    dc1idx_sub = dc1idx + 1;
    dc1idx_sub( dc1idx_sub == 5) = 1;

    dc2idx_sub = dc2idx + 1;
    dc2idx_sub( dc2idx_sub == 5) = 1;

    xc11 = Pcorner( sub2ind( size(Pcorner), (1:nump)', ones(nump,1), dc1idx ) );
    xc12 = Pcorner( sub2ind( size(Pcorner), (1:nump)', ones(nump,1), dc1idx_sub ) );
    yc11 = Pcorner( sub2ind( size(Pcorner), (1:nump)', 2*ones(nump,1), dc1idx ) );
    yc12 = Pcorner( sub2ind( size(Pcorner), (1:nump)', 2*ones(nump,1), dc1idx_sub ) );
    xc10 = Corner1(1) * ones( 1, nump, 'gpuArray' );
    yc10 = Corner1(2) * ones( 1, nump, 'gpuArray' );


    xc21 = Pcorner( sub2ind( size(Pcorner), (1:nump)', ones(nump,1), dc2idx ) );
    xc22 = Pcorner( sub2ind( size(Pcorner), (1:nump)', ones(nump,1), dc2idx_sub ) );
    yc21 = Pcorner( sub2ind( size(Pcorner), (1:nump)', 2*ones(nump,1), dc2idx ) );
    yc22 = Pcorner( sub2ind( size(Pcorner), (1:nump)', 2*ones(nump,1), dc2idx_sub ) );
    xc20 = Corner2(1) * ones(1,nump, 'gpuArray');
    yc20 = Corner2(2) * ones(1,nump, 'gpuArray');

    x0_sub = zeros(nw-2,nump, 'gpuArray');
    y0_sub = zeros(nw-2,nump, 'gpuArray');
    x1_sub = zeros(nw-2,nump, 'gpuArray');
    y1_sub = zeros(nw-2,nump, 'gpuArray');
    x2_sub = zeros(nw-2,nump, 'gpuArray');
    y2_sub = zeros(nw-2,nump, 'gpuArray');
    d_sub = zeros(nw-2,nump, 'gpuArray');

    x0_sub(1:3,:) = x0(1:3,:);
    y0_sub(1:3,:) = y0(1:3,:);

    x1_sub(1:3,:) = x1(1:3,:);
    y1_sub(1:3,:) = y1(1:3,:);

    x2_sub(1:3,:) = x2(1:3,:);
    y2_sub(1:3,:) = y2(1:3,:);
    
    d_sub(1:3,:) = d(1:3,:);
    % 
    % x0_sub(end,:) = x0(end,:);
    % y0_sub(end,:) = y0(end,:);
    % 
    % x1_sub(end,:) = x1(end,:);
    % y1_sub(end,:) = y1(end,:);
    % 
    % x2_sub(end,:) = x2(end,:);
    % y2_sub(end,:) = y2(end,:);

    d_sub(end,:) = d(end,:);

    [dcat1 , dcat1idx]= min( cat(3, d(4,:), d(5,:), dc1'), [], 3);
    
    x0sub = cat(3, x0(4,:), x0(5,:), xc10);
    x1sub = cat(3, x1(4,:), x1(5,:), xc11');
    x2sub = cat(3, x2(4,:), x2(5,:), xc12');
    y0sub = cat(3, y0(4,:), y0(5,:), yc10);
    y1sub = cat(3, y1(4,:), y1(5,:), yc11');
    y2sub = cat(3, y2(4,:), y2(5,:), yc12');


    c1idx = sub2ind( size(x0sub), ones(nump, 1), (1:nump)', dcat1idx');

    d_sub(4,:) = dcat1;
    x0_sub(4,:)  = x0sub(c1idx);
    y0_sub(4,:) = y0sub(c1idx);

    x1_sub(4,:)  = x1sub(c1idx);
    y1_sub(4,:) = y1sub(c1idx);
    
    x2_sub(4,:)  = x2sub(c1idx);
    y2_sub(4,:) = y2sub(c1idx);
     
    

    [dcat2 , dcat2idx]= min( cat(3, d(6,:), d(7,:), dc2'), [], 3);

    x0sub = cat(3, x0(6,:), x0(7,:), xc20);
    x1sub = cat(3, x1(6,:), x1(7,:), xc21');
    x2sub = cat(3, x2(6,:), x2(7,:), xc22');
    y0sub = cat(3, y0(6,:), y0(7,:), yc20);
    y1sub = cat(3, y1(6,:), y1(7,:), yc21');
    y2sub = cat(3, y2(6,:), y2(7,:), yc22');

    c2idx = sub2ind( size(x0sub), ones(nump, 1), (1:nump)', dcat2idx');

    d_sub(5,:) = dcat2;
    x0_sub(5,:)  = x0sub(c2idx);
    y0_sub(5,:) = y0sub(c2idx);

    x1_sub(5,:)  = x1sub(c2idx);
    y1_sub(5,:) = y1sub(c2idx);
    
    x2_sub(5,:)  = x2sub(c2idx);
    y2_sub(5,:) = y2sub(c2idx);

    tx = (x1_sub - x2_sub) ./ sqrt( (x1_sub - x2_sub).^2 + (y1_sub - y2_sub).^2   );
    ty = (y1_sub - y2_sub) ./ sqrt( (x1_sub - x2_sub).^2 + (y1_sub - y2_sub).^2   );

    nx = (x1_sub - x0_sub) - ( (x1_sub-x0_sub).*tx + (y1_sub-y0_sub).*ty ) .* tx;
    ny = (y1_sub - y0_sub) - ( (x1_sub-x0_sub).*tx + (y1_sub-y0_sub).*ty ) .* ty;

    qx = x0_sub + nx;
    qy = y0_sub + ny;
     
    d_sub(d_sub<1e+4) = (nx(d_sub<1e+4).^2 + ny(d_sub<1e+4).^2).^0.5;

   VF = V;
   omegaF = omega;


    rx = x0_sub - Pcenter(:,1)';
    ry = y0_sub - Pcenter(:,2)';

    dcat1corn = dcat1idx == 3;
    rx(4,dcat1corn) = qx(4,dcat1corn) - Pcenter(dcat1corn,1)';
    ry(4,dcat1corn) = qy(4,dcat1corn) - Pcenter(dcat1corn,2)';
    
    dcat2corn = dcat2idx == 3;
    rx(5,dcat2corn) = qx(5,dcat2corn) - Pcenter(dcat2corn,1)';
    ry(5,dcat2corn) = qy(5,dcat2corn) - Pcenter(dcat2corn,2)';
    

    % 
    Chaix = V(:,1)' - omega' .* ry;
    Chaiy = V(:,2)' + omega' .* rx;


    % paral = abs( d_ab - d_ba ) < 0.1*1e-6;

    % d(d<1e+4) = sqrt( nx(d<1e+4).^2 + ny(d<1e+4).^2 );

    rm_idx = d_sub < d_wcol;
    % [rw_idx, rp_idx] = find(rm_idx);

    b = - (f_wcol) .*( (Chaix .* nx + Chaiy .* ny)   );
    
    r_n = rx.*ny - ry.*nx;
    

    alpha = (nx.^2 + ny.^2) ./ Mass' + ( -r_n.*ry.*nx + r_n.*rx.*ny) ./ I' ; 

    f = b ./ alpha;
    
    % nx(4,dc1P) = -nx(4,dc1P);
    % ny(4,dc1P) = -ny(4,dc1P);
    % nx(5,dc1P) = 0;
    % ny(5,dc1P) = 0;
    % nx(6,dc2P) = -nx(6,dc2P);
    % ny(6,dc2P) = -ny(6,dc2P);
    % nx(7,dc2P) = 0;
    % ny(7,dc2P) = 0;

    % f(4,dc1P) = -f(4,dc1P);
    % f(5,dc1P) = 0;
    % % f(6,dc2P) = -f(6,dc2P);
    % f(7,dc2P) = 0;



    % checkdir = (rx.*nx + ry.*ny) >= 0 ;

    % nx(checkdir) = -nx(checkdir);
    % ny(checkdir) = -ny(checkdir);
    % f = abs(f);

    
    VF(:,1) = VF(:,1) + sum(f.*nx.*rm_idx, 1)' ./ Mass;

    VF(:,2) = VF(:,2) + sum(f.*ny.*rm_idx, 1)' ./ Mass;

    omegaF = omegaF + sum( f.*(rx.*ny - ry.*nx).*rm_idx,1 )' ./ I;




end

function [VF, omegaF] = Collision2(d_col, f_col,Pcorner, nump, Pcenter, Mass, I, V, omega, uptidx, lotidx, overIdx )


    Pcorn_1234 = reshape( repmat( Pcorner , [1,4]), [nump,2,16]); 

    Pcorner = reshape( Pcorner, [],2,4);

    Pcorn_1111 = reshape( repmat( Pcorner, [1,4,1]) , [nump,2,16]);



    Q2Q1 = cat(3, Pcorner(:,:,2:end), Pcorner(:,:,1) ) - Pcorner;

    Q2Q1_sub = reshape( repmat(Q2Q1, [1,4,1]), [nump,2,16]) ;



    nQ2Q1 = vecnorm( Q2Q1, 2, 2) ;
    nQ2Q1_sub = reshape( repmat(nQ2Q1, [1,4,1]), [nump,1,16]) ;

    PQx = pagetranspose( Pcorn_1234(:,1,:) ) - Pcorn_1111(:,1,:);
    PQy = pagetranspose( Pcorn_1234(:,2,:) ) - Pcorn_1111(:,2,:);


    dnorm = 0 < ( Q2Q1_sub(:,1,:) .* PQx + Q2Q1_sub(:,2,:) .* PQy ) ./ (nQ2Q1_sub.^2) & ( Q2Q1_sub(:,1,:) .* PQx + Q2Q1_sub(:,2,:) .* PQy ) ./ (nQ2Q1_sub.^2) < 1;


    dmat = abs( Q2Q1_sub(:,1,:) .* PQy - Q2Q1_sub(:,2,:) .* PQx ) ./ nQ2Q1_sub .* ( dnorm ) + sqrt(PQx.^2 + PQy.^2) .* (~dnorm)  ;
    % dmat = abs( Q2Q1_sub(:,1,:) .* PQy - Q2Q1_sub(:,2,:) .* PQx ) ./ nQ2Q1_sub .* ( dnorm ) + 1e+4 .* (~dnorm)  ;

    [d, didx] = min(dmat, [], 3);

    subidx = gpuArray( repmat( 1:nump,[nump,1]) );

    D2 = reshape( dnorm( sub2ind( size(dmat), reshape( subidx',[],1), reshape( subidx,[],1), reshape( didx,[],1) ) ), [nump, nump]);

    dcorner = rem(didx,4);
    dcorner(dcorner==0) = 4;
    dsurf = fix( (didx-1) /4)+1;

    x0 = zeros(nump, nump, 'gpuArray');
    y0 = zeros(nump, nump, 'gpuArray');
    x1 = zeros(nump, nump, 'gpuArray');
    x2 = zeros(nump, nump, 'gpuArray');
    y1 = zeros(nump, nump, 'gpuArray');
    y2 = zeros(nump, nump, 'gpuArray');


    [~,dcC1] = find( dcorner == 1);
    x0( dcorner == 1) = Pcorner(dcC1, 1,1); 
    y0( dcorner == 1) = Pcorner(dcC1, 2,1); 

    [dcS1,~] = find( dsurf == 1);
    x1( dsurf == 1 ) = Pcorner(dcS1, 1,1);
    x2( dsurf == 1 ) = Pcorner(dcS1, 1,2);
    y1( dsurf == 1 ) = Pcorner(dcS1, 2,1);
    y2( dsurf == 1 ) = Pcorner(dcS1, 2,2);

    [~,dcC2] = find( dcorner == 2);
    x0( dcorner == 2) = Pcorner(dcC2, 1,2);
    y0( dcorner == 2) = Pcorner(dcC2, 2,2);

    [dcS2,~] = find( dsurf == 2);
    x1( dsurf == 2 ) = Pcorner(dcS2, 1,2);
    x2( dsurf == 2 ) = Pcorner(dcS2, 1,3);
    y1( dsurf == 2 ) = Pcorner(dcS2, 2,2);
    y2( dsurf == 2 ) = Pcorner(dcS2, 2,3);


    [~,dcC3] = find( dcorner == 3);
    x0( dcorner == 3) = Pcorner(dcC3, 1,3);
    y0( dcorner == 3) = Pcorner(dcC3, 2,3);

    [dcS3,~] = find( dsurf == 3);
    x1( dsurf == 3 ) = Pcorner(dcS3, 1,3);
    x2( dsurf == 3 ) = Pcorner(dcS3, 1,4);
    y1( dsurf == 3 ) = Pcorner(dcS3, 2,3);
    y2( dsurf == 3 ) = Pcorner(dcS3, 2,4);



    [~,dcC4] = find( dcorner == 4);
    x0( dcorner == 4 ) = Pcorner(dcC4, 1,4);
    y0( dcorner == 4 ) = Pcorner(dcC4, 2,4);

    [dcS4, ~] = find( dsurf == 4);
    x1( dsurf == 4 ) = Pcorner(dcS4, 1,4);
    x2( dsurf == 4 ) = Pcorner(dcS4, 1,1);
    y1( dsurf == 4 ) = Pcorner(dcS4, 2,4);
    y2( dsurf == 4 ) = Pcorner(dcS4, 2,1);

    tx = (x1 - x2) ./ sqrt( (x1 - x2).^2 + (y1 - y2).^2   );
    ty = (y1 - y2) ./ sqrt( (x1 - x2).^2 + (y1 - y2).^2   );

    nx = (x1 - x0) - ( (x1-x0).*tx + (y1-y0).*ty ) .* tx;
    ny = (y1 - y0) - ( (x1-x0).*tx + (y1-y0).*ty ) .* ty;

    nx(~D2) = x1(~D2) - x0(~D2);
    ny(~D2) = y1(~D2) - y0(~D2);

    d = ( nx.^2 + ny.^2 ).^0.5;

    qx = x0 + nx;
    qy = y0 + ny;

    ijmat = gpuArray( repmat( 1:nump, [nump,1]) );
    

   VF = V;
   omegaF = omega;

   d_ab = d(lotidx);
   d_ba = d(uptidx);

   D2_sub = [D2(lotidx), D2(uptidx)];


   ij = [ijmat(lotidx), ijmat(uptidx)];

   rx1 = x0 - Pcenter(:,1)';
   rx2 = qx - Pcenter(:,1);
   ry1 = y0 - Pcenter(:,2)';
   ry2 = qy - Pcenter(:,2);


   rx = [ rx1(lotidx), rx2(uptidx)];
   ry = [ ry1(lotidx), ry2(uptidx)];

   rx_sub = [rx2(lotidx), rx1(uptidx)];
   ry_sub = [ry2(lotidx), ry1(uptidx)];

   vxmat = repmat( V(:,1), [1,nump]);
   vymat = repmat( V(:,2), [1,nump]);

   omegamat = repmat( omega, [1, nump]);

   Chaix = vxmat(uptidx) - omegamat(uptidx) .* ry;
   Chaiy = vymat(uptidx) + omegamat(uptidx) .* rx;

   Chaix_sub = vxmat(lotidx) - omegamat(lotidx) .* ry_sub;
   Chaiy_sub = vymat(lotidx) + omegamat(lotidx) .* rx_sub;

   qx_sub = [x0(lotidx), qx(uptidx)];
   qy_sub = [y0(lotidx), qy(uptidx)];

   qx_sub2 = [qx(lotidx), x0(uptidx)];
   qy_sub2 = [qy(lotidx), y0(uptidx)];

   nx_sub = [nx(lotidx), nx(uptidx)];
   ny_sub = [ny(lotidx), ny(uptidx)];

   [d_min, d_midx] = min( [d_ab, d_ba], [], 2);    

   minIdx = sub2ind( size(rx), (1:size(rx,1))', d_midx);

   D2_min = D2_sub(minIdx);

   nx_min = nx_sub(minIdx);
   ny_min =  ny_sub(minIdx);

   rx_min = rx(minIdx);
   ry_min = ry(minIdx);
   rx_sub_min = rx_sub(minIdx);
   ry_sub_min = ry_sub(minIdx);


   qx_min =  qx_sub(minIdx) ;
   qy_min =  qy_sub(minIdx) ;

   qx_min2 = qx_sub2(minIdx) ;
   qy_min2 = qy_sub2(minIdx) ;


   Chaix_min = Chaix(minIdx) ;
   Chaiy_min = Chaiy(minIdx) ;

   Chaix_sub_min = Chaix_sub(minIdx) ;
   Chaiy_sub_min = Chaiy_sub(minIdx) ;


   sign = zeros( length(d_min), 1, 'gpuArray');
   sign(d_midx == 1) = -1;
   sign(d_midx == 2) = 1;

   coldirCheck = -sign.*( nx_min .* Chaix_min + ny_min .* Chaiy_min ) <= 0 & sign.*( nx_min .* Chaix_sub_min + ny_min .* Chaiy_sub_min ) <= 0 ;

   % rm_idx = ( ( (d_min.* D2_min) < d_col & (d_min.* D2_min) > 0 ) | (d_min.*(~D2_min) < d_col & d_min.*(~D2_min) > 0 ) ) & (~coldirCheck);

   rm_idx = d_min < d_col & ~coldirCheck;


    nx_min = nx_min(rm_idx);
    ny_min = ny_min(rm_idx);

    Chaix_min = Chaix_min(rm_idx);
    Chaiy_min = Chaiy_min(rm_idx);

    Chaix_sub_min = Chaix_sub_min(rm_idx);
    Chaiy_sub_min = Chaiy_sub_min(rm_idx);

    qx_min = qx_min(rm_idx );
    qy_min = qy_min(rm_idx);
    qx_min2 = qx_min2(rm_idx);
    qy_min2 = qy_min2(rm_idx);

    rx_min = rx_min(rm_idx);
    ry_min = ry_min(rm_idx);
    rx_sub_min = rx_sub_min(rm_idx);
    ry_sub_min = ry_sub_min(rm_idx);

    sign = sign(rm_idx);        



    ij = ij(rm_idx,:);


    b = - (f_col) * sign.*( (Chaix_min .* nx_min + Chaiy_min .* ny_min)  -  (Chaix_sub_min .* nx_min + Chaiy_sub_min .* ny_min)   );

    ij_list = ij;

    b_list =  b;

    n_list = [nx_min, ny_min];

    q1_list = [qx_min, qy_min];
    q2_list = [qx_min2, qy_min2];

    r1_list = [rx_min, ry_min];
    r2_list = [rx_sub_min, ry_sub_min];

    
%%

%%
    H = size(ij_list,1);
    n = [n_list; [n_list(:,2), n_list(:,1)] ];
    ijIdx = [ij_list; [ij_list(:,2), ij_list(:,1)] ];
    q_list = [q1_list; q2_list];
    % 
    % hold on
    % scatter( q1_list(:,1), q1_list(:,2), 'r','.')
    % 
    % scatter( q2_list(:,1), q2_list(:,2), 'g','.')

    ri_qij = q_list - Pcenter(ijIdx(:,1), :); 

    ri_qij_nij = ri_qij(:,1) .* n(:,2) - ri_qij(:,2) .* n(:,1);

    ri_qkl_x =  q_list(:,1) - Pcenter( ijIdx(:,1), 1)';
    ri_qkl_y =  q_list(:,2) - Pcenter( ijIdx(:,1), 2)';

    alphai = repmat( ijIdx(:,1)', [2*H,1]);
    alphaj = repmat( ijIdx(:,2)', [2*H,1]);
    alphak = repmat( ijIdx(:,1), [1,2*H]);
    alphal = repmat( ijIdx(:,2), [1,2*H]);
    
    overalpha = sum( alphai == reshape(overIdx,1,1,[]), 3) ~= 0;

    alpha = n * n' ./ ( Mass( ijIdx(:,1) )' ) + ri_qij_nij' .* ( -ri_qkl_y .* n(:,1) + ri_qkl_x .* n(:,2)  ) ./ ( I( ijIdx(:,1) )' );
    
    alpha(overalpha) = 0;

    ai = repmat( ij_list(:,1)' , [H,1]);
    aj = repmat( ij_list(:,2)' , [H,1]);
    ak = repmat( ij_list(:,1), [1,H] );
    al = repmat( ij_list(:,2), [1,H] );

    a = zeros(H,H, 'gpuArray');

    aidx1 = (ai == ak) & (aj == al);
    aidx2 = (ai == ak) & (aj ~= al);
    aidx3 = (aj == ak) & (ai ~= al);
    aidx4 = (aj ~= ak) & (ai == al);
    aidx5 = (ai ~= ak) & (aj == al);
    aidx6 = (ak ~= ai) & (al ~= aj);

    av1_1 =   alphai == reshape( ai( aidx1 ), 1,1,[]) & alphak == reshape( ai( aidx1 ), 1,1,[]) & alphaj == reshape( aj( aidx1 ), 1,1,[]) & alphal == reshape( aj( aidx1 ), 1,1,[]); 

    av1_2 =   alphai == reshape( aj( aidx1 ), 1,1,[]) & alphak == reshape( aj( aidx1 ), 1,1,[]) & alphaj == reshape( ai( aidx1 ), 1,1,[]) & alphal == reshape( ai( aidx1 ), 1,1,[]);  

    av2 =  alphai == reshape( ai( aidx2 ), 1,1,[]) & alphak == reshape( ai( aidx2 ), 1,1,[]) & alphaj == reshape( aj( aidx2 ), 1,1,[]) & alphal == reshape( al( aidx2 ), 1,1,[]);  

    av3 =  alphai == reshape( aj( aidx3 ), 1,1,[]) & alphak == reshape( aj( aidx3 ), 1,1,[]) & alphaj == reshape( ai( aidx3 ), 1,1,[]) & alphal == reshape( al( aidx3 ), 1,1,[]);  

    av4 =  alphai == reshape( ai( aidx4 ), 1,1,[]) & alphak == reshape( ak( aidx4 ), 1,1,[]) & alphaj == reshape( aj( aidx4 ), 1,1,[]) & alphal == reshape( ai( aidx4 ), 1,1,[]);  

    av5 =  alphai == reshape( aj( aidx5 ), 1,1,[]) & alphak == reshape( ak( aidx5 ), 1,1,[]) & alphaj == reshape( ai( aidx5 ), 1,1,[]) & alphal == reshape( aj( aidx5 ), 1,1,[]);  


    [rav1_1,cav1_1,~] = ind2sub(size(av1_1),find(av1_1));
    [rav1_2,cav1_2,~] = ind2sub(size(av1_2),find(av1_2));
    [rav2,cav2,~] = ind2sub(size(av2),find(av2));
    [rav3,cav3,~] = ind2sub(size(av3),find(av3));
    [rav4,cav4,~] = ind2sub(size(av4),find(av4));
    [rav5,cav5,~] = ind2sub(size(av5),find(av5));

    a(aidx1) = diag( alpha( rav1_1, cav1_1) ) +  diag( alpha( rav1_2, cav1_2) );
    a(aidx2) = diag( alpha( rav2, cav2 ) );
    a(aidx3) = diag( alpha( rav3, cav3 ) );
    a(aidx4) = diag( alpha( rav4, cav4 ) );
    a(aidx5) = diag( alpha( rav5, cav5 ) );
    a(aidx6) = 0;

    f = a\b_list;

    Pidx = unique( ij_list );

    ijsum = sum(ij_list,2);

    % coldirCheck = ( rx_min .* Chaix_min + ry_min .* Chaiy_min ) >= 0 & ( rx_sub_min .* Chaix_sub_min + ry_sub_min .* Chaiy_sub_min ) >= 0 ;
    % 
    % f(coldirCheck) = 0;

    % hold on
%%
    % VF = V;
    % omegaF = omega;
    for i = 1:1:length(Pidx)
    % for i = 1:1:1
        idx1 = ij_list(:,1) == Pidx(i);
        idx2 = ij_list(:,2) == Pidx(i);
        idx = idx1 | idx2;

        n_sub = n_list(idx,:);
        % n_sub2 = - n_sub;

        rsub = zeros( sum(idx), 2,  'gpuArray');


        rsub(  idx1(idx),: ) = r1_list(idx1,:);
        rsub( idx2(idx),: ) = r2_list(idx2,:);

        dvect =dot(n_sub, Pcenter( ( ijsum(idx) - Pidx(i) ), :) - Pcenter( Pidx(i),:), 2 );

        n_sub( dvect >= 0,:) = - n_sub(dvect >= 0,:);

        
        % n_sub2(dvect<0, :) = -n_sub2(dvect<0,:);
        % 
        % n_sub
        % n_sub2
        
        VF(Pidx(i),:) = VF(Pidx(i),:) + sum( abs(f(idx)) .* n_sub, 1) / Mass(Pidx(i));
        % Pidx(i)
        % ( f(idx) .* n_sub) / Mass(Pidx(i))
        % (  f(idx).* (rsub(:,1).*n_sub(:,2) - rsub(:,2).*n_sub(:,1)) ) / I(Pidx(i))


        omegaF(Pidx(i)) = omegaF(Pidx(i)) + sum(  abs(f(idx)).* (rsub(:,1).*n_sub(:,2) - rsub(:,2).*n_sub(:,1)), 1 ) / I(Pidx(i)); 


    end




end


function idx = upperTriangularIndexing(nump)
    % Check if the input is a square matrix

    
    idx = reshape( (1:1:nump^2), nump, nump);
    idx = triu( idx, 1);
    idx = reshape( idx', [], 1);
    idx(idx==0) = [];
    

end

function idx = lowerTriangularIndexing(nump)

    idx = reshape( (1:1:nump^2), nump, nump);
    idx = tril( idx, -1);
    idx = reshape( idx, [], 1);
    idx(idx==0) = [];

end