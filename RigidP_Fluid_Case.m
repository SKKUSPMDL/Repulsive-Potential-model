afunction [ uelem, velem, pelem, dx, dy, xlimit, ylimit,  L, H, Lp, Hp, density, gamma, dt, Node, Element, U, Wall, viscosity] = RigidP_Fluid_Case(c)



if c == 1

    Node = readmatrix("50um_101325_20000.xlsx",'Sheet', 'Nodes');
    nnode = size(Node,1);
    
    Element = readmatrix("50um_101325_20000.xlsx",'Sheet', 'Elements');
    Element = Element + 1;
    nelem = size(Element,1);
    
    Data = readmatrix("50um_101325_20000.xlsx",'Sheet', 'Data');
    
    % dw = 1e-6;
    


    Wall1 = [0, -150*1e-6, -250*1e-6, -150*1e-6, 0, 1, 250*1e-6];

    Wall2 = [-250*1e-6, -150*1e-6, -250*1e-6, 250*1e-6, 1, 0, 400*1e-6];

    Wall3 = [0, 250*1e-6, -250*1e-6, 250*1e-6, 0, -1, 100*1e-6];
    
    Wall4 = [0, 75*1e-6, 0, 250*1e-6, -1, 0, 175*1e-6];

    Wall5 = [0, 75*1e-6, 100*1e-6, 75*1e-6, 0, -1, 100*1e-6];

    Wall6 = [0, 25*1e-6, 100*1e-6, 25*1e-6, 0, 1, 100*1e-6];

    Wall7 = [0, 25*1e-6, 0, -150*1e-6, -1, 0, 175*1e-6];

    Wall = [Wall1; Wall2; Wall3; Wall4; Wall5; Wall6; Wall7];

    





elseif c == 2
    
    Node = readmatrix("100um_101325_20000.xlsx",'Sheet', 'Nodes');
    nnode = size(Node,1);
    
    Element = readmatrix("100um_101325_20000.xlsx",'Sheet', 'Elements');
    Element = Element + 1;
    nelem = size(Element,1);
    
    Data = readmatrix("100um_101325_20000.xlsx",'Sheet', 'Data');

    

    Wall1 = [0, -150*1e-6, -250*1e-6, -150*1e-6, 0, 1, 100*1e-6];

    Wall2 = [-250*1e-6, -150*1e-6, -250*1e-6, 250*1e-6, 1, 0, 400*1e-6];

    Wall3 = [0, 250*1e-6, -250*1e-6, 250*1e-6, 0, -1, 100*1e-6];
    
    Wall4 = [0, 100*1e-6, 0, 250*1e-6, -1, 0, 150*1e-6];

    Wall5 = [0, 100*1e-6, 100*1e-6, 100*1e-6, 0, -1, 100*1e-6];

    Wall6 = [0, 0*1e-6, 100*1e-6, 0*1e-6, 0, 1, 100*1e-6];

    Wall7 = [0, 0*1e-6, 0, -150*1e-6, -1, 0, 150*1e-6];

    % Wall8 = [100e-6, 0, 100e-6, 100e-6, -1, 0, 150*1e-6];

    Wall = [Wall1; Wall2; Wall3; Wall4; Wall5; Wall6; Wall7];

elseif c == 3

    Node = readmatrix("150um_101325_20000.xlsx",'Sheet', 'Nodes');
    nnode = size(Node,1);
    
    Element = readmatrix("150um_101325_20000.xlsx",'Sheet', 'Elements');
    Element = Element + 1;
    nelem = size(Element,1);
    
    Data = readmatrix("150um_101325_20000.xlsx",'Sheet', 'Data');
    
    
    
    dw = 1e-6;

    

    Wall1 = [0, -150*1e-6, -250*1e-6, -150*1e-6, 0, 1, 100*1e-6];

    Wall2 = [-250*1e-6, -150*1e-6, -250*1e-6, 250*1e-6, 1, 0, 400*1e-6];

    Wall3 = [0, 250*1e-6, -250*1e-6, 250*1e-6, 0, -1, 100*1e-6];
    
    Wall4 = [0, 125*1e-6, 0, 250*1e-6, -1, 0, 125*1e-6];

    Wall5 = [0, 125*1e-6, 100*1e-6, 125*1e-6, 0, -1, 100*1e-6];

    Wall6 = [0, -25*1e-6, 100*1e-6, -25*1e-6, 0, 1, 100*1e-6];

    Wall7 = [0, -25*1e-6, 0, -150*1e-6, -1, 0, 125*1e-6];

    Wall = [Wall1; Wall2; Wall3; Wall4; Wall5; Wall6; Wall7];


elseif c == 4

    Node = readmatrix("200um_101325_20000.xlsx",'Sheet', 'Nodes');
    nnode = size(Node,1);
    
    Element = readmatrix("200um_101325_20000.xlsx",'Sheet', 'Elements');
    Element = Element + 1;
    nelem = size(Element,1);
    
    Data = readmatrix("200um_101325_20000.xlsx",'Sheet', 'Data');
    
    
    
    dw = 1e-6;

    

    Wall1 = [0, -150*1e-6, -250*1e-6, -150*1e-6, 0, 1, 100*1e-6];

    Wall2 = [-250*1e-6, -150*1e-6, -250*1e-6, 250*1e-6, 1, 0, 400*1e-6];

    Wall3 = [0, 250*1e-6, -250*1e-6, 250*1e-6, 0, -1, 100*1e-6];
    
    Wall4 = [0, 150*1e-6, 0, 250*1e-6, -1, 0, 125*1e-6];

    Wall5 = [0, 150*1e-6, 100*1e-6, 150*1e-6, 0, -1, 100*1e-6];

    Wall6 = [0, -50*1e-6, 100*1e-6, -50*1e-6, 0, 1, 100*1e-6];

    Wall7 = [0, -50*1e-6, 0, -150*1e-6, -1, 0, 125*1e-6];

    Wall = [Wall1; Wall2; Wall3; Wall4; Wall5; Wall6; Wall7];
%%%%%%%%%%%%%%%%%%%%%%%%%%55
elseif c == 15

    Node = readmatrix("15um_101325_20000.xlsx",'Sheet', 'Nodes');
    nnode = size(Node,1);
    
    Element = readmatrix("15um_101325_20000.xlsx",'Sheet', 'Elements');
    Element = Element + 1;
    nelem = size(Element,1);
    
    Data = readmatrix("15um_101325_20000.xlsx",'Sheet', 'Data');
    
    
    
    dw = 1e-6;

    

    Wall1 = [0, -150*1e-6, -250*1e-6, -150*1e-6, 0, 1, 100*1e-6];

    Wall2 = [-250*1e-6, -150*1e-6, -250*1e-6, 250*1e-6, 1, 0, 400*1e-6];

    Wall3 = [0, 250*1e-6, -250*1e-6, 250*1e-6, 0, -1, 100*1e-6];
    
    Wall4 = [0, 57.5*1e-6, 0, 250*1e-6, -1, 0, 125*1e-6];

    Wall5 = [0, 57.5*1e-6, 100*1e-6, 57.5*1e-6, 0, -1, 100*1e-6];

    Wall6 = [0, 42.5*1e-6, 100*1e-6, 42.5*1e-6, 0, 1, 100*1e-6];

    Wall7 = [0, 42.5*1e-6, 0, -150*1e-6, -1, 0, 125*1e-6];

    Wall = [Wall1; Wall2; Wall3; Wall4; Wall5; Wall6; Wall7];

elseif c == 25

    Node = readmatrix("25um_101325_20000.xlsx",'Sheet', 'Nodes');
    nnode = size(Node,1);
    
    Element = readmatrix("25um_101325_20000.xlsx",'Sheet', 'Elements');
    Element = Element + 1;
    nelem = size(Element,1);
    
    Data = readmatrix("25um_101325_20000.xlsx",'Sheet', 'Data');
    
    
    
    dw = 1e-6;

    

    Wall1 = [0, -150*1e-6, -250*1e-6, -150*1e-6, 0, 1, 100*1e-6];

    Wall2 = [-250*1e-6, -150*1e-6, -250*1e-6, 250*1e-6, 1, 0, 400*1e-6];

    Wall3 = [0, 250*1e-6, -250*1e-6, 250*1e-6, 0, -1, 100*1e-6];
    
    Wall4 = [0, 62.5*1e-6, 0, 250*1e-6, -1, 0, 125*1e-6];

    Wall5 = [0, 62.5*1e-6, 100*1e-6, 62.5*1e-6, 0, -1, 100*1e-6];

    Wall6 = [0, 37.5*1e-6, 100*1e-6, 37.5*1e-6, 0, 1, 100*1e-6];

    Wall7 = [0, 37.5*1e-6, 0, -150*1e-6, -1, 0, 125*1e-6];

    Wall = [Wall1; Wall2; Wall3; Wall4; Wall5; Wall6; Wall7]; 
elseif c == 35

    Node = readmatrix("35um_101325_20000.xlsx",'Sheet', 'Nodes');
    nnode = size(Node,1);
    
    Element = readmatrix("35um_101325_20000.xlsx",'Sheet', 'Elements');
    Element = Element + 1;
    nelem = size(Element,1);
    
    Data = readmatrix("35um_101325_20000.xlsx",'Sheet', 'Data');
    
    
    
    dw = 1e-6;

    

    Wall1 = [0, -150*1e-6, -250*1e-6, -150*1e-6, 0, 1, 100*1e-6];

    Wall2 = [-250*1e-6, -150*1e-6, -250*1e-6, 250*1e-6, 1, 0, 400*1e-6];

    Wall3 = [0, 250*1e-6, -250*1e-6, 250*1e-6, 0, -1, 100*1e-6];
    
    Wall4 = [0, 67.5*1e-6, 0, 250*1e-6, -1, 0, 125*1e-6];

    Wall5 = [0, 67.5*1e-6, 100*1e-6, 67.5*1e-6, 0, -1, 100*1e-6];

    Wall6 = [0, 32.5*1e-6, 100*1e-6, 32.5*1e-6, 0, 1, 100*1e-6];

    Wall7 = [0, 32.5*1e-6, 0, -150*1e-6, -1, 0, 125*1e-6];

    Wall = [Wall1; Wall2; Wall3; Wall4; Wall5; Wall6; Wall7];    
elseif c == 300

    Node = readmatrix("300um_101325_20000.xlsx",'Sheet', 'Nodes');
    nnode = size(Node,1);
    
    Element = readmatrix("300um_101325_20000.xlsx",'Sheet', 'Elements');
    Element = Element + 1;
    nelem = size(Element,1);
    
    Data = readmatrix("300um_101325_20000.xlsx",'Sheet', 'Data');
    
    
    
    dw = 1e-6;

    

    Wall1 = [0, -150*1e-6, -250*1e-6, -150*1e-6, 0, 1, 100*1e-6];

    Wall2 = [-250*1e-6, -150*1e-6, -250*1e-6, 250*1e-6, 1, 0, 400*1e-6];

    Wall3 = [0, 250*1e-6, -250*1e-6, 250*1e-6, 0, -1, 100*1e-6];
    
    Wall4 = [0, 200*1e-6, 0, 250*1e-6, -1, 0, 125*1e-6];

    Wall5 = [0, 200*1e-6, 100*1e-6, 200*1e-6, 0, -1, 100*1e-6];

    Wall6 = [0, -100*1e-6, 100*1e-6, -100*1e-6, 0, 1, 100*1e-6];

    Wall7 = [0, -100*1e-6, 0, -150*1e-6, -1, 0, 125*1e-6];

    Wall = [Wall1; Wall2; Wall3; Wall4; Wall5; Wall6; Wall7];  
elseif c == 350

    Node = readmatrix("350um_101325_20000.xlsx",'Sheet', 'Nodes');
    nnode = size(Node,1);
    
    Element = readmatrix("350um_101325_20000.xlsx",'Sheet', 'Elements');
    Element = Element + 1;
    nelem = size(Element,1);
    
    Data = readmatrix("350um_101325_20000.xlsx",'Sheet', 'Data');
    
    
    
    dw = 1e-6;

    

    Wall1 = [0, -150*1e-6, -250*1e-6, -150*1e-6, 0, 1, 100*1e-6];

    Wall2 = [-250*1e-6, -150*1e-6, -250*1e-6, 250*1e-6, 1, 0, 400*1e-6];

    Wall3 = [0, 250*1e-6, -250*1e-6, 250*1e-6, 0, -1, 100*1e-6];
    
    Wall4 = [0, 225*1e-6, 0, 250*1e-6, -1, 0, 125*1e-6];

    Wall5 = [0, 225*1e-6, 100*1e-6, 225*1e-6, 0, -1, 100*1e-6];

    Wall6 = [0, -125*1e-6, 100*1e-6, -125*1e-6, 0, 1, 100*1e-6];

    Wall7 = [0, -125*1e-6, 0, -150*1e-6, -1, 0, 125*1e-6];

    Wall = [Wall1; Wall2; Wall3; Wall4; Wall5; Wall6; Wall7];  
end


dt = 0.7e-10;

L = 100*1e-6;
H = 400*1e-6;

rho_Total = 1410;
rho_BN = 2100;

Lp = 40*1e-6;
Hp = 2.105*1e-6;

Hf_BN = Lp;

rf_BN = sqrt( ( pi * (Lp/2)^2 * Hp) / (pi*Hf_BN) );


viscosity = 3.2;


I = pi/2 * rf_BN^4;
E = 0.865 * 1e+12;
gamma = E*I;

density = abs(rho_BN - rho_Total) * pi * (rf_BN^2);

IdxMatrix = ( ( Node(:,2) - Data(:,1)' ).^2 + ( Node(:,3) - Data(:,2)' ).^2 ).^0.5;

[rw, col] = find( IdxMatrix < 1e-12 );


Data = [rw, Data];

xelem = zeros(nelem, 4);
yelem = zeros(nelem, 4);
uelem = zeros(nelem, 4);
velem = zeros(nelem, 4);
pelem = zeros(nelem, 4);
Uelem = zeros(nelem, 4);

for i = 1:1:4
 [dataidx, ~] = find( Data(:,1) == Element(:,i)' ); 

 xelem(:,i) = Data(dataidx, 2);
 yelem(:,i) = Data(dataidx, 3);
 uelem(:,i) = Data(dataidx, 4);
 velem(:,i) = Data(dataidx, 5);
 pelem(:,i) = Data(dataidx, 6);
 Uelem(:,i) = Data(dataidx, 7);

end

xlimit = [min(xelem,[], 2), max(xelem,[], 2)];
ylimit = [min(yelem,[], 2), max(yelem,[], 2)];

dx = xlimit(:,2) - xlimit(:,1);
dy = ylimit(:,2) - ylimit(:,1);

U = zeros(nnode, 1);

U( Data(:,1), 1) = Data(:,7);



end