function a_d = drag_value(h,v,n_1,n_2,n_3)
    
    %required inputs (3u sattelite)
    %h=height in meters
    %v=velocity in m/s (3d vector)
    %n_1 = normal vector of side 1 in body frame (10cm x 10 cm) (3d vector)
    %n_2 = normal vector of side 2 in body frame (10cm x 30 cm) (3d vector)
    %n_3 = normal vector of side 3 in body frame (10cm x 30 cm) (3d vector)
    %mass of sattelite (taken to be 4 kgs here by default) 

    consts_vals=consts(h);
    A=0.0001*abs(dot(n_1,v))+0.0003*(abs(dot(n_2,v))+abs(dot(n_3,v)));
    rho=consts_vals(2)*exp((consts_vals(2)-h)/consts_vals(3));
    a_d = 0.125*rho*norm(v)*2.2*A.*v;
    return;
end    


function consts_out = consts(h)
    vals=[0,25,0,1.225,0,8.44;
        25,30,25,3.899,-2,6.49;
        30,35,30,1.774,-2,6.75;
        35,40,35,8.279,-3,7.07;
        40,45,40,3.972,-3,7.47;
        45,50,45,1.995,-3,7.83;
        50,55,50,1.057,-3,7.95;
        55,60,55,5.821,-4,7.73;
        60,65,60,3.206,-4,7.29;
        65,70,65,1.718,-4,6.81;
        70,75,70,8.77,-5,6.33;
        75,80,75,4.178,-5,6;
        80,85,80,1.905,-5,5.7;
        85,90,85,8.337,-6,5.41;
        90,95,90,3.396,-6,5.38;
        95,100,95,1.343,-6,5.74;
        100,110,100,5.297,-7,6.15;
        110,120,110,9.661,-8,8.06;
        120,130,120,2.438,-8,11.6;
        130,140,130,8.484,-9,16.1;
        140,150,140,3.845,-9,20.6;
        150,160,150,2.07,-9,24.6;
        160,180,160,1.224,-9,26.3;
        180,200,180,5.464,-10,33.2;
        200,250,200,2.789,-10,38.5;
        250,300,250,7.248,-11,46.9;
        300,350,300,2.418,-11,52.5;
        350,400,350,9.158,-12,56.4;
        400,450,400,3.725,-12,59.4;
        450,500,450,1.585,-12,62.2;
        500,600,500,6.967,-13,65.8;
        600,700,600,1.454,-13,79;
        700,800,700,3.614,-14,109;
        800,900,800,1.17,-14,164;
        900,1000,900,5.245,-15,225;
        1000,10000000,1000,3.019,-15,268];
       for i = 1:length(vals)
            if h<=vals(i,2) && h>=vals(i,1) 
                consts_out = [vals(i,3),(vals(i,4))*10^(vals(i,5)),vals(i,6)];
                break;
    
            else 
                consts_out=[-1,-1,-1];    
            end    
       end     
    return
end