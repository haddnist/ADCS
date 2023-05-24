%This is the low fidelity model for the Solar Radiation pressure
%Here we are assuming sun facing area as constant(independent of attitude)
%for shadow region we are assuming earth as spherical
%Work to be done -  add SRP variation based on attitude
%                   
%                   
%Inputs
%           s = time in seconds after mission has started(on Jan 1,2024)
%           u_EC = unit vector from the earth's center to Cube-sat 
%           r_EC = distance of cubesat from earth
%           m = mass of cubesat
function a = srp_low(s,u_EC,r_EC,m)
%defining some constants
sc=1362; %solar constant
c = 299792458; %speed of light
c_SRP =1.8; %Solar drag coefficient for satellites
area = 3 * 0.01;
r_E =  6.3781366 * 10^6; %radius of earth
r_S = 6.957 *10^8; %radius of sun
pi = 3.141592653; 
au = 149597870691; %astromical unit in metre

%calculating the unit vector in the Earth-to-Sun direcion and distance
%between earth and sun
T_jul = JD(s);
T = (T_jul-2451545)/36525;
phi =  280.460 + 36000.771*T; %mean longitude of the sun
M = 357.5277233 + 35999.05034*T; %mean anomaly of the sun
phi_ec = phi + (1.914666471 * sin(M)) + (0.019994643 * sin(2*M)); %Longitude of the ecliptic
e =23.439291 - (0.0130042 * T); %Obliquity of the ecliptic
u_ES = [cos(phi_ec),(cos(e)*sin(phi_ec)),(sin(e)*sin(phi_ec))];
%for distance, reference Jean Meeus - Astronomical Algorithm pg 164
%centre of Sun
cen = (1.914602 - (0.004817 * T) - (0.000014 * T *T))*sin(M) + (0.01993 - 0.000101*T)*sin(2*M) + 0.000289 * sin(3*M);

%true anomaly is
neu = cen + M;
%eccentricity of earth orbit
ex = 0.016708634 - 0.000042037*T - 0.0000001267*T*T;
r_ES = ((1.000001018 *(1- (ex)^2))/(1+ex*cos(neu)))*au; %distance of earth from the sun

%calculating shadow region
sh = 0;
umbra = (r_E/(r_S-r_E))*r_ES;
theta = asin(r_E/umbra); %half angle of the shadow cone
t = quadsolve(((umbra/r_EC)+1),(2/tan(theta)),((umbra/r_EC)-1));
if (t(1)>=0)
    alpha = 2*atan(t(1)); %alpha is max angle range for a given r_EC
else
    alpha = 2*atan(t(2));
end
if (r_EC<=umbra)
    beta = pi - acos(dot(u_ES,u_EC));
    if(beta<=alpha)
        sh=-1;
    end
end


%calclating u_CS, unit vector from centre of satellite to centre of the
%sun
R_CS = (r_ES * u_ES) - (r_EC * u_EC); %direction vector from cubesat to centre of the sun
u_CS = R_CS/norm(R_CS);

%calculating Solar radiation pressure
P = sc/(c*(norm(R_CS)^2));


%calculating the acceleration
F = -(P*c_SRP*u_CS)*area;

if(sh~=-1)
    a = F/m;
else
    a=0;
end
end

function roots = quadsolve(a,b,c)
roots = zeros(2,1);
delta = (b.^2)-(4*a*c);
if (delta == 0)
   roots(1) = (-b./(2*a));
   roots(2) = (-b./(2*a));
else
  roots(1) = ((-b+sqrt(delta))./(2*a));
  roots(2) = ((-b-sqrt(delta))./(2*a));
end
end
