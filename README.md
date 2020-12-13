# Mathematical-computationsclear variables;
clc;

rhop = 2700;  % mass density
E = 70e9;     % youngs modulus
nu = 0.33;    % poisson ratio
xi = 0.04;    %damping coefficient
D = 150e-3;   %panel position
a = 0.5;    % panel width (x)
b = 0.5;    % panel length (y)
t = 3e-3; % panel thickness
Dp = 3.1030e+03;
%Dp = (E*t^3)/12*(1-nu);% rigidity
d = 0.2e-3; % hole diameter
sigma = 0.01;    % perforation ratio
rhoa = 1.204;   % density of air
ca = 343;       % speed of sound in air
p = 1;          %external pressure
rhopt = 1;
% rhopt = rhop*t;
%f = 21.2114;    % frequency of excitation
f_space =  600:5:800;
%len = size(f_space);
alpha = zeros(size(f_space));
for i_f = 1:length(f_space)   
    disp(i_f);
    f = f_space(i_f);

    omega = f*2*pi;
    % find Z0r and Z0i
    Z0r = rhoa*ca*(0.147*t/(d^2))*(sqrt(9+(100*d^2*f/32)) + 1.768*sqrt(f)*(d^2/t));
    Z0i = rhoa*ca*1.847*f*t*(1+(1/sqrt(9+50*d^2*f))+0.85*d/t);
    Z0 = Z0r+1i*Z0i;
    Z0_bar = Z0/sigma;



    M = 2;
    N = 2;
    A = zeros(M*N);
    C = zeros((M*N),1);

    i = 1;


    for mp = 1:M %row of A
        for np = 1:N
            j = 1;
            for m = 1:M %column of A
                for n = 1:N

                    % T1 

                    lambda_m = m*pi/a;
                    lambda_n = n*pi/b;
                    omega_mn = sqrt(Dp/rhopt)*(lambda_m^2+lambda_n^2);
                    %find eta_mpnp
                    
                    
                    eta_mpnp = a*b/4;
                    %fun_eta = @(x,y) ((sin(mp*pi*x/a)).^2).*((sin(np*pi*y/b)).^2);
                    
                    %eta_mpnp = (integral2(fun_eta,0,a,0,b));


                    T1 = rhopt*(xi*omega_mn*omega+1i*(omega^2-omega_mn^2))/omega*eta_mpnp;

                    % T2
                    % find Z0_bar
                    %Z0_bar = Z0/sigma; 
                    sum_uw = 0;
                    for u = 1:M
                        for w = 1:N
                            % find Za_uw, beta_uw, gamma_mnuw, gamma_mpnpuw
                            mu_uw = sqrt(((u*pi/a)^2 + (w*pi/b)^2 - (omega/ca)^2));
                            Za_uw = 1i*rhoa*omega*(coth(mu_uw*D))/mu_uw;
                             
                            beta_uw = a*b/4;
%                             if u==w
%                                 beta_uw = a*b/4;
%                             else
%                                 beta_uw = 0;
%                             end
%                             fun_beta = @(x,y) ((cos(u*pi*x/a)).^2).*((cos(w*pi*y/b)).^2);
%                             beta_uw = integral2(fun_beta,0,a,0,b);
                            fun_gamma_mn = @(x,y) cos(u*pi*x/a).*cos(w*pi*y/b).*sin(lambda_m*x).*sin(lambda_n*y);
                            fun_gamma_mpnp = @(x,y) cos(u*pi*x/a).*cos(w*pi*y/b).*sin(mp*pi*x/a).*sin(np*pi*y/b);
                            gamma_mnuw = integral2(fun_gamma_mn,0,a,0,b);
                            gamma_mpnpuw = integral2(fun_gamma_mpnp,0,a,0,b);
                            sum_uw = sum_uw + Z0_bar*Za_uw/(beta_uw*(Z0_bar+Za_uw))*gamma_mnuw*gamma_mpnpuw;
                        end
                    end
                    T2 = sum_uw;
                    A(i,j) = T1 - T2;
                    j = j+1;
                end
            end
            i = i+1;
        end
    end

    k=1;
    for mp = 1:M                %RHS constant term
        for np = 1:N 
            Z0_bar = Z0/sigma;
            mu_00 = sqrt(((0*pi/a)^2 + (0*pi/b)^2 - (omega/ca)^2));
            Za_00 = 1i*rhoa*omega*(coth(mu_00*D))/mu_00;
            if mod(mp,2)*mod(np,2)==0
                epsilon_mpnp = 0;
            else
                epsilon_mpnp = 4*a*b/(pi*pi*mp*np);
            end
            
            %fun_epsilon = @(x,y) sin(mp*pi*x/a).*sin(np*pi*y/b);
           % epsilon_mpnp = integral2(fun_epsilon,0,a,0,b);
            T3 = (p*epsilon_mpnp*Z0_bar)/(Z0_bar+Za_00);
            C(k,1) = T3;
            k=k+1;
        end
    end


    B = A\C;             % B Matrix

    %find V and PD 


        v_fun = @(x,y) v_spatial(B,x,y,M,N,a,b);
        pd_fun = @(xp,yp) pd_spatial(v_fun, xp, yp, p, Z0_bar, a, b, omega, ca, rhoa, M, N, Z0, D);
    %v_fun(0.1,0.1)
   

    %pd_fun(xp,yp)
    %PD(i_x,j_y)

    integral_v = integral2(v_fun,0,a,0,b);

    p_fun = @(xp,yp) (p-pd_fun(xp,yp))/Z0_bar;
    integral_p = integral2(p_fun,0,a,0,b);

    vd_bar = (integral_v + integral_p)/(a*b);

    Z_bar = p/(rhoa*ca*vd_bar);

    %disp(Z_bar);

    alpha(i_f) = 4*abs(real(Z_bar))/((1+real(Z_bar))^2 + (imag(Z_bar))^2);
    
    plot(f_space,alpha);
    

end


function PD = pd_spatial(v_fun, xp, yp, p, Z0_bar, a, b, omega, ca, rhoa, M, N, Z0, D) 
    
    [x_m,y_n] = size(xp);
    PD = zeros(x_m,y_n);
    for i_x = 1:x_m
        for j_y = 1:y_n
            sum = 0;
            for u = 1:M

                for w = 1:N
                    integrand = @(x,y) (v_fun(x,y)+(p/Z0)).*cos(u.*pi.*x./a).*cos(w.*pi.*y./b);
                    integral_xy = integral2(integrand,0,a,0,b);
                    %integral_xy

                    mu_uw = sqrt((u*pi/a)^2 + (w*pi/b)^2 + (omega/ca)^2);
                    Za_uw = 1i*rhoa*omega*(coth(mu_uw*D))/mu_uw;
                    
                     beta_uw = a*b/4;
%                     if u==w
%                         beta_uw = a*b/4;
%                     else
%                         beta_uw = 0;
%                     end

%                      fun_beta = @(x,y) ((cos(u.*pi.*x./a)).^2).*((cos(w.*pi.*y./b)).^2);
%                    beta_uw = integral2(fun_beta,0,a,0,b);

                    sum = sum + ((Z0_bar*Za_uw/(beta_uw*(Z0_bar+Za_uw)))*integral_xy*cos(u*pi*xp(i_x,j_y)/a)*cos(w*pi*yp(i_x,j_y)/b));
                end
            end
            PD(i_x,j_y) = sum;
        end
    end
end




function v = v_spatial(B, x, y, M, N, a, b)
   
	X = @(m,x) sin(m.*pi./a.*x);
	Y = @(n,y) sin(n.*pi./b.*y);
    %size(x)
    [x_m,y_n] = size(x);
    v = zeros(x_m,y_n);
    for i_x = 1:x_m
        for j_y = 1:y_n
            
            XY = zeros(1,M*N);
            i = 1;
            for m = 1:M
                for n = 1:N
                    %X(m,x)
                   % Y(n,y)
                    XY(1,i) = X(m,x(i_x,j_y)).*Y(n,y(i_x,j_y));
                    i = i+1;

                end
            end
            v(i_x,j_y) = XY*B;
        end
    end
    %disp(v);
           
	
end

