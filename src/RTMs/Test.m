%%Test

function Lb = Planck(wl,Tb,em)

    c1 = 1.191066e-22;
    c2 = 14388.33;
    if nargin<3
        em = ones(size(Tb));
    end
    
    Lb = em.* c1*(wl*1e-9).^(-5)./(exp(c2./(wl*1e-3*Tb))-1);
    
end    

function Lb = PlanckSI(wl,Tb,em)

    c1 = 1.191066e-16;
    c2 = 0.01438833;
    if nargin<3
        em = ones(size(Tb));
    end
    
    Lb = em.* c1*(wl).^(-5)./(exp(c2./(wl.*Tb))-1);
    
end    

function int = Sint(y,x)

    % Simpson integration
    % x and y must be any vectors (rows, columns), but of the same length
    % x must be a monotonically increasing series
    
    % WV Jan. 2013, for SCOPE 1.40
    
    nx   = length(x);
    if size(x,1) == 1
        x = x';
    end
    if size(y,1) == size(x,1)
        y = y';
    end

    step = x(2:nx) - x(1:nx-1);
    mean = .5 * (y(:,1:nx-1) + y(:,2:nx));
    int  = mean * step;
end




Tu = 30.8;
Tc_K = Tu + 273.15;
em = 0.98;
%wl = linspace(1000, 50000, 49000); % creates a range with no of total values in between as specified by the last value
wl_a = 2500:1:50000;
%wl_m = wl_a*1e-9;
lb = pi*Planck(wl_a,Tc_K,em);
int_model = 1E-3*Sint(lb,wl_a) 
%int_h_model = pi*int_model
%int_trap_model = trapz (wl_a,lb)
%int_h_model = pi*int_trap_model

%lb2 = PlanckSI(wl_m,Tc_K,em);
%int2 = Sint(lb,wl_m) 
%int_h2 = pi*int



function H      =   Stefan_Boltzmann(T_C,constants)

C2K     = constants.C2K;
sigmaSB = constants.sigmaSB;

H       = sigmaSB*(T_C + C2K).^4;
return

end




Tu = 30.82;
Tc_K = Tu + 273.15;
em = 0.98;
wl_b = 2500:1:500000;
wl_m = wl_b*1e-9;
lb2 = PlanckSI(wl_m,Tc_K,em);
%int_trap = trapz (wl_m,lb);
int2 = Sint(lb2,wl_m) 
%int3 = simps(wl_m,lb2)
int_h2 = pi*int2
%int_h3 = pi*int3

emb = 1;
H = Stefan_Boltzmann(Tu,constants);
R = em*H
Rb= emb*H

test_emis = int_h2/Rb 


Ts = 23.08;
Ts_K = Ts + 273.15;
wl_c = 1000:1:500000;

em_s = 0.94;
R_s = Planck(wl_c,Ts_K,em_s);
int_s = 1E-3*Sint(R_s,wl_c) 
soil = pi*int_s

Ts = 23.08;
Ts_K = Ts + 273.15;
wl_c = 2500:1:500000;
wl_m_si = wl_b*1e-9;
em_s = 0.94;
R_s = PlanckSI(wl_m_si,Ts_K,em_s);
int_si = Sint(R_s,wl_m_si) 
soil_si = pi*int_si


emb_s = 1;
H_s = Stefan_Boltzmann(Ts,constants);
Rb_s= emb_s*H_s
R_s = em_s*H_s
test_emis_soil_bb = soil_si/Rb_s 
test_emiss_soil =  soil_si/R_s 


%Per Wavelength Estimation

wl_start = 2500;
wl_m_si_st = wl_start*1e-9;
R_start_si = pi*PlanckSI(wl_m_si_st,Ts_K,em_s)
R_start_s = pi*Planck(wl_start,Ts_K,em_s)

wl_start1 = 8000;
wl_m_si_st1 = wl_start1*1e-9;
R_start_si1 = pi*PlanckSI(wl_m_si_st1,Ts_K,em_s)
R_start_s1 = pi*Planck(wl_start1,Ts_K,em_s)

wl_start2 = 500000;
wl_m_si_st2 = wl_start2*1e-9;
R_start_si2 = pi*PlanckSI(wl_m_si_st2,Ts_K,em_s)
R_start_s2 = pi*Planck(wl_start2,Ts_K,em_s)