function R_rs = lee(X,wl)

    C_chl = X(1);
    C_nap = X(2);
    C_cdom=X(3);
    tetaw = 30;
    tetav = 0;
    z = X(4);

    S_cdom = 0.0157;
    S_nap = 0.0106;
    a_nap_e = 0.0048; %440
    b_bphy_e =0.00038; %542
    Y_phy = 0.681;
    Y_nap = 0.681;
    b_bnap_e = 0.0054; %542

    %pf = load(wl);
    pf = evalin('base','hico_grey_sand');   %importation de de hico_aw
    pf = pf(:,2);
    %pf = [ pf ; ones(100,1)*pf(length(pf))];

    %data = load('hico_aw.txt');
    data = evalin('base','hico_aw');   %importation de de hico_aw
    lambda = data(:,1); % longueur d'onde 
    a_w = data(:,2);     % absorption de l'eau

    %data = load('hico_aphya.txt');
    data2 = evalin('base','hico_aphya');
    a_phy_e = data2(:,2);

    a_phy = C_chl*a_phy_e;
    a_nap = C_nap*a_nap_e.*exp((-S_nap*(lambda-440)));
    a_cdom = C_cdom.*exp(-S_cdom*(lambda - 440));

    a_tot = a_nap + a_cdom + a_w + a_phy;

    %b_bw = 0.00144*(lambda/500).^(-4.32);
    %b_bw = load('hico_bbw.txt');
    b_bw = evalin('base','hico_bbw');

    b_bw = b_bw(:,2);
    b_bp = C_chl*b_bphy_e*(542./lambda).^(Y_phy) + C_nap*b_bnap_e*(542./lambda).^(Y_nap);
    b_b = b_bw + b_bp;

    K = a_tot + b_b; % coefficient d'attténuation diffuse 
    u = b_b./(a_tot + b_b);
    u_p = b_bp./(a_tot + b_b);
    D_uc = 1.03*(1+2.4.*u).^(0.5);
    D_ub = 1.04*(1+5.4*u).^(0.5);
    g_p = 0.184.*(1-0.602.*exp(-3.852.*u_p));
    g_w = 0.115;
    r_rs_dp = g_w.*b_bw./(a_tot + b_b) + g_p.*b_bp./(a_tot +b_b);
    r_rsC = r_rs_dp.*(1-exp(-(1./cosd(tetaw) + D_uc./cosd(tetav)).*K.*z));
    r_rsB = 1/pi .* pf .* exp( -(1./cosd(tetaw) + D_ub./cos(tetav)).*K.*z);
    
    r_rs = r_rsC + r_rsB;

    R_rs = (0.52.*r_rs)./(1-1.56.*r_rs); % Réflectance de télédection

end