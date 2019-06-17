function output = calceff(L1,m_max,wavelength)
    % L1 - number of grating strata (for "staircase" approximation)
    % m_max - max diffraction order index
    % Wavelength in nm
    
    ctr_sect=true;
    % Define parameters for grating structure and incident field:
    d=160;              % grating period (nm)
%     d = 444.4444;
    blaze = 29.5;    
%     blaze = 4;
    %   Gold = 'Au'
    %   Silicon = 'Si'
    %   Nickel = 'Ni'

    element_grating = 'Au';
    element_substrate = 'Si';
    % [grating_pmt, wave_choice] = opticalConstantCode(element_grating, wavelength);
    grating_pmt=4.; % grating permittivity

    h=d.*tand(blaze); % grating height
    
    yaw = 0.85 * pi/180;
    graze = 1.5 * pi/180;
%     yaw = 1.07*pi/180;
%     graze = 15*pi/180;
 
    % phi3 = angle between the grating normal (x1 axis) and the incident wave
    % vector's projection normal to the grating lines (i.e., x1, x2 plane
    % projection)
    phi3=atan(sin(yaw)/tan(graze));
    
    % theta3 = angle between the grating lines (x3 axis) and the incident wave
    % vector
    theta3=asin(sin(graze)/cos(phi3));

    %%
    if ctr_sect
        % Define stratum half-widths (in units of d) by center-sectioning.
        c1=acos(-1+2*((1:L1)-0.5)/L1)/(2*pi);    
    %     c1 = (((1:L1)./cosd(blaze)))/L1;
    %     c1=c1./max(c1);
    else
        % Define stratum half-widths using slightly more accurate averaging
        % method.
        c1=2*(1:L1-1)/L1-1;
        c1=diff([-0.25, (c1.*acos(c1)-sqrt(1-c1.^2))/(4*pi), 0])*L1;
    end

    % Construct grating. (Note: Two grating periods must be specified
    % because the grating can generally be biperiodic, although for this
    % example the second period (d22,d32) is irrelevant.)
    clear grating
    grating.pmt={1,grating_pmt}; % grating material permittivities

    % [grating_substrate, wave_choice] = opticalConstantCode(element_substrate, wavelength);
    % grating.pmt_sub_index=grating_substrate; % substrate permittivity index
    zeta = 30*pi/180;
    grating.pmt_sub_index=2; % substrate permittivity index
    grating.pmt_sup_index=1; % superstrate permittivity index (1 is a vacuum)
    grating.d21=d; % first grating period: x2 projection - d
    grating.d31=0; % first grating period: x3 projection - 0
    grating.d22=0; % second grating period: x2 projection - 0
    grating.d32=d; % second grating period: x3 projection - d
    grating.stratum={};
    clear stratum stripe
    stratum.type=1; % uniperiodic stratum
    stratum.thick=h/L1; % stratum thickness
    % The following h11, h12 spec indicates that the stratum's period
    % vector matches the first grating period (GD-Calc.pdf, equations 3.22
    % and 3.23).
    stratum.h11=1; % - 1
    stratum.h12=0; % - 0
    for l1=1:L1
    %     stripe.pmt_index=1; % first stripe's permittivity index
    %     stripe.c1=-c1(l1); % first stripe's boundary on positive side
    %     stratum.stripe{1}=stripe;
    %     stripe.pmt_index=2; % second stripe's permittivity index
    %     stripe.c1=c1(l1); % second stripe's boundary on positive side
    %     stratum.stripe{2}=stripe;
    %     grating.stratum{end+1}=stratum;
        stripe.pmt_index=1; % first stripe's permittivity index
        stripe.c1=(l1-.5)/L1; %-c1(l1); % first stripe's boundary on positive side
        stratum.stripe{1}=stripe;
        stripe.pmt_index=2; % second stripe's permittivity index
        stripe.c1=1; %c1(l1); % second stripe's boundary on positive side
        stratum.stripe{2}=stripe;
        grating.stratum{end+1}=stratum;
    end
    clear c1 stratum stripe

    % Define the indicent field.
    clear inc_field
    inc_field.wavelength=wavelength;
    % f2 and f3 are the grating-tangential (x2 and x3) coordinate projections
    % of the incident field's spatial-frequency vector. The grating-normal (x1)
    % projection is implicitly f1=-sin(theta3)*cos(phi3)/wavelength.
    inc_field.f2=sin(theta3)*sin(phi3)/wavelength;
    inc_field.f3=cos(theta3)/wavelength;

    % Specify which diffracted orders are to be retained in the calculations.
    % (m2 is zero because this is a uniperiodic grating - all diffraction
    % orders for m2~=0 have zero amplitude.)
    clear order
    order(1).m2=0;
    order(1).m1=-m_max:m_max;

    % Run the diffraction calculations.
    tic
        [param_size,scat_field,inc_field]=gdc(grating,inc_field,order,false);
    toc

    % Compute the diffraction efficiencies.
    [R,T]=gdc_eff(scat_field,inc_field);
    % Discard diffracted waves that decay exponentially with distance from the
    % grating. (These include evanescent waves and, if the substrate's
    % permittivity is not real-valued, all transmitted waves.)
    % (Note: "[scat_field.f1r]" is MATLAB syntax for
    % "scat_field(1).f1r, scat_field(2).f1r, ...".)
    R=R(imag([scat_field.f1r])==0);
    T=T(imag([scat_field.f1t])==0);
    % Extract the diffraction order indices for the reflected and transmitted
    % waves. (Only the m1 indices matter; the m2's are all zero.)
    R_m1=[R.m1].';
    T_m1=[T.m1].';
    % Extract the reflection efficiencies (R1...R4) and transmission
    % efficiencies (T1...T4) corresponding to the four incident polarization
    % states defined in gdc_eff.m.
    R1=[R.eff1].';
    R2=[R.eff2].';
    R3=[R.eff3].';
    R4=[R.eff4].';
    T1=[T.eff1].';
    T2=[T.eff2].';
    T3=[T.eff3].';
    T4=[T.eff4].';
    % Tabulate the diffraction order indices and diffraction efficiencies. Also
    % tabulate the fractional energy loss in the grating.
%     disp(' ');
%     disp('Diffraction efficiencies (m1, eff1, eff2, eff3, eff4)');
%     disp('R:');
%     disp(num2str([R_m1 R1 R2 R3 R4]));
%     disp('T:');
%     disp(num2str([T_m1 T1 T2 T3 T4]));
%     disp('Energy loss:');
%     disp(num2str(1-sum([[R1 R2 R3 R4]; [T1 T2 T3 T4]])));

    % Compute the diffraction efficiencies for the incident E field's
    % polarization normal to the grating lines (x3 axis). Assume an incident E
    % field complex amplitude of the form A*s+B*p, as described in gdc_eff.m,
    % with A and B defined so that the field is orthogonal to the x3 axis.
    A=-inc_field.p3;
    B=inc_field.s3;
    % Normalize abs(A)^2+abs(B)^2 to 1.
    C=1/sqrt(abs(A)^2+abs(B)^2);
    A=A*C;
    B=B*C;
    % Calculate the reflection and transmission efficiencies, denoted R5 and
    % T5, for the above-defined incident field.
    R5=abs(A)^2*R1+abs(B)^2*R2...
        +real(conj(A)*B)*(2*R3-R1-R2)-imag(conj(A)*B)*(2*R4-R1-R2);
    T5=abs(A)^2*T1+abs(B)^2*T2...
        +real(conj(A)*B)*(2*T3-T1-T2)-imag(conj(A)*B)*(2*T4-T1-T2);
    % Do the same for the orthogonal incident polarization (H field orthogonal
    % to the x3 axis); denote the corresponding efficiencies as R6, T6:
    [A,B]=deal(B,-A);
    R6=abs(A)^2*R1+abs(B)^2*R2...
        +real(conj(A)*B)*(2*R3-R1-R2)-imag(conj(A)*B)*(2*R4-R1-R2);
    T6=abs(A)^2*T1+abs(B)^2*T2...
        +real(conj(A)*B)*(2*T3-T1-T2)-imag(conj(A)*B)*(2*T4-T1-T2);
    % Tabulate the results.
%     disp(' ');
%     disp('Diffraction efficiencies (with H3=0, E3=0)');
%     disp('R:');
%     disp(num2str([R_m1 R6 R5]));
%     disp('T:');
%     disp(num2str([T_m1 T6 T5]));
%     disp('Energy loss:');
%     disp(num2str(1-sum([[R6 R5]; [T6 T5]])));
    
    
    output = R;
end
