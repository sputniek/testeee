tic
clear all
close all
%% INPUTS

%% testando 1 2 3

% ASA - [cma, b, S_w, cdf, cl_w0, cd_w0, cm_w0, cl_w2, cd_w2]
% EH - [S_eh, d_CA, cl0_eh, cd0_eh, cl0def_eh, cd0def_eh, cl2def_eh, cd2def_eh, deflexao]
% GMP - [at, bt, ct]
% GEO - [x_tp, dist_principal, d_Wz, h_boom]
% FT - [FCL, FCD, FAT, FBT, mi, Iyy, lim, alt]
[~,~,~,rho] = atmosisa(1250);
% ASA = [0.27, 2.7, 0.729, 0.0219, 0.88, 0.042, -0.286, 1.049, 0.054];
% EH = [0.112, 0.85, -0.475, 0.047, -0.985, 0.158, -0.883, 0.138, -15.0];
% GMP=[-0.019639512*(1.225/rho), -0.178909728*(1.225/rho),  1.996*9.80665];
% GEO = [0.5, 0.18827, 0.11973, 0];
% FT = [0.98, 1.37, 0.9091946, 0.60196, 0.065, 0.27230731605, 50, 1250.0];

ASA = [0.265, 2.48, 0.6572, 0.02597, 1.168, 0.0445, -0.327, 1.331, 0.05162];
EH15 = [0.12, 0.78, -0.45, 0.041, -1.052, 0.166, -0.944, 0.146, -15.0];
EH16 = [0.12, 0.78, -0.45, 0.041, -1.092, 0.177, -0.985, 0.155, -16.0];
EH17 = [0.12, 0.78, -0.45, 0.041, -1.132, 0.190, -1.025, 0.168, -17.0];
EH18 = [0.12, 0.78, -0.45, 0.041, -1.172, 0.202, -1.065, 0.180, -18.0];
EH19 = [0.12, 0.78, -0.45, 0.041, -1.212, 0.215, -1.105, 0.192, -19.0];
EH20 = [0.12, 0.78, -0.45, 0.041, -1.252, 0.228, -1.146, 0.204, -20.0];
EH21 = [0.12, 0.78, -0.45, 0.041, -1.291, 0.241, -1.186, 0.216, -21.0];

cd0 = importdata('C:\Users\Lusca\Desktop\MATLAB\raz_descida_relatorio\cd_viscoso_0_graus.txt');
cd0 = cd0.data();
cd2 = importdata('C:\Users\Lusca\Desktop\MATLAB\raz_descida_relatorio\cd_viscoso_2_graus.txt');
cd2 = cd2.data();
pf0 = polyfit(cd0(:,1),cd0(:,2),5);
pf2 = polyfit(cd2(:,1),cd2(:,2),4);

h17x8N = [-0.02993182338323175, -0.5657868430012268, 74.34598978231658];

GEO = [0.5, 0.21761, 0.113, 0.22];
% FT = [0.98, 1.37,  0.9091946, 0.60196,0.065,0.44387282335, 50, 1250.0];
FT = [0.98, 1.37,  1, 1,0.065,0.44387282335, 50, 1250.0];
FT = [0.98, 1.37,  1, 1,0.065,0.65866458948, 50, 1250.0];
% a = decola(14.77,36,ASA,EH18,h17x8N,GEO,FT,pf0,pf2)
h17x8N = [-1.39401495e-02 -9.12258133e-01  7.61008921e+01];
%{
[m15,~,a15,~] = MTOW(ASA, EH15, h17x8N, GEO, FT,pf0,pf2);
[m16,~,a16,~] = MTOW(ASA, EH16, h17x8N, GEO, FT,pf0,pf2);
[m17,~,a17,~] = MTOW(ASA, EH17, h17x8N, GEO, FT,pf0,pf2);
[m18,~,a18,~] = MTOW(ASA, EH18, h17x8N, GEO, FT,pf0,pf2);
[m19,~,a19,~] = MTOW(ASA, EH19, h17x8N, GEO, FT,pf0,pf2);
[m20,~,a20,~] = MTOW(ASA, EH20, h17x8N, GEO, FT,pf0,pf2);
[m21,~,a21,~] = MTOW(ASA, EH21, h17x8N, GEO, FT,pf0,pf2);

%% Grafico
yyaxis left
title('MTOW e Ângulo de Decolagem x Deflexão do Profundor','FontSize',16,'FontName','Times','FontWeight','bold')
bar([15,16,17,18,19],[m15 0;m16 0;m17 0;m18 0;m19 0],'FaceColor',[0 0 1])
hold on
xlabel('Deflexão do Profundor','FontSize',14,'FontName','Times','FontWeight','bold');
ylabel('MTOW','FontSize',14,'FontName','Calibri','FontWeight','bold');
ylim([13.5 15])
yticks([13.5:0.25:15])

yyaxis right
bar([15,16,17,18,19],[0 a15;0 a16;0 a17;0 a18;0 a19],'FaceColor',[1 0 0])
ylabel('Ângulo de Decolagem','FontSize',14,'FontName','Times','FontWeight','bold')
grid
%}
% ylim([3 6])
% legend('MTOW','Angulo de Decolagem')
% hold on
% bar([15,16,17,18,19],[a15,a16,a17,a18,a19],'b')


[x1,x2,x3,x4] = MTOW(ASA, EH18, h17x8N, GEO, FT,pf0,pf2);
fprintf("\nMTOW: %.2f", x1)
fprintf("\nProfundor: %.2f", x2)
fprintf("\nAngulo: %.2f", x3)
fprintf("\nVelocidade: %.2f\n", x4)
toc


## mudanca de teste

xx = 33;

function deg = rad2deg(num)
    deg = num*180/pi;
end
    
%% CALCULA COMPRIMENTO DE PISTA
function [metros,angulo,vel] = decola(m,p_acionamento,ASA,EH,GMP,GEO,FT,CD0,CD2)
    %{
    Inputs:
        m - Massa [kg]
        p_acionamento - Posicao de acionamento do profundor [m]
        ASA - [cma, b, S_w, cdf, cl0_w, cd0_w, cm0_w, cl2_w, cd2_w]
        EH - [S_eh, d_CA, cl0_eh, cd0_eh, cl0def_eh, cd0def_eh, cl2def_eh, 
          cd2def_eh, deflexao]
        GMP - [at, bt, ct]
        GEO - [x_tp, dist_principal, d_Wz, h_boom]
        FT - [FCL, FCD, FAT, FBT, mi, Iyy, lim, alt]
        
    Outputs:
        [metros,angulo,vel]
        metros - Comprimento de pista utilizado [m]
        angulo - Angulo de decolagem [graus]
        vel - Velocidade de decolagem [m/s]
    %}
    
%% PARaMETROS FIXOS INICIAIS
    alt_densidade = FT(8);
    [~, ~, ~, rho] = atmosisa(alt_densidade);

    FCL = FT(1);
    FCD = FT(2);
    FAT = FT(3);
    FBT = FT(4);
    g = 9.80665;                                                                % Gravidade
    u = 0;                                                                      % Velocidade do vento

%% DADOS DA ASA
    cma = ASA(1);                                                                % Corda media aerodinamica
    % b = ASA[2] 
    S_w = ASA(3);                                                              % angulo de estol

    CD_f = ASA(4);
    CL0_w = ASA(5)*FCL;                                                          % CL para alpha 0 da asa 
    CD0_w = ASA(6)*FCD + CD_f;                                                  % CD para alpha 0 da asa
    CM0_w = ASA(7);                                                             % Cm para alpha 0 da asa
    CDI_w = 0.03584*FCD;
    CDI_w2 = 0.04516*FCD;
    CDI_alpha_w = (CDI_w2-CDI_w)/2;
    
    CL0_w2 = ASA(8)*FCL;                                                         
    CD0_w2 = ASA(9)*FCD + CD_f; 

    CL_alpha_w = (CL0_w2-CL0_w)/2;
    CD_alpha_w = (CD0_w2-CD0_w)/2;

%% DADOS EMPENAGEM HORIZONTAL
    S_ht = EH(1);  
    % cma_ht = 0.2
    d_CA = EH(2); %Distancia entre CA's
    
    CL0_ht = EH(3)*FCL;
    CD0_ht = EH(4)*FCD;

    deflexao = EH(9);

    CLmax_ht = EH(5)*FCL;
    CDmax_ht = EH(6)*FCD;

    CLtrim_ht = EH(7)*FCL;
    CDtrim_ht = EH(8)*FCD;

    CL_alpha_ht = (CLtrim_ht - CLmax_ht)/(2);
    CD_alpha_ht = abs((CDtrim_ht - CDmax_ht)/(2));

%% DADOS GMP
    % Coeficientes da curva de tracao do GMP (Polinomio 2 grau) MaXIMA ROTAcaO

    at = GMP(1)*FAT*(rho/1.225);
    bt = GMP(2)*FBT*(rho/1.225);
    ct = GMP(3)*(rho/1.225);
    GMP = [at,bt,ct];

%% DADOS GEOMETRICOS
    x_tp = GEO(1)*cma; %Bordo de Ataque ao Trem de Pouso
    d_CG = 0.25*cma; %Bordo de Ataque ao CG
    %D_helice = 18
    
    %dist_principal = 188.27E-3
    dist_principal = GEO(2);
    d_Lw = x_tp - 0.25*cma;
    d_Lht = -(d_CA - x_tp + 0.25*cma);
    d_Wx = -(x_tp - d_CG);
    %d_Wz = 119.73E-3
    d_Wz = GEO(3);
    d_T = -(dist_principal);
    d_Dw = dist_principal; 
    d_Dht = dist_principal+GEO(4);
    %rRoda = 2.5E-3
    %d_h = dist_principal + rRoda
    
    mi = FT(5);                                                          % Coeficiente de atrito dinamico
    Iyy_CG = FT(6);                 %Momento de inercia
    t_acionamento = 0.5;
    
    var_tempo = 0.001;
    rot = 0;

    W = -m*g;
    Iyy = Iyy_CG + m*(d_Wx^2 + d_Wz^2);
    
%% INICIA VARIAVEIS
    cont = 1;
    t = [0];
    s = [0];
    v = [0];
    vx = [0];
    defl = [0];
    a_alpha = [0];
    v_alpha = [0];
    alpha = [0];
    alphar = [0];
    T = [polyval(GMP,u)];
    Tx = [T(cont)];
    Tz = [0];
    CL_w = CL0_w;
    CD_w = polyval(CD0,v(cont));
%     CD_w = CD0_w;
    CL_ht = CL0_ht;
    CD_ht = CD0_ht;
    D_w = [0.5*rho*(u)^2*CD_w*S_w];
    D_ht = [0.5*rho*(u)^2*CD_ht*S_ht];
    D = [D_w(cont) + D_ht(cont)];

    L_w = [0.5*rho*(u)^2*CL_w*S_w];
    L_ht = [0.5*rho*(u)^2*CL_ht*S_ht];
    L = [L_w(cont) + L_ht(cont)];

    N = [-L(cont)-W-Tz(cont) ];
    Fat = [mi*N(cont) ];
    a = [(T(cont) -D(cont) -Fat(cont) )/m];
    
    M_CM_w = [ 0.5*(u)^2*CM0_w*S_w*cma ];
    M_Lw = [ L_w(cont)*d_Lw ];
    M_Lht = [ L_ht(cont)*d_Lht ];
    M_Dw = [ D_w(cont)*d_Dw ];
    M_Dht = [ D_ht(cont)*d_Dht ];
    M_T = [ T(cont)*d_T ];
    M_W = [ -W*d_Wx ];
    M_a = [ m*a(cont)*(d_Wz) ];
    M_res = [ M_CM_w(cont) + M_Lw(cont) + M_Lht(cont) + M_Dw(cont)+ M_Dht(cont) + M_T(cont) + M_W(cont) + M_a(cont)];

%% ITERACOES DE TEMPO
    while 1
        
        cont = cont + 1;
        t(cont) = (t(cont-1)+var_tempo);
        v(cont) = (v(cont-1) + a(cont-1)*var_tempo);
        
        T(cont) = (polyval(GMP,v(cont-1)+u));
        Tx(cont) = (T(cont)*cos(alphar(cont-1)));
        Tz(cont) = (T(cont)*sin(alphar(cont-1)));

        CL_w = CL0_w + CL_alpha_w*alpha(cont-1);
        CD0_w = polyval(CD0,v(cont))*FCD;
        CD0_w2 = polyval(CD2,v(cont))*FCD;
        CD_alpha_w = (CD0_w2-CD0_w)/2;
        CD_w = CDI_w + CD0_w + CDI_alpha_w*alpha(cont-1) + CD_alpha_w*alpha(cont-1) + CD_f;
%         CD_w = CD0_w + CD_alpha_w*alpha(cont-1);

        if s(cont-1) < p_acionamento
            CL_ht = CL0_ht;
            CD_ht = CD0_ht;
            defl(cont) = (0);
        elseif s(cont-1) >= p_acionamento && defl(cont-1) < abs(deflexao)
            CL_def = (CLmax_ht - CL0_ht)/(abs(deflexao));
            CD_def = (CDmax_ht - CD0_ht)/(abs(deflexao));
            defl(cont) = (defl(cont-1) + var_tempo*abs(deflexao)/t_acionamento);
            CL_ht = CL0_ht + CL_def*defl(cont);
            CD_ht = CD0_ht + CD_def*defl(cont);
        else
            CL_ht = CLmax_ht + CL_alpha_ht*alpha(cont-1);
            CD_ht = CDmax_ht + CD_alpha_ht*alpha(cont-1);
            defl(cont) = (defl(cont-1));
        end
                

        D_w(cont) = (0.5*rho*(v(cont-1)+u)^2*CD_w*S_w);
        D_ht(cont) = (0.5*rho*(v(cont-1)+u)^2*CD_ht*S_ht);
        D(cont) = (D_w(cont) + D_ht(cont));

        L_w(cont) = (0.5*rho*(v(cont-1)+u)^2*CL_w*S_w);
        L_ht(cont) = (0.5*rho*(v(cont-1)+u)^2*CL_ht*S_ht);
        L(cont) = (L_w(cont) + L_ht(cont));

        N(cont) = (-L(cont-1)-W-Tz(cont-1));
        Fat(cont) = (mi*N(cont-1));

        a(cont) = ((T(cont-1)-D(cont-1)-Fat(cont-1))/m);
        
        vx(cont) = (vx(cont-1) + a(cont-1)*var_tempo);
        s(cont) = (s(cont-1) + vx(cont-1)*var_tempo + 0.5*a(cont-1)*var_tempo^2);

        M_CM_w(cont) = (0.5*rho*(v(cont-1)+u)^2*CM0_w*S_w*cma);
        M_Lw(cont) = (L_w(cont-1)*d_Lw);
        M_Lht(cont) = (L_ht(cont-1)*d_Lht);
        M_Dw(cont) = (D_w(cont-1)*d_Dw);
        M_Dht(cont) = (D_ht(cont-1)*d_Dht);
        M_T(cont) = ( T(cont-1)*d_T);
        M_W(cont) = ( -W*d_Wx);
        M_a(cont) = ( m*a(cont-1)*(d_Wz));
        M_res(cont) = (M_CM_w(cont-1) + M_Lw(cont-1) + M_Lht(cont-1) + M_Dw(cont-1)+ M_Dht(cont-1) + M_T(cont-1) + M_W(cont-1) + M_a(cont-1));

        if M_res(cont-1)>0 | rot == 1
            a_alpha(cont) = (M_res(cont-1)/Iyy);
            v_alpha(cont) = (v_alpha(cont-1) + a_alpha(cont-1)*var_tempo);
            alphar(cont) = (alphar(cont-1) + v_alpha(cont-1)*var_tempo + 0.5*a_alpha(cont-1)*var_tempo^2);
            alpha(cont) = (rad2deg(alphar(cont)));
            rot = 1;
        else
            a_alpha(cont) = (0);
            v_alpha(cont) = (0);
            alpha(cont) = (0);
            alphar(cont) = (0);
        end
        

        if N(cont) < 0
            break
        elseif s(cont) > FT[7]+10
            s(cont) = (100);
            break
        end
    end
        
    metros = s(cont);
    angulo = alpha(cont);
    vel = v(cont);
end
       

%% DETERMINA MTOW
function [TOW,Profundor,Angulo,Velocidade] = MTOW(ASA,EH,GMP,GEO,FT,CD0,CD2)
    %{
    Inputs:
        ASA - [cma, b, S_w, cdf, cl0_w, cd0_w, cm0_w, cl2_w, cd2_w]
        EH - [S_eh, d_CA, cl0_eh, cd0_eh, cl0def_eh, cd0def_eh, cl2def_eh, 
          cd2def_eh, deflexao]
        GMP - [at, bt, ct]
        GEO - [x_tp, dist_principal, d_Wz, h_boom]
        FT - [FCL, FCD, FAT, FBT,mi, Iyy, lim, alt]
        
    Outputs:
        [MTOW,profundor,angulo,vel]
        MTOW - Maximum TakeOff Weight [kg]
        profundor - Posicao de acionamento do profundor [m]
        angulo - Angulo de decolagem [graus]
        vel - Velocidade de decolagem [m/s]
    %}
    
    final = [];
    res = 0.001;
    lim = FT(7);
    cont1 = 1;
    for aciona = [lim-20:lim-10]
        mini = 1.;
        maxi = 30.;
        while (maxi-mini)>res
            m = (maxi+mini)/2;
            [comp,ang,vel] = decola(m,aciona,ASA,EH,GMP,GEO,FT,CD0,CD2);
            if comp > lim
                maxi = m;
            else
                mini = m;
            end
        end
        final = [final;aciona,mini,ang,vel];
        cont1 = cont1 + 1;
    end
      
    TOW = max(final(:,2));
    [x,~] = find(final==TOW);
    Profundor = final(x,1);
    Angulo = final(x,3);
    Velocidade = final(x,4);

end
