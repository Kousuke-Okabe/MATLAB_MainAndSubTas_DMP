%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   ***
%
%                                                       20.03.16 by.OKB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all

type = '2d_RRRR'

% set parameter
x = [0.2; 0.15; 0.1; 0.1];
dx = [0; 0; 0; 0];
g = [0;0;0];

 
% DMP計算
[Graph,V] = get_DMPgraph(type, x);
% Gr = 1/10000;
% Gz = 1/100000;
% Gx = eye(size(dx,1));
%     Gx(1,1) = Gr;
%     Gx(2,2) = Gr;
%     Gx(3,3) = Gz;
% 
% V = Gx*V;
% Vt = get_DMPtranlate(type, x,dx,g,f);
% Vt = Gx*Vt;
% 
% 
% for i = 1:size(V,2)
%     V(:,i) = V(:,i) + Vt;% + [x(1:2);0];
% end

% マニピュレータ描画
FH = 2;
figure(FH)
clf(FH)
q = fIKinematics(type,x);
fRoboAnimation(type,FH,q, x(1:2),0);

%%
% DMP描画
figure(FH)
% clf(FH2)
    xlabel('ddx [m/s^2]')
    ylabel('ddy [m/s^2]')
    zlabel('ddz [m/s^2]')
    
%     view([-10,42])
%     rotate3d on    


fDraw_Graph3(FH,Graph,V+x*ones(1,8),3);
quiver(x(1),x(2), Vt(1),Vt(2),'r', 'AutoScaleFactor',1)

xlim([0, 0.33])
ylim([-0.05,0.25])
box on
hold on

% %% 測定データ描画 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% data = csvread('test.csv');
% res_0 = Gx(1:2,1:2)*data(:,1:2)' + x(1:2)*ones(1,length(data));
% res_10 = Gx(1:2,1:2)*data(:,3:4)' + x(1:2)*ones(1,length(data)); 
% 
% figure(FH)
% hold on
% % plot(res_0(1,:),res_0(2,:),'b.')
% plot(res_10(1,:),res_10(2,:),'ro')