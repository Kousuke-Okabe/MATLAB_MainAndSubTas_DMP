function[G,V] = get_DMPgraph(type,x)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   DMP�̒��_�O���t��Return����֐�
%
%   [G,V] = get_DMPgraph(type, x)
%       G : �O���t
%       V : ���_�f�[�^
%       x : ��Ƌ�ԏ�̈ʒu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

draw = 0;

% %%
% clear all
% 
% x = [ 0.3;
%       0.2;
%       0.15;
%       0.1];
% draw = 1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% �A�[���^�C�v
type1 = '2d_RRRR';
type2 = '2d_RR';


% �ϐ��錾
param = get_parameter(type1);
    Link = param.Link;
    Tlim = param.Tlim;

% ���_�f�[�^�쐬
V = nan(Link,2^Link);
for i = 1:Link
    for j = 1:size(V,2)
        if j-floor((j-1)/(2^(Link+1)/2^i))*(2^(Link+1)/2^i) <= 2^Link/2^i
            V(i,j) = Tlim(i);
        else
            V(i,j) = -Tlim(i);
        end
    end
end

% �O���t�쐬
G = zeros(2^Link);    % �O���t

for i = 1:size(G,2)
    for j = i+1:size(G,2)
        temp = 0;
        for k = 1:Link
            if V(k,i)*V(k,j) < 0
                temp = temp+1;
            end
        end
        
        if temp == 1
            G(i,j) = 1;
            G(j,i) = 1;
        end
    end
end

% MFP���_�쐬
q = fIKinematics(type1,x);
[M,g] = get_matrix_minus(type1, q);

for i = 1:size(V,2)
    V(:,i) = V(:,i)-g;
end

[J1,Jp1] = fJacobi_minus_q(type1,q);
[J2,~] = fJacobi_minus_q(type2,q(1:2));
J2 = [J2,zeros(2)];
Jp2 = J2'/(J2*J2');
Ji = [Jp1, (eye(Link)-Jp1*J1)*Jp2];
V = Ji\M\V;

% [J,Jp,U] = fJacobi(type, r,Rx);
% Je = [J;U'];
% V = Je/M*V;

if draw == 1
    FH = 1;
    figure(FH)
    clf(FH)
    fDraw_Graph3(FH,G,V(1:2,:),2)
end
