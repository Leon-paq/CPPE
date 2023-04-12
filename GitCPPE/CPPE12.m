function [Gbest,Fit_min,fit_PEO]= CPPE9(fhd,Dimension,Particle_Number,Max_Gen,VRmin,VRmax,varargin)
%fhd数据,Dimension粒子维度,Particle_Number粒子数,Max_Gen迭代次数,VRmin最小值,VRmax最大值,
% fhd=str2func('cec14_func'); Fnum = 30;

% func_num =30;

% Max_Gen = Max_Gen;
Np = Particle_Number;
Lb = VRmin;
Ub = VRmax;
Dim = Dimension;

sigma = 0.1*(Ub-Lb);   
TT = sigma;
% 初始化k个最优解
K = floor(log(Np))+1; %floor朝负无穷大方向取整
L_pha = zeros(K,Dim); %zeros创建一个全为零元素的数组
L_Pha_fit = zeros(1,K);


% 初始化生成x
Pha=initialization(Np,Dim,Ub,Lb,12);
% 计算fitness
fitness = feval(fhd,Pha',varargin{:});
% k个最优解赋值
[~,f_index] = sort(fitness);
% 初始化全局最优
Fit_min = fitness(f_index(1));
Gbest = Pha(f_index(1),:);
% 初始化K个局部最优
for i = 1:K
    L_pha(i,:) = Pha(f_index(i),:);
    L_Pha_fit(i) =  fitness(f_index(i));
end
% 初始化趋势种群趋势
A = zeros(Np,Dim);
% 初始化种群数量x 
Px= (1/Np)*ones(1,Np);
% 初始化种群增长率
Pa = 1.1*ones(1,Np);

count = 1;
fit_PEO(count) = Fit_min;
% pic_num = 1;
for t = 2:Max_Gen
    count = count + 1;
    % 计算新的位置
    new_Pha = Pha + A;
    %边界约束
    new_Pha=space_bound(new_Pha,Ub,Lb);
    
    % 计算所有的fitness
    new_fitness = feval(fhd,new_Pha',varargin{:});
    new_fit_best = min(new_fitness); 
    new_fit_worst = max(new_fitness);
    % 更新全局最优
    [~,new_f_index] = sort(new_fitness);
    if new_fitness(new_f_index(1)) <= Fit_min
        Fit_min = new_fitness(new_f_index(1));
        Gbest = new_Pha(new_f_index(1),:);
    end
    % 更新局部最优
    for j = 1:K
        if new_fitness(new_f_index(j))<= L_Pha_fit(j)
            L_Pha_fit(j) = new_fitness(new_f_index(j));
            L_pha(j,:) = new_Pha(new_f_index(j),:);
        end
    end
    % 更新解的位置
    for p = 1:Np
        % 位置是否更新
        if new_fitness(p) <= fitness(p)
            Pha(p,:) = new_Pha(p,:);
            fitness(p) = new_fitness(p);
            Px(p) = Pa(p)*Px(p)*(1-Px(p));
            A1 = choosebest(L_pha,Pha(p,:)); % 趋同进化，最近最优
            A1=A1.*0.2;
            A3 = tubian(Dim,sigma);
            A(p,:) =  (1-Px(p)).*A1+ Px(p).*(A(p,:) + A3);
        else
            if rand <= (Px(p))
                Pha(p,:) = new_Pha(p,:);
                fitness(p) = new_fitness(p);
                Px(p) = Pa(p)*Px(p)*(1-Px(p));
            end
            A1 = choosebest(L_pha,Pha(p,:)); 
            A(p,:) = rand(1,Dim).*A1+TT.*randn(1,Dim); 
            TT = TT*0.99;
        end
        
        %计算种群间的竞争与共存
        temp_hab = sigma; 
        temp = randperm(Np);
        if temp(1) ~=p, tp = temp(1);
        else 
            tp = temp(2);
        end
        % 是否会出现竞争
        if dist(Pha(p,:),Pha(tp,:)') < temp_hab*((Max_Gen+1-t)/Max_Gen)
            d_p = Pa(p)*Px(p)*(1-Px(p)-(fitness(tp)/fitness(p))*Px(tp));
            Px(p) = Px(p) + d_p;
            A(p,:) = A(p,:) + ((fitness(tp)-fitness(p))/fitness(tp)).*(Pha(tp,:)-Pha(p,:));
            
        end

        
        % 约束a的取值
        if Pa(p) <=0.1 || Pa(p) >=4 || Px(p) <=0.001
            Px(p) = (1/Np);
            Pa(p) = 1.1;
            Pha(p,:) = Lb.*ones(1,Dim) + (Ub-Lb).*ones(1,Dim).*rand(1,Dim);
            A(p,:) = zeros(1,Dim);
            fitness(p) = Inf;
        end
    end
    fit_PEO(count) = Fit_min;

end          
end
function A3 = tubian(Dim,sigma)
S = floor(rand/(1/Dim));
A3 = zeros(1,Dim);
if S >=1 && S<=Dim
    J = randperm(Dim,S);
    A3(J) = 1;
    A3 = A3.*2.*sigma.*randn(1,Dim);
%     A3 = (S/Dim).*A3;
end

end

function A1 = choosebest(L_pha,newpha)
[K, dim] = size(L_pha);
temp_dsit = zeros(1,K);
for i = 1:K
    temp_dsit(i) = dist(L_pha(i,:),newpha');
end
[~,index]=min(temp_dsit);
% 第一种 
A1 = (L_pha(index,:)- newpha);
% 第二种 
end
% This function initialize the first population of search agents
% function Positions=initialization(SearchAgents_no,dim,ub,lb)
% 
% Boundary_no= size(ub,2); % numnber of boundaries
% if Boundary_no==1
%     Positions=rand(SearchAgents_no,dim).*(ub-lb)+lb;
% end
% 
% if Boundary_no>1
%     for i=1:dim
%         ub_i=ub(i);
%         lb_i=lb(i);
%         Positions(:,i)=rand(SearchAgents_no,1).*(ub_i-lb_i)+lb_i;
%     end
% end
% end
%This function checks the search space boundaries for agents.
function  X=space_bound(X,up,low)

[N,dim]=size(X);
for i=1:N 
    Tp=X(i,:)>up;Tm=X(i,:)<low;X(i,:)=(X(i,:).*(~(Tp+Tm)))+((rand(1,dim).*(up-low)+low).*(Tp+Tm));

end
end
