function X = initialization(SearchAgents_no,dim,ub,lb,type)
switch type
    case 0
        Boundary_no = size(ub,2); % numnber of boundaries
        if Boundary_no == 1
            X=rand(SearchAgents_no,dim).*(ub-lb)+lb;
        end

        if Boundary_no>1
            for i=1:dim
                ub_i=ub(i);
                lb_i=lb(i);
                X(:,i)=rand(SearchAgents_no,1).*(ub_i-lb_i)+lb_i;
            end
        end
    case 1   
        disp('Logistic 混沌映射')
        Z = zeros(SearchAgents_no, dim);
        % 随机生成一个d维向量
        Z(1, :) = rand(1, dim);
        % 利用logistic生成n_pop个向量
        for i=2:SearchAgents_no
            Z(i,:) = 4.0*Z(i-1,:).*(1-Z(i-1,:));
        end
        % 将z的各个分量载波到对应变量的取值区间
        X = zeros(SearchAgents_no, dim);
        for i=1:SearchAgents_no
            X(i,:) = lb + (ub - lb)*Z(i,:);
        end
    case 2
        disp('Piecewise 混沌映射')
        Z = rand(SearchAgents_no, dim);
        p = 0.4;
        for i=1:SearchAgents_no
            for j=1:dim
                if Z(i,j)>=0 & Z(i,j)<p
                    Z(i,j) = Z(i,j)/p;
                elseif Z(i,j)>=p & Z(i,j)<0.5
                    Z(i,j) = (Z(i,j) - p)/(0.5-p);
                elseif Z(i,j)>=0.5 & Z(i,j)<1-p
                    Z(i,j) = (1 - p - Z(i,j))/(0.5-p);
                elseif Z(i,j)>=1-p & Z(i,j)<1
                    Z(i,j) = (1 - Z(i,j))/p;
                end
            end
        end
        X = lb + Z*(ub - lb);
    case 3
        disp('Singer 混沌映射')
        Z = rand(SearchAgents_no, dim);
        u = 1.07;
        for i=1:SearchAgents_no
            for j=1:dim
                Z(i,j) = u*(7.86*Z(i,j)-23.31*Z(i,j).^2+28.75*Z(i,j).^3-13.302875*Z(i,j).^4);
            end
        end
        X = lb + Z*(ub - lb);
    case 4
        disp('Sine 混沌映射')
        Z = zeros(SearchAgents_no, dim);
        % 随机生成一个d维向量
        Z(1, :) = rand(1, dim);
        a = 1;
        for i=2:SearchAgents_no
            Z(i, :) = a.*sin(pi.*Z(i-1,:));
        end
        X = zeros(SearchAgents_no, dim);
        for i=1:SearchAgents_no
            X(i,:) = lb + (ub - lb)*Z(i,:);
        end
    case 5
        disp('Gauss 混沌映射')
        Z = rand(SearchAgents_no, dim);
        mu =1;
        for i=1:SearchAgents_no
            for j=1:dim
                if Z(i,j) == 0
                    Z(i,j)=0;
                else
                    Z(i,j)=mod(mu\Z(i,j),1);
                end
            end
        end
        X = lb + Z*(ub - lb);
    case 6  
        disp('Tent 混沌映射')
        Z = rand(SearchAgents_no, dim);
        alpha = 0.4;
        for i=1:SearchAgents_no
            for j=1:dim
                if Z(i,j)<=alpha
                    Z(i,j)=Z(i,j)/alpha;
                elseif Z(i,j)>alpha
                    Z(i,j)=(1- Z(i,j))/(1 - alpha);
                end
            end
        end
        X = lb + Z*(ub - lb);
    case 7
        disp('Bernoulli 混沌映射')
        Z = rand(SearchAgents_no, dim);
        lamda = 0.4;
        for i=1:SearchAgents_no
            for j=1:dim
                if Z(i,j)<=(1-lamda) & Z(i,j)>0
                    Z(i,j)=Z(i,j)/(1-lamda);
                else
                    Z(i,j)=(Z(i,j)-1+lamda)/lamda;
                end
                if Z(i,j)<0
                    Z(i,j) = -1*Z(i,j);
                end
            end
        end
        X = lb + Z*(ub - lb);
    case 8
        disp('Chebyshev 混沌映射')
        Z = rand(SearchAgents_no, dim);
        fai = 4;
        for i=1:SearchAgents_no
            for j=1:dim
                Z(i,j)= cos(fai.*acos(Z(i,j)));
                if Z(i,j)<0
                    Z(i,j) = -1*Z(i,j);
                end
            end
        end
        X = lb + Z*(ub - lb);
    case 9
        disp('Circle 混沌映射')
        Z = rand(SearchAgents_no, dim);
        a=0.5;
        b=2.2;
        for i=1:SearchAgents_no
            for j=1:dim
                Z(i,j)= Z(i,j) + a - mod(b/(2*pi)*(sin(2*pi*Z(i,j))),1);
                %Z(i,j)= mod(Z(i,j) + a - b/(2*pi)*(sin(2*pi*Z(i,j))),1);
            end
        end
        X = lb + Z*(ub - lb);
    case 10
        disp('Cubic 混沌映射')
        Z = rand(SearchAgents_no, dim);
        p=2.595;
        for i=1:SearchAgents_no
            for j=1:dim                
                Z(i,j) = p.*Z(i,j).*(1-Z(i,j).^2);
            end
        end
        X = lb + Z*(ub - lb);
    case 11
        disp('Sinusoidal 混沌映射')
        Z = rand(SearchAgents_no, dim);
        a=2.3;
        for i=1:SearchAgents_no
            for j=1:dim                
                Z(i,j) = a*Z(i,j).^2*sin(pi*Z(i,j));
            end
        end
        X = lb + Z*(ub - lb);
    case 12
        disp('ICMIC 混沌映射')
        Z = rand(SearchAgents_no, dim);
        a=2;
        for i=1:SearchAgents_no
            for j=1:dim                
                Z(i,j) = sin(a/Z(i,j));
                if Z(i,j)<0
                    Z(i,j) = -1*Z(i,j);
                end
            end
        end
        X = lb + Z*(ub - lb);
end       
end