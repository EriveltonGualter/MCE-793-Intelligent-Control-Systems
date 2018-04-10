% Erivelton Gualter
% 26/03/2018

clc
clear all
close all

hold on;

rec = [5 5 7 7; 4 6 6 4];

rec1 = [-1 -1 1 1; -1 1 1 -1] + [6*ones(1,4); 5*ones(1,4);];
rec2 = [-1 -1 1 1; -1 1 1 -1] + [8.5*ones(1,4); 1.5*ones(1,4);];
rec3 = [-1 -1 1 1; -1 1 1 -1] + [10.5*ones(1,4); 2*ones(1,4);];

obstacle1 =  rectangle(...
        'Position',[rec1(1,1), rec1(2,1), rec1(1,3)-rec1(1,1), rec1(2,3)-rec1(2,1)],...
        'LineWidth',2,...
        'FaceColor',[1, 0, 0],...
        'EdgeColor',0.5*[1, 0, 0]); %Darker version of color

obstacle2 =  rectangle(...
        'Position',[rec2(1,1), rec2(2,1), rec2(1,3)-rec2(1,1), rec2(2,3)-rec2(2,1)],...
        'LineWidth',2,...
        'FaceColor',[1, 0, 0],...
        'EdgeColor',0.5*[1, 0, 0]); %Darker version of color
    
obstacle3 =  rectangle(...
        'Position',[rec3(1,1), rec3(2,1), rec3(1,3)-rec3(1,1), rec3(2,3)-rec3(2,1)],...
        'LineWidth',2,...
        'FaceColor',[1, 0, 0],...
        'EdgeColor',0.5*[1, 0, 0]); %Darker version of color
   
axis([0 12 0 10]);
axis equal
drawnow


q0 = [4.5 6.5];
qf = [10.5 0.5];

plot([q0(1), qf(1)], [q0(2), qf(2)],'ok','LineWidth',2);

three.id = 1;
three.x = q0(1);
three.y = q0(2);
three.parent = [];

threeplot = q0;

plot(threeplot(:,1), threeplot(:,2),'-o','LineWidth',2,...
                       'MarkerEdgeColor','b',...
                       'MarkerFaceColor','r',...
                       'MarkerSize',4)
%
k = 2;
Nmax = Inf;
while (k<Inf)
    
    three(k).id = k;
    three(k).x = 12*rand(1);
    three(k).y = 10*rand(1);
    
%     three(k).x = 9;
%     three(k).y = 0;
    
    [mIdx,mD] = knnsearch([three(1:end-1).x; three(1:end-1).y]', [three(k).x three(k).y]);
    three(k).parent = mIdx;
    
    idx_qplot = find([threeplot(:)] == [three(mIdx).x]);

    xl = linspace(three(mIdx).x, [three(k).x], 500);
    yl = linspace(three(mIdx).y, [three(k).y], 500);
    inl1 = inpolygon(xl, yl, rec1(1,:),rec1(2,:));
    inl2 = inpolygon(xl, yl, rec2(1,:),rec2(2,:));
    inl3 = inpolygon(xl, yl, rec3(1,:),rec3(2,:));

    inl = inl1 + inl2 + inl3;
    if sum(inl) > 1
        idx_l = find(inl,1)-1;
        if idx_l < 1
            flag = 0;
        else
            three(k).x = xl(idx_l);
            three(k).y = yl(idx_l);
        end
    end    
        
    if flag
        threeplot = [threeplot(1:idx_qplot,:); [three(k).x three(k).y]; threeplot(idx_qplot:end,:)];
        plot(threeplot(:,1), threeplot(:,2),'-b');

        qk = [three(k).x three(k).y];

        if norm(qk-qf) < 0.1
            break;
        end
        k = k + 1;
    else
        flag = 1;
    end
    pause(1/15)
end
%%
pos = k;
path(1,1) = three(pos).x;
path(1,2) = three(pos).y;
path(1,3) = pos;
pos = three(pos).parent;
for i=2:k
    path(i,3) = pos;
    path(i,1) = three(pos).x;
    path(i,2) = three(pos).y;
    pos = three(pos).parent;
    
    if pos == 1
        path(i+1,1) = three(1).x;
        path(i+1,2) = three(1).y;
        break;
    end
end
plot(path(:,1), path(:,2), '-r', 'LineWidth',2);

vertices = k;
edges = length(path);
display(['vertices = ',num2str(vertices)]);
display(['edges = ',num2str(edges)]);

axis([0 12 0 10]);
axis equal
